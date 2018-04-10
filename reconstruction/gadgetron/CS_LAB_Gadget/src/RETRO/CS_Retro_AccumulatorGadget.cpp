/*
file name	:	CS_Retro.cpp
author		:	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
version		:	v1.0
date		:	08.02.2015
description	:
references	:	 AccumulatorGadget.cpp (origial Gadgetron - mri_core project)
changes		:
*/

#include "CS_Retro_AccumulatorGadget.h"

using namespace Gadgetron;

// class constructor
CS_Retro_AccumulatorGadget::CS_Retro_AccumulatorGadget()
{
	GlobalVar::instance()->vPE_.clear();
	GlobalVar::instance()->vPA_.clear();
}

// class destructor - delete temporal buffer/memory
CS_Retro_AccumulatorGadget::~CS_Retro_AccumulatorGadget()
{
	// free bufferkSpace_
	while (bufferkSpace_.size() > 0) {
		if (bufferkSpace_.at(bufferkSpace_.size()-1)) {
			bufferkSpace_.at(bufferkSpace_.size()-1)->release();
		}

		bufferkSpace_.pop_back();
	}
}

// read flexible data header
int CS_Retro_AccumulatorGadget::process_config(ACE_Message_Block *mb)
{
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	ISMRMRD::IsmrmrdHeader h;
	ISMRMRD::deserialize(mb->rd_ptr(),h);

	ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
	ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
	ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

	// get matrix size
	dimensionsIn_[0] = r_space.matrixSize.x;
	dimensionsIn_[1] = e_space.matrixSize.y;
	dimensionsIn_[2] = e_space.matrixSize.z;

	// get FOV
	field_of_view_[0] = r_space.fieldOfView_mm.x;
	field_of_view_[1] = e_space.fieldOfView_mm.y;
	field_of_view_[2] = e_space.fieldOfView_mm.z;

	// get echo line and echo partition
	iEchoLine_ = e_limits.kspace_encoding_step_1.get().center;
	iEchoPartition_ = e_limits.kspace_encoding_step_2.get().center;

	// repetition time
	GlobalVar::instance()->fTR_ = h.sequenceParameters.get().TR.get().at(0);
#else
	// read xml header file
	boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

	// create sequence encoding parameters
	ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
	ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
	ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
	ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

	// get matrix size
	dimensionsIn_[0] = r_space.matrixSize().x();
	dimensionsIn_[1] = e_space.matrixSize().y();
	dimensionsIn_[2] = e_space.matrixSize().z();

	// get FOV
	field_of_view_[0] = r_space.fieldOfView_mm().x();
	field_of_view_[1] = e_space.fieldOfView_mm().y();
	field_of_view_[2] = e_space.fieldOfView_mm().z();

	// get echo line and echo partition
	iEchoLine_ = e_limits.kspace_encoding_step_1().get().center();
	iEchoPartition_ = e_limits.kspace_encoding_step_2().get().center();

	// repetition time
	GlobalVar::instance()->fTR_ = cfg->sequenceParameters().get().TR().at(0);
#endif

	// some debug output
	GDEBUG("Matrix size: %d, %d, %d\n", dimensionsIn_[0], dimensionsIn_[1], dimensionsIn_[2]);
	GDEBUG("FOV: %f, %f, %f\n", field_of_view_[0], field_of_view_[1], field_of_view_[2]);
	GDEBUG("echo line: %i, echo partition: %i\n", iEchoLine_, iEchoPartition_);
	GDEBUG("TR: %f\n", GlobalVar::instance()->fTR_);

	// set properties
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	GlobalVar::instance()->iNavPeriod_			= NavPeriod.value();
	GlobalVar::instance()->iNavPERes_			= NavPERes.value();
	GlobalVar::instance()->iMeasurementTime_	= MeasurementTime.value();
	iNPhases_									= Phases.value();
#else
	GlobalVar::instance()->iNavPeriod_			= *(get_int_value("NavPeriod").get());
	GlobalVar::instance()->iNavPERes_			= *(get_int_value("NavPERes").get());
	GlobalVar::instance()->iMeasurementTime_	= *(get_int_value("MeasurementTime").get());
	iNPhases_									= *(get_int_value("Phases").get());
#endif

	int iESPReSSoY = 0;
	int iESPReSSoZ = 0;
	iBodyRegion_ = 0;
	iVDMap_ = 0;
	iSamplingType_ = 0;
	fCSAcc_ = 0;
	fFullySa_ = 0;
	try {
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		if (h.encoding[0].trajectoryDescription) {
			GDEBUG("\n\nTrajectory description present!\n\n");
			ISMRMRD::TrajectoryDescription traj_desc = *h.encoding[0].trajectoryDescription;

			for (std::vector<ISMRMRD::UserParameterLong>::iterator i (traj_desc.userParameterLong.begin()); i != traj_desc.userParameterLong.end(); ++i) {
				if (i->name == "WIP_PF_y") {
					iESPReSSoY = i->value;
				} else if (i->name == "WIP_PF_z") {
					iESPReSSoZ = i->value;
				} else if (i->name == "SamplingType") {
					iSamplingType_ = i->value;
				} else if (i->name == "VDMap") {
					iVDMap_ = i->value;
				} else if (i->name == "BodyRegion") {
					iBodyRegion_ = i->value;
				} else if (i->name == "OuterIterations") {
					GlobalVar::instance()->iNOuter_ = i->value;
				} else if (i->name == "InnerIterations") {
					GlobalVar::instance()->iNInner_ = i->value;
				} else if (i->name == "fftSparseDim") {
					GlobalVar::instance()->iDimFFT_ = i->value;
				} else if (i->name == "dctSparseDim") {
					GlobalVar::instance()->iDimDCTSparse_ = i->value;
				} else if (i->name == "pcaSparseDim") {
					GlobalVar::instance()->iDimPCASparse_ = i->value;
				} else if (i->name == "kernelFftDim") {
					GlobalVar::instance()->iDimKernelFFT_ = i->value;
				} else if (i->name == "transformFftBaDim") {
					GlobalVar::instance()->iTransformFFTBA_ = i->value;
				} else if (i->name == "kSpaceOutDim") {
					GlobalVar::instance()->ikSpaceOut_ = i->value;
				} else if (i->name == "CSESPReSSo") {
					GlobalVar::instance()->bESPRActiveCS_ = i->value;
				} else if (i->name == "NavPeriod") {
					GlobalVar::instance()->iNavPeriod_ = i->value;
				} else if (i->name == "NavPERes") {
					GlobalVar::instance()->iNavPERes_ = i->value;
				} else if (i->name == "MeasurementTime") {
					GlobalVar::instance()->iMeasurementTime_ = i->value;
				} else if (i->name == "Phases") {
					iNPhases_ = i->value;
				} else if (i->name == "PopulationMode") {
					GlobalVar::instance()->iPopulationMode_ = i->value;
				} else if (i->name == "GatingMode") {
					GlobalVar::instance()->iGatingMode_ = i->value;
				}
			}

			for (std::vector<ISMRMRD::UserParameterDouble>::iterator i (traj_desc.userParameterDouble.begin()); i != traj_desc.userParameterDouble.end(); ++i) {
				if (i->name == "CS_Accel") {
					fCSAcc_ = i->value;
				} else if (i->name == "FullySampled") {
					fFullySa_ = i->value;
				} else if (i->name == "lambdaESPReSSo") {
					GlobalVar::instance()->cfLambdaESPReSSo_ = i->value;
				} else if (i->name == "lambda") {
					GlobalVar::instance()->cfLambda_ = i->value;
				}
			}
		}
#else
		if ((*e_seq.begin()).trajectoryDescription().present()) {
			GINFO("Trajectory description present!");
			ISMRMRD::trajectoryDescriptionType traj_desc = (*e_seq.begin()).trajectoryDescription().get();

			for (ISMRMRD::trajectoryDescriptionType::userParameterLong_sequence::iterator i (traj_desc.userParameterLong().begin()); i != traj_desc.userParameterLong().end(); ++i) {
				if (std::strcmp(i->name().c_str(),"WIP_PF_y") == 0) {
					iESPReSSoY = i->value();
				} else if (std::strcmp(i->name().c_str(),"WIP_PF_z") == 0) {
					iESPReSSoZ = i->value();
				} else if (std::strcmp(i->name().c_str(),"SamplingType") == 0) {
					iSamplingType_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"VDMap") == 0) {
					iVDMap_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"BodyRegion") == 0) {
					iBodyRegion_ = i->value();
				} else if (std::strcmp(i->name().c_str(), "OuterIterations") == 0) {
					GlobalVar::instance()->iNOuter_ = i->value();
				} else if (std::strcmp(i->name().c_str(), "InnerIterations") == 0) {
					GlobalVar::instance()->iNInner_ = i->value();
				} else if (std::strcmp(i->name().c_str(), "fftSparseDim") == 0) {
					GlobalVar::instance()->iDimFFT_ = i->value();
				} else if (std::strcmp(i->name().c_str(), "dctSparseDim") == 0) {
					GlobalVar::instance()->iDimDCTSparse_ = i->value();
				} else if (std::strcmp(i->name().c_str(), "pcaSparseDim") == 0) {
					GlobalVar::instance()->iDimPCASparse_ = i->value();
				} else if (std::strcmp(i->name().c_str(), "kernelFftDim") == 0) {
					GlobalVar::instance()->iDimKernelFFT_ = i->value();
				} else if (std::strcmp(i->name().c_str(), "transformFftBaDim") == 0) {
					GlobalVar::instance()->iTransformFFTBA_ = i->value();
				} else if (std::strcmp(i->name().c_str(), "kSpaceOutDim") == 0) {
					GlobalVar::instance()->ikSpaceOut_ = i->value();
				} else if (std::strcmp(i->name().c_str(), "CSESPReSSo") == 0) {
					GlobalVar::instance()->bESPRActiveCS_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"NavPeriod") == 0) {
					GlobalVar::instance()->iNavPeriod_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"NavPERes") == 0) {
					GlobalVar::instance()->iNavPERes_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"MeasurementTime") == 0) {
					GlobalVar::instance()->iMeasurementTime_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"Phases") == 0) {
					iNPhases_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"PopulationMode") == 0) {
					GlobalVar::instance()->iPopulationMode_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"GatingMode") == 0) {
					GlobalVar::instance()->iGatingMode_ = i->value();
				}
			}

			for (ISMRMRD::trajectoryDescriptionType::userParameterDouble_sequence::iterator i (traj_desc.userParameterDouble().begin()); i != traj_desc.userParameterDouble().end(); ++i) {
				if (std::strcmp(i->name().c_str(),"CS_Accel") == 0) {
					fCSAcc_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"FullySampled") == 0) {
					fFullySa_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"lambdaESPReSSo") == 0) {
					GlobalVar::instance()->cfLambdaESPReSSo_ = i->value();
				} else if (std::strcmp(i->name().c_str(),"lambda") == 0) {
					GlobalVar::instance()->cfLambda_ = i->value();
				}
			}
		}
#endif
		else {
			GWARN("No trajectory description present!\n");
		}

		//-------------------------------------------------------------------------
		//----------------------- Interpret Integer Data  -------------------------
		//-------------------------------------------------------------------------
		switch (iBodyRegion_) {
		case 14:
			GDEBUG("Body region is none\n");
			break;
		case 15:
			GDEBUG("Body region is head\n");
			break;
		case 16:
			GDEBUG("Body region is thorax\n");
			break;
		case 17:
			GDEBUG("Body region is abdomen\n");
			break;
		case 18:
			GDEBUG("Body region is pelvis\n");
			break;
		}

		switch (iVDMap_) {
		case 1:
			GDEBUG("VDMap is none\n");
			break;
		case 2:
			GDEBUG("VDMap is point\n");
			break;
		case 3:
			GDEBUG("VDMap is block\n");
			break;
		case 4:
			GDEBUG("VDMap is ellipse\n");
			break;
		case 5:
			GDEBUG("VDMap is ring\n");
			break;
		}

		switch (iSamplingType_) {
		case 6:
			GDEBUG("sampling type is poisson\n");
			break;
		case 7:
			GDEBUG("sampling type is random\n");
			break;
		case 8:
			GDEBUG("sampling type is proba\n");
			break;
		}

		iESPReSSoDirection_ = 10;
		fPartialFourierVal_ = 1.0;

		if ((iESPReSSoY > 9 && iESPReSSoY < 14) || (iESPReSSoZ > 9 && iESPReSSoZ < 14)) {
			GDEBUG("Partial Fourier data..\n");
			GDEBUG("ESPReSSo Y: %f, ESPReSSo Z: %f\n", iESPReSSoY, iESPReSSoZ);
			// get Partial Fourier dimension
			if (iESPReSSoY > 9) {
				iESPReSSoDirection_ = 1;
				// get Partial Fourier value
				switch (iESPReSSoY) {
				case 10:
					fPartialFourierVal_ = 0.5;
					break;
				case 11:
					fPartialFourierVal_ = 0.625;
					break;
				case 12:
					fPartialFourierVal_ = 0.75;
					break;
				case 13:
					fPartialFourierVal_ = 0.875;
					break;
				default:
					fPartialFourierVal_ = 1.0;
					break;
				}
			} else if (iESPReSSoZ > 9) {
				iESPReSSoDirection_ = 2;
				// get Partial Fourier value
				switch (iESPReSSoZ) {
				case 10:
					fPartialFourierVal_ = 0.5;
					break;
				case 11:
					fPartialFourierVal_ = 0.625;
					break;
				case 12:
					fPartialFourierVal_ = 0.75;
					break;
				case 13:
					fPartialFourierVal_ = 0.875;
					break;
				default:
					fPartialFourierVal_ = 1.0;
					break;
				}
			}
		}

		GDEBUG("Partial Fourier is %f \n", fPartialFourierVal_);
		GDEBUG("ESPReSSo Direction is %i \n", iESPReSSoDirection_);
	} catch(...) {
		GERROR("Cannot find CS entries in trajectory description..\n");
	}

	// concat higher and lower bytes from total measurement variable
	lNoScans_ = std::ceil(GlobalVar::instance()->iMeasurementTime_*1000/GlobalVar::instance()->fTR_);

	return GADGET_OK;
}

// process data - incoming unordered k-space(RO,t,c) --> ordered k-space(x,y,z,Phases,c)
int CS_Retro_AccumulatorGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	// set scan counter value
	const uint32_t current_scan = m1->getObjectPtr()->scan_counter-1;

	// protect Gadget from more inputs than expected
	if (current_scan > lNoScans_) {
		GDEBUG("Drop scan no. %d (unexpected)\n", current_scan);
		m1->release();

		return GADGET_OK;
	}

	// only init the buffers in case of real data acquisition (otherwise wrong values (e.g. base resolution) can occur)
	if (!(is_content_dataset(*m1->getObjectPtr()) || is_navigator_dataset(*m1->getObjectPtr()))) {
		GDEBUG("Reject scan with idx.set=%d, scan no. %d\n", m1->getObjectPtr()->idx.set, current_scan);
		m1->release();

		return GADGET_OK;
	}

	if (GlobalVar::instance()->AcqVec_.size() == 0) {
		// copy header information of first acquisition to global variable (create new header, copy header, push onto global vector)
		GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *tmp_m1 = new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();
		fCopyAcqHeader(tmp_m1, m1->getObjectPtr());
		GlobalVar::instance()->AcqVec_.push_back(tmp_m1->getObjectPtr());
	}

	/*---------------------------------------------------*/
	/*----------- init buffer for k-space data ----------*/
	/*---------------------------------------------------*/
	// add pointer to big pointer vector
	bufferkSpace_.push_back(m2);

	/*---------------------------------------------------*/
	/*----------- store incoming k-space data -----------*/
	/*---------------------------------------------------*/
	// navigator flag
	bool bNavigator = false;
	if (is_navigator_dataset(*m1->getObjectPtr())) {
		bNavigator = true;
	}

	// push current loop counters on according vector (temporal)
	uint16_t line = m1->getObjectPtr()->idx.kspace_encode_step_1;
	GlobalVar::instance()->vPE_.push_back(line);

	uint16_t partition = m1->getObjectPtr()->idx.kspace_encode_step_2;
	GlobalVar::instance()->vPA_.push_back(partition);

	if (bNavigator == true) {
		// add pointer to big pointer vector
		buffer_nav_.push_back(m2->getObjectPtr());

		iNoNavLine_++;

		if (iNoNavLine_ == GlobalVar::instance()->iNavPERes_) {
			iNoNav_++;
			iNoNavLine_ = 0;
		} else if (iNoNavLine_ == (GlobalVar::instance()->iNavPERes_/2)) {
			GlobalVar::instance()->vNavInd_.push_back(static_cast<float>(current_scan));
		}
	}

	/*---------------------------------------------------*/
	/*--------------- process sampled data --------------*/
	/*---------------------------------------------------*/
	if (current_scan == lNoScans_) {
		GINFO("data received.. try to process data\n");

		// crop non-empty data from navigator array
		// check if last measurement is in between a navigator block
		if (iNoNavLine_ != 0) {
			iNoNav_--;
			GlobalVar::instance()->vNavInd_.pop_back();
			buffer_nav_.pop_back();
		}

		// bring GlobalVar::instance()->vNavInd_ to length of iNoNav_ if necessary
		if (iNoNav_ > GlobalVar::instance()->vNavInd_.size()) {
			iNoNav_ = GlobalVar::instance()->vNavInd_.size();
			GINFO("iNoNav_ (=%d) larger than vNavInd_.size() (=%d). Correcting to equality.\n", iNoNav_, GlobalVar::instance()->vNavInd_.size());
		} else {
			// this could also be done by:
			// if (iNoNavLine_ >= (GlobalVar::instance()->iNavPERes_/2)) GlobalVar::instance()->vNavInd_.pop_back();
			while (iNoNav_ < GlobalVar::instance()->vNavInd_.size()) {
				GINFO("iNoNav_ (=%d) smaller than vNavInd_.size() (=%d). Correcting to equality.\n", iNoNav_, GlobalVar::instance()->vNavInd_.size());
				GlobalVar::instance()->vNavInd_.pop_back();
			}
		}

		GINFO("%i navigator data found..\n", iNoNav_);

		// create new ContainerMessages for header, navigator and kSpace data header
		GadgetContainerMessage<ISMRMRD::ImageHeader> *tmp_m1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();

		// initialize the image header
		memset(tmp_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));

		// initialize flags
		tmp_m1->getObjectPtr()->flags = 0;

		//tmp_m1->getObjectPtr()->user_int[0] = 7;
		tmp_m1->getObjectPtr()->user_int[0]			= iNPhases_;
		tmp_m1->getObjectPtr()->user_int[1]			= iBodyRegion_;
		tmp_m1->getObjectPtr()->user_int[2]			= iSamplingType_;
		tmp_m1->getObjectPtr()->user_int[3]			= iVDMap_;
		tmp_m1->getObjectPtr()->user_int[4]			= iESPReSSoDirection_;
		tmp_m1->getObjectPtr()->user_int[5]			= iNoNav_;
		tmp_m1->getObjectPtr()->user_float[0]		= fCSAcc_;
		tmp_m1->getObjectPtr()->user_float[1]		= fFullySa_/100;
		tmp_m1->getObjectPtr()->user_float[2]		= fPartialFourierVal_;
		tmp_m1->getObjectPtr()->matrix_size[0]		= dimensionsIn_[0];
		tmp_m1->getObjectPtr()->matrix_size[1]		= dimensionsIn_[1];
		tmp_m1->getObjectPtr()->matrix_size[2]		= dimensionsIn_[2];
		tmp_m1->getObjectPtr()->field_of_view[0]	= field_of_view_[0];
		tmp_m1->getObjectPtr()->field_of_view[1]	= field_of_view_[1];
		tmp_m1->getObjectPtr()->field_of_view[2]	= field_of_view_[2];
		tmp_m1->getObjectPtr()->channels			= static_cast<uint16_t>(m1->getObjectPtr()->active_channels);
		tmp_m1->getObjectPtr()->slice				= m1->getObjectPtr()->idx.slice;
		memcpy(tmp_m1->getObjectPtr()->position,m1->getObjectPtr()->position,sizeof(float)*3);
		memcpy(tmp_m1->getObjectPtr()->read_dir,m1->getObjectPtr()->read_dir,sizeof(float)*3);
		memcpy(tmp_m1->getObjectPtr()->phase_dir,m1->getObjectPtr()->phase_dir,sizeof(float)*3);
		memcpy(tmp_m1->getObjectPtr()->slice_dir,m1->getObjectPtr()->slice_dir, sizeof(float)*3);
		memcpy(tmp_m1->getObjectPtr()->patient_table_position,m1->getObjectPtr()->patient_table_position, sizeof(float)*3);
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		tmp_m1->getObjectPtr()->data_type		= ISMRMRD::ISMRMRD_CXFLOAT;
#else
		tmp_m1->getObjectPtr()->image_data_type	= ISMRMRD::DATA_COMPLEX_FLOAT;
#endif
		tmp_m1->getObjectPtr()->image_index = 1;
		tmp_m1->getObjectPtr()->image_series_index = 0;

		// delete header - it is not needed anymore
		m1->cont(NULL);
		m1->release();

		// navigator
		hoNDArray<std::complex<float> > total_nav_array;

		// create output
		try {
			std::vector<size_t> nav_dims;
			nav_dims.push_back(m1->getObjectPtr()->number_of_samples);	// get number of samples in acquisition (equals base resolution)
			nav_dims.push_back(m1->getObjectPtr()->active_channels);
			nav_dims.push_back(GlobalVar::instance()->iNavPERes_);
			nav_dims.push_back(buffer_nav_.size()/GlobalVar::instance()->iNavPERes_);
			total_nav_array.create(&nav_dims);
		} catch (std::runtime_error &err) {
			GEXCEPTION(err, "Unable to allocate new image array\n");

			tmp_m1->release();

			return GADGET_FAIL;
		}

		// copy data
		for (size_t i = 0; i < buffer_nav_.size(); i++) {
			size_t offset = (i/total_nav_array.get_size(2))*total_nav_array.get_size(0)*total_nav_array.get_size(1)*total_nav_array.get_size(2)
				+ (i % total_nav_array.get_size(2))*total_nav_array.get_size(0)*total_nav_array.get_size(1);
			memcpy(total_nav_array.get_data_ptr()+offset, buffer_nav_.at(i)->get_data_ptr(), buffer_nav_.at(i)->get_number_of_bytes());
		}

		// permute nav: [baseRes channels NavPERes scans] -> [baseRes scans NavPEREs channels]
		std::vector<size_t> new_nav_dim;
		new_nav_dim.push_back(0);
		new_nav_dim.push_back(3);
		new_nav_dim.push_back(2);
		new_nav_dim.push_back(1);
		total_nav_array = *permute(&total_nav_array, &new_nav_dim, false);
		total_nav_array.delete_data_on_destruct(false);

		// create nav output message
		GadgetContainerMessage<hoNDArray<std::complex<float> > > *tmp_m2 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();
		tmp_m2->getObjectPtr()->create(total_nav_array.get_size(0), total_nav_array.get_size(1), total_nav_array.get_size(2), total_nav_array.get_size(3), total_nav_array.get_data_ptr(), true);

		// concatenate data with header
		tmp_m1->cont(tmp_m2);

		// kSpace
		hoNDArray<std::complex<float> > total_kspace_array;

		// create output
		try {
			std::vector<size_t> kspace_dims;
			for (size_t i = 0; i < bufferkSpace_.at(0)->getObjectPtr()->get_number_of_dimensions(); i++) {
				kspace_dims.push_back(bufferkSpace_.at(0)->getObjectPtr()->get_size(i));
			}
			kspace_dims.push_back(bufferkSpace_.size());
			total_kspace_array.create(&kspace_dims);
		} catch (std::runtime_error &err) {
			GEXCEPTION(err, "Unable to allocate new image array\n");

			tmp_m1->release();

			return GADGET_FAIL;
		}

		// copy data
		size_t kspace_elements_copied = 0;
		size_t kspace_elements_to_copy = 0;
		for (size_t i = 0; i < bufferkSpace_.size(); i++) {
			kspace_elements_to_copy = bufferkSpace_.at(i)->getObjectPtr()->get_number_of_elements();
			memcpy(total_kspace_array.get_data_ptr()+kspace_elements_copied, bufferkSpace_.at(i)->getObjectPtr()->get_data_ptr(), bufferkSpace_.at(i)->getObjectPtr()->get_number_of_bytes());
			kspace_elements_copied += kspace_elements_to_copy;

			// free memory
			bufferkSpace_.at(i)->release();
			bufferkSpace_.at(i) = NULL;
		}

		// permute kspace: [baseRes channels scans] -> [baseRes scans channels]
		std::vector<size_t> new_kspace_dim;
		new_kspace_dim.push_back(0);
		new_kspace_dim.push_back(2);
		new_kspace_dim.push_back(1);
		total_kspace_array = *permute(&total_kspace_array, &new_kspace_dim, false);
		total_kspace_array.delete_data_on_destruct(false);

		// create kspace output message
		GadgetContainerMessage<hoNDArray<std::complex<float> > > *tmp_m3 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();
		tmp_m3->getObjectPtr()->create(total_kspace_array.get_size(0), total_kspace_array.get_size(1), total_kspace_array.get_size(2), total_kspace_array.get_size(3), total_kspace_array.get_data_ptr(), true);

		// concatenate data
		tmp_m2->cont(tmp_m3);

		// put on stream
		if (this->next()->putq(tmp_m1) < 0) {
			return GADGET_FAIL;
		}

		GDEBUG("global PE: %i, PA: %i\n", GlobalVar::instance()->vPE_.size(), GlobalVar::instance()->vPA_.size());

		return GADGET_OK;
	}

	// free memory (only m1, we need m2 later on (saved in bufferkSpace_)
	m1->cont(NULL);
	m1->release();

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_AccumulatorGadget)
