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

namespace Gadgetron{
	// class constructor
	CS_Retro_AccumulatorGadget::CS_Retro_AccumulatorGadget() : bufferkSpace_(0), bufferNav_(0), iBaseRes_(0), fFullySa_(.065), iEchoLine_(0), iEchoPartition_(0), lNoScans_(0), iNoSamples_(0), iNoChannels_(0) {
		GlobalVar::instance()->vPE_.clear();
		GlobalVar::instance()->vPA_.clear();
	}

	// class destructor - delete temporal buffer/memory
	CS_Retro_AccumulatorGadget::~CS_Retro_AccumulatorGadget(){
		if (bufferkSpace_)
			delete bufferkSpace_;

		if (bufferNav_)
			delete bufferNav_;
	}

	// read flexible data header
	int CS_Retro_AccumulatorGadget::process_config(ACE_Message_Block* mb)
	{
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		ISMRMRD::IsmrmrdHeader h;
		ISMRMRD::deserialize(mb->rd_ptr(),h);

		ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
		ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
		ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

		// get FOV
		field_of_view_.push_back(r_space.fieldOfView_mm.x);
		field_of_view_.push_back(e_space.fieldOfView_mm.y);
		field_of_view_.push_back(e_space.fieldOfView_mm.z);

		// get matrix size
		dimensionsIn_.push_back(r_space.matrixSize.x);
		dimensionsIn_.push_back(e_space.matrixSize.y);
		dimensionsIn_.push_back(e_space.matrixSize.z);

		iEchoLine_ = e_limits.kspace_encoding_step_1.get().center;
		iEchoPartition_ = e_limits.kspace_encoding_step_2.get().center;

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
		dimensionsIn_.push_back(r_space.matrixSize().x());
		dimensionsIn_.push_back(e_space.matrixSize().y());
		dimensionsIn_.push_back(e_space.matrixSize().z());
		GADGET_DEBUG2("Matrix size: %d, %d, %d\n", dimensionsIn_[0], dimensionsIn_[1], dimensionsIn_[2]);

		// get FOV
		field_of_view_.push_back(r_space.fieldOfView_mm().x());
		field_of_view_.push_back(e_space.fieldOfView_mm().y());
		field_of_view_.push_back(e_space.fieldOfView_mm().z());
		GADGET_DEBUG2("FOV: %f, %f, %f\n", r_space.fieldOfView_mm().x(), e_space.fieldOfView_mm().y(), e_space.fieldOfView_mm().z());

		// get echo line and echo partition
		iEchoLine_ = e_limits.kspace_encoding_step_1().get().center();
		iEchoPartition_ = e_limits.kspace_encoding_step_2().get().center();
		GADGET_DEBUG2("echo line: %i, echo partition: %i", iEchoLine_, iEchoPartition_);

		// repetition time
		GlobalVar::instance()->fTR_ = cfg->sequenceParameters().get().TR().at(0);
	#endif

		// set properties
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		GlobalVar::instance()->iNavPeriod_			= NavPeriod.value();
		GlobalVar::instance()->iNavPERes_			= NavPERes.value();
		GlobalVar::instance()->iMeasurementTime_	= MeasurementTime.value();
		GlobalVar::instance()->iNPhases_			= Phases.value();
		GlobalVar::instance()->iPopulationMode_		= PopulationMode.value();
		GlobalVar::instance()->iGatingMode_			= GatingMode.value();
	#else
		GlobalVar::instance()->iNavPeriod_			= *(get_int_value("NavPeriod").get();
		GlobalVar::instance()->iNavPERes_			= *(get_int_value("NavPERes").get();
		GlobalVar::instance()->iMeasurementTime_	= *(get_int_value("MeasurementTime").get();
		GlobalVar::instance()->iNPhases_			= *(get_int_value("Phases").get();
		GlobalVar::instance()->iPopulationMode_		= *(get_int_value("PopulationMode").get();
		GlobalVar::instance()->iGatingMode_			= *(get_int_value("GatingMode").get();
	#endif

		int iESPReSSoY = 0;
		int iESPReSSoZ = 0;
		iBodyRegion_ = 0;
		GlobalVar::instance()->iVDMap_ = 0;
		iSamplingType_ = 0;
		GlobalVar::instance()->fCSAcc_ = 0;
		fFullySa_ = 0;
		try {
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
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
						GlobalVar::instance()->iVDMap_ = i->value;
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
					} else if (i->name == "MeasurementTime"){
						GlobalVar::instance()->iMeasurementTime_ = i->value;
					} else if (i->name == "Phases") {
						GlobalVar::instance()->iNPhases_ = i->value;
					} else if (i->name == "PopulationMode") {
						GlobalVar::instance()->iPopulationMode_ = i->value;
					} else if (i->name == "GatingMode") {
						GlobalVar::instance()->iGatingMode_ = i->value;
					}
				}

				for (std::vector<ISMRMRD::UserParameterDouble>::iterator i (traj_desc.userParameterDouble.begin()); i != traj_desc.userParameterDouble.end(); ++i) {
					if (i->name == "CS_Accel") {
						GlobalVar::instance()->fCSAcc_ = i->value;
					} else if (i->name == "FullySampled") {
						GlobalVar::instance()->fFullySampled_ = i->value;
					} else if (i->name == "lambdaESPReSSo") {
						GlobalVar::instance()->cfLambdaESPReSSo_ = i->value;
					} else if (i->name == "lambda") {
						GlobalVar::instance()->cfLambda_ = i->value;
					}
				}
			}
	#else
			if ((*e_seq.begin()).trajectoryDescription().present()) {
				GADGET_DEBUG1("\n\nTrajectory description present!\n\n");
				ISMRMRD::trajectoryDescriptionType traj_desc = (*e_seq.begin()).trajectoryDescription().get();

				for (ISMRMRD::trajectoryDescriptionType::userParameterLong_sequence::iterator i (traj_desc.userParameterLong().begin ()); i != traj_desc.userParameterLong().end(); ++i) {
					if (std::strcmp(i->name().c_str(),"WIP_PF_y") == 0) {
						iESPReSSoY = i->value();
					} else if (std::strcmp(i->name().c_str(),"WIP_PF_z") == 0) {
						iESPReSSoZ = i->value();
					} else if (std::strcmp(i->name().c_str(),"SamplingType") == 0) {
						iSamplingType_ = i->value();
					} else if (std::strcmp(i->name().c_str(),"VDMap") == 0) {
						GlobalVar::instance()->iVDMap_ = i->value();
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
						GlobalVar::instance()->iNPhases_ = i->value();
					} else if (std::strcmp(i->name().c_str(),"PopulationMode") == 0) {
						GlobalVar::instance()->iPopulationMode_ = i->value();
					} else if (std::strcmp(i->name().c_str(),"GatingMode") == 0) {
						GlobalVar::instance()->iGatingMode_ = i->value();
					}
				}

				for (ISMRMRD::trajectoryDescriptionType::userParameterDouble_sequence::iterator i (traj_desc.userParameterDouble().begin ()); i != traj_desc.userParameterDouble().end(); ++i) {
					if (std::strcmp(i->name().c_str(),"CS_Accel") == 0) {
						GlobalVar::instance()->fCSAcc_ = i->value();
					} else if (std::strcmp(i->name().c_str(),"FullySampled") == 0) {
						GlobalVar::instance()->fFullySampled_ = i->value();
					} else if (std::strcmp(i->name().c_str(),"lambdaESPReSSo") == 0) {
						GlobalVar::instance()->cfLambdaESPReSSo_ = i->value();
					} else if (std::strcmp(i->name().c_str(),"lambda") == 0) {
						GlobalVar::instance()->cfLambda_ = i->value();
					}
				}
			}
	#endif
			else{
				GADGET_DEBUG1("\n\nNo trajectory description present!\n\n");
			}

			//-------------------------------------------------------------------------
			//----------------------- Interpret Integer Data  -------------------------
			//-------------------------------------------------------------------------
			switch (iBodyRegion_){
				case 14:
					GADGET_DEBUG1("Body region is none\n");
					break;
				case 15:
					GADGET_DEBUG1("Body region is head\n");
					break;
				case 16:
					GADGET_DEBUG1("Body region is thorax\n");
					break;
				case 17:
					GADGET_DEBUG1("Body region is abdomen\n");
					break;
				case 18:
					GADGET_DEBUG1("Body region is pelvis\n");
					break;
			}

			switch (GlobalVar::instance()->iVDMap_){
				case 1:
					GADGET_DEBUG1("VDMap is none\n");
					break;
				case 2:
					GADGET_DEBUG1("VDMap is point\n");
					break;
				case 3:
					GADGET_DEBUG1("VDMap is block\n");
					break;
				case 4:
					GADGET_DEBUG1("VDMap is ellipse\n");
					break;
				case 5:
					GADGET_DEBUG1("VDMap is ring\n");
					break;
			}

			switch (iSamplingType_){
				case 6:
					GADGET_DEBUG1("sampling type is poisson\n");
					break;
				case 7:
					GADGET_DEBUG1("sampling type is random\n");
					break;
				case 8:
					GADGET_DEBUG1("sampling type is proba\n");
					break;
			}

			GlobalVar::instance()->iESPReSSoDirection_ = 10;
			GlobalVar::instance()->fPartialFourierVal_ = 1.0;

			if ((iESPReSSoY > 9 && iESPReSSoY < 14) || (iESPReSSoZ > 9 && iESPReSSoZ < 14)) {
				GADGET_DEBUG1("Partial Fourier data..\n");
				GADGET_DEBUG2("ESPReSSo Y: %f, ESPReSSo Z: %f\n", iESPReSSoY, iESPReSSoZ);
				// get Partial Fourier dimension
				if (iESPReSSoY > 9){
					GlobalVar::instance()->iESPReSSoDirection_ = 1;
					// get Partial Fourier value
					switch (iESPReSSoY){
						case 10:
							GlobalVar::instance()->fPartialFourierVal_ = 0.5;
							break;
						case 11:
							GlobalVar::instance()->fPartialFourierVal_ = 0.625;
							break;
						case 12:
							GlobalVar::instance()->fPartialFourierVal_ = 0.75;
							break;
						case 13:
							GlobalVar::instance()->fPartialFourierVal_ = 0.875;
							break;
						default:
							GlobalVar::instance()->fPartialFourierVal_ = 1.0;
							break;
					}
				}

				else if (iESPReSSoZ > 9){
					GlobalVar::instance()->iESPReSSoDirection_ = 2;
					// get Partial Fourier value
					switch (iESPReSSoZ){
						case 10:
							GlobalVar::instance()->fPartialFourierVal_ = 0.5;
							break;
						case 11:
							GlobalVar::instance()->fPartialFourierVal_ = 0.625;
							break;
						case 12:
							GlobalVar::instance()->fPartialFourierVal_ = 0.75;
							break;
						case 13:
							GlobalVar::instance()->fPartialFourierVal_ = 0.875;
							break;
						default:
							GlobalVar::instance()->fPartialFourierVal_ = 1.0;
							break;
					}
				}
			}

			GADGET_DEBUG2("Partial Fourier is %f \n", GlobalVar::instance()->fPartialFourierVal_);
			GADGET_DEBUG2("ESPReSSo Direction is %i \n", GlobalVar::instance()->iESPReSSoDirection_);
		} catch(...) {
			GADGET_DEBUG1("Error occured - cannot find CS entries in trajectory description..\n");
		}

		// concat higher and lower bytes from total measurement variable
		lNoScans_ = std::ceil(GlobalVar::instance()->iMeasurementTime_*1000/GlobalVar::instance()->fTR_);

		return GADGET_OK;
	}

	// process data - incoming unordered k-space(RO,t,c) --> ordered k-space(x,y,z,Phases,c)
	int CS_Retro_AccumulatorGadget::process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2) {
		/*---------------------------------------------------*/
		/*----------- init buffer for k-space data ----------*/
		/*---------------------------------------------------*/
		if (!bufferkSpace_) {
			iNoChannels_ = m1->getObjectPtr()->active_channels;

			// initialize k-space buffer
			if (!(bufferkSpace_ = new hoNDArray< std::complex<float> >())) {
				GADGET_DEBUG1("Failed to create k-space buffer\n");
				return GADGET_FAIL;
			}

			// get number of samples in acquisition (equals base resolution)
			iBaseRes_ = m1->getObjectPtr()->number_of_samples;

			GADGET_DEBUG2("base res.: %d, no. scans: %lu, no. channel: %u\n", iBaseRes_, lNoScans_, iNoChannels_);

			// dimension vector of k-space array
			dimkSpace_.push_back(iBaseRes_);
			dimkSpace_.push_back(lNoScans_);
			dimkSpace_.push_back(iNoChannels_);

			// create buffer array for incoming k-space data (readout, time, channel)
			try {
				bufferkSpace_->create(&dimkSpace_);
			} catch (std::runtime_error &err) {
				GADGET_DEBUG_EXCEPTION(err, "Failed to allocate k-space buffer array\n");
				return GADGET_FAIL;
			}

			// copy header information of first acquisition to global variable (create new header, set memory zero, copy header, push onto global vector)
			GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* tmp_m1 = new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();
			memset(tmp_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));
			fCopyAcqHeader(tmp_m1, m1);
			GlobalVar::instance()->AcqVec_.push_back(tmp_m1);
			GADGET_DEBUG1("Receiving data..\n");
			GADGET_DEBUG1("\n \n bufferkSpace_ \n \n:");
		}

		/*---------------------------------------------------*/
		/*--------- init buffer for navigator data ----------*/
		/*---------------------------------------------------*/
		if (!bufferNav_) {
			// initialize k-space buffer
			if (!(bufferNav_ = new hoNDArray< std::complex<float> >())) {
				GADGET_DEBUG1("Failed to create navigator buffer\n");
				return GADGET_FAIL;
			}

			// dimension vector of navigator array
			dimNav_.push_back(iBaseRes_);
			dimNav_.push_back(lNoScans_);
			dimNav_.push_back(GlobalVar::instance()->iNavPERes_);
			dimNav_.push_back(iNoChannels_);
			iNoNav_ = 0;
			iNoNavLine_ = 0;

			GADGET_DEBUG2("navigator dimensions: base res: %d, no. scans: %lu, PE resolution: %d, no. channels: %u\n", iBaseRes_, lNoScans_, GlobalVar::instance()->iNavPERes_, iNoChannels_);

			// create buffer array for incoming navigator data (readout, time, PE, channel)
			try {
				bufferNav_->create(&dimNav_);
			} catch (std::runtime_error &err) {
				GADGET_DEBUG_EXCEPTION(err, "Failed to allocate navigator buffer array\n");
				return GADGET_FAIL;
			}

			GADGET_DEBUG1("\n \n bufferNav \n \n:");

			bufferNav_->print(std::cout);
			lCurrentScan_ = 0;
		}

		/*---------------------------------------------------*/
		/*----------- store incoming k-space data -----------*/
		/*---------------------------------------------------*/
		// navigator flag
		bool  bNavigator = false;
		if (m1->getObjectPtr()->user_int[1] & 0x1 != 0){
			bNavigator = true;
		}

		// get current loop counters
		int samples		= m1->getObjectPtr()->number_of_samples;
		int line		= m1->getObjectPtr()->idx.kspace_encode_step_1;
		int partition	= m1->getObjectPtr()->idx.kspace_encode_step_2;
		int slice		= m1->getObjectPtr()->idx.slice;
		int repetition	= m1->getObjectPtr()->idx.repetition;
		int set			= m1->getObjectPtr()->idx.set;

		// push current loop counters on according vector (temporal)
		GlobalVar::instance()->vPE_.push_back(line); GlobalVar::instance()->vPA_.push_back(partition);

		// get data pointer
		std::complex<float> *pkSpace	= bufferkSpace_->get_data_ptr();
		std::complex<float> *pNav		= bufferNav_->get_data_ptr();
		std::complex<float>	*pIncoming	= m2->getObjectPtr()->get_data_ptr();

		for (int c = 0; c < m1->getObjectPtr()->active_channels; c++) {
			size_t offset_kSpace = c*dimkSpace_[0]*dimkSpace_[1] + lCurrentScan_*dimkSpace_[0] + (dimkSpace_[0]>>1)-m1->getObjectPtr()->center_sample;
			memcpy(pkSpace + offset_kSpace, pIncoming+c*samples, sizeof(std::complex<float>)*samples);

			if (bNavigator == true) {
				size_t offset_Nav = c*dimNav_[0]*dimNav_[1]*dimNav_[2] + iNoNavLine_*dimNav_[0]*dimNav_[1] + iNoNav_*dimNav_[0]+(dimNav_[0]>>1)-m1->getObjectPtr()->center_sample;
				memcpy(pNav + offset_Nav, pIncoming+c*samples, sizeof(std::complex<float>)*samples);
			}
		}

		if (bNavigator == true) {
			iNoNavLine_++;

			if (iNoNavLine_ == GlobalVar::instance()->iNavPERes_) {
				iNoNav_++;
				iNoNavLine_ = 0;
			} else if (iNoNavLine_ == (GlobalVar::instance()->iNavPERes_/2)){
				GlobalVar::instance()->vNavInd_.push_back(static_cast<float>(lCurrentScan_));
			}
		}

		lCurrentScan_++;

		/*---------------------------------------------------*/
		/*--------------- process sampled data --------------*/
		/*---------------------------------------------------*/
		if (lCurrentScan_ == lNoScans_) {
			GADGET_DEBUG1("data received.. try to process data\n");

			// crop non-empty data from navigator array
			// check if last measurement is in between a navigator block
			if (iNoNavLine_!=0) {
				iNoNav_--;
				GlobalVar::instance()->vNavInd_.pop_back();
			}

			GADGET_DEBUG2("%i navigator data found..\n", iNoNav_);

			std::vector<size_t> vStart, vSize;
			vStart.push_back(0);  vStart.push_back(0); vStart.push_back(0); vStart.push_back(0);
			vSize.push_back(bufferNav_->get_size(0)); vSize.push_back(iNoNav_); vSize.push_back(bufferNav_->get_size(2)); vSize.push_back(bufferNav_->get_size(3));
			bufferNav_->print(std::cout);
			get_subarray(*bufferNav_, vStart, vSize, *bufferNav_);
			bufferNav_->print(std::cout);

			GADGET_DEBUG1("Flag - last in measurement detected..\n");
			GADGET_DEBUG1("\n------------- navigator buffer ------------------\n");
			bufferNav_->print(std::cout);
			GADGET_DEBUG1("\n------------- kSpace buffer -----------------\n");
			bufferkSpace_->print(std::cout);

			// create new ContainerMessages for header, navigator and kSpace data header
			GadgetContainerMessage<ISMRMRD::ImageHeader>* tmp_m1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();

			// initialize the image header
			memset(tmp_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));

			// initialize flags
			tmp_m1->getObjectPtr()->flags = 0;

			//tmp_m1->getObjectPtr()->user_int[0] = 7;
			tmp_m1->getObjectPtr()->user_int[0]			= GlobalVar::instance()->iNPhases_;
			tmp_m1->getObjectPtr()->user_int[1]			= m1->getObjectPtr()->user_int[1];
			tmp_m1->getObjectPtr()->user_int[2]			= m1->getObjectPtr()->user_int[2];
			tmp_m1->getObjectPtr()->user_int[3]			= m1->getObjectPtr()->user_int[3];
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
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			tmp_m1->getObjectPtr()->data_type		= ISMRMRD::ISMRMRD_CXFLOAT;
	#else
			tmp_m1->getObjectPtr()->image_data_type	= ISMRMRD::DATA_COMPLEX_FLOAT;
	#endif
			tmp_m1->getObjectPtr()->image_index = static_cast<uint16_t>(++image_counter_);
			tmp_m1->getObjectPtr()->image_series_index = static_cast<uint16_t>(image_series_);

			// set user values, if Compressed Sensing is active
			//if(this->get_bool_value("CS_on") == true){
			tmp_m1->getObjectPtr()->user_float[0] = GlobalVar::instance()->fCSAcc_;
			tmp_m1->getObjectPtr()->user_float[1] = fFullySa_/100;
			tmp_m1->getObjectPtr()->user_float[2] = GlobalVar::instance()->fPartialFourierVal_;
			tmp_m1->getObjectPtr()->user_float[3] = fLQ_;
			tmp_m1->getObjectPtr()->user_float[4] = fLESPReSSo_;

			tmp_m1->getObjectPtr()->user_int[1] = iBodyRegion_;
			tmp_m1->getObjectPtr()->user_int[2] = iSamplingType_;
			tmp_m1->getObjectPtr()->user_int[3] = GlobalVar::instance()->iVDMap_;
			tmp_m1->getObjectPtr()->user_int[4] = GlobalVar::instance()->iESPReSSoDirection_;
			//}

			tmp_m1->getObjectPtr()->user_int[5] = iNoNav_;

			// navigator
			GadgetContainerMessage<hoNDArray<std::complex<float>>>* tmp_m2 = new GadgetContainerMessage<hoNDArray<std::complex<float>>>();

			// concatenate data with header
			tmp_m1->cont(tmp_m2);

			// create output
			try {
				tmp_m2->getObjectPtr()->create(bufferNav_->get_dimensions());
			} catch (std::runtime_error &err) {
				GADGET_DEBUG_EXCEPTION(err, "Unable to allocate new image array\n");

				tmp_m1->release();

				return -1;
			}

			// copy data
			memcpy(tmp_m2->getObjectPtr()->get_data_ptr(), bufferNav_->get_data_ptr(), sizeof(std::complex<float>)*bufferNav_->get_number_of_elements());

			// kSpace
			GadgetContainerMessage<hoNDArray<std::complex<float> > > *tmp_m3 = new GadgetContainerMessage<hoNDArray<std::complex<float>>>();

			// concatenate data
			tmp_m2->cont(tmp_m3);

			// create output
			try {
				tmp_m3->getObjectPtr()->create(bufferkSpace_->get_dimensions());
			} catch (std::runtime_error &err) {
				GADGET_DEBUG_EXCEPTION(err, "Unable to allocate new image array\n");

				tmp_m1->release();

				return -1;
			}

			// copy data
			memcpy(tmp_m3->getObjectPtr()->get_data_ptr(), bufferkSpace_->get_data_ptr(), sizeof(std::complex<float>)*bufferkSpace_->get_number_of_elements());

			// put on stream
			if (this->next()->putq(tmp_m1) < 0) {
				return GADGET_FAIL;
			}

			GADGET_DEBUG2("global PE: %i, PA: %i\n", GlobalVar::instance()->vPE_.size(), GlobalVar::instance()->vPA_.size());

			return GADGET_OK;
		}
	}

	GADGET_FACTORY_DECLARE(CS_Retro_AccumulatorGadget)
} // close namespace Gadgetron
