﻿/*
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
	// Deletion of memory causes problems.
	// TODO: backtrace them and delete all the allocated memory!
// 		if (bufferkSpace_)
// 			delete bufferkSpace_;
//
// 		if (bufferNav_)
// 			delete bufferNav_;
}

// print AcquisitionHeader (for debugging purposes only)
void print_header(int counter, ISMRMRD::AcquisitionHeader h)
{
	static ISMRMRD::AcquisitionHeader last_header;
	static bool header_set = false;

	std::cout << "print header " << counter << ":" << std::endl;
	if (!header_set) {
		// print all
		std::cout << "version=" << h.version << ", ";
		std::cout << "flags=" << h.flags << ", ";
		std::cout << "measurement_uid=" << h.measurement_uid << ", ";
		std::cout << "scan_counter=" << h.scan_counter << ", ";
		std::cout << "acquisition_time_stamp=" << h.acquisition_time_stamp << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++)
			std::cout << "physiology_time_stamp[" << i << "]=" << h.physiology_time_stamp[i] << ", ";

		std::cout << "number_of_samples=" << h.number_of_samples << ", ";
		std::cout << "available_channels=" << h.available_channels << ", ";
		std::cout << "active_channels=" << h.active_channels << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_CHANNEL_MASKS; i++)
			std::cout << "channel_mask[" << i << "]=" << h.channel_mask[i] << ", ";

		std::cout << "discard_pre=" << h.discard_pre << ", ";
		std::cout << "discard_post=" << h.discard_post << ", ";
		std::cout << "center_sample=" << h.center_sample << ", ";
		std::cout << "encoding_space_ref=" << h.encoding_space_ref << ", ";
		std::cout << "trajectory_dimensions=" << h.trajectory_dimensions << ", ";
		std::cout << "sample_time_us=" << h.sample_time_us << ", ";

		for (int i = 0; i < 3; i++)
			std::cout << "position[" << i << "]=" << h.position[i] << ", ";
		for (int i = 0; i < 3; i++)
			std::cout << "read_dir[" << i << "]=" << h.read_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			std::cout << "phase_dir[" << i << "]=" << h.phase_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			std::cout << "slice_dir[" << i << "]=" << h.slice_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			std::cout << "patient_table_position[" << i << "]=" << h.patient_table_position[i] << ", ";

		std::cout << "idx.kspace_encode_step_1=" << h.idx.kspace_encode_step_1 << ", ";
		std::cout << "idx.kspace_encode_step_2=" << h.idx.kspace_encode_step_2 << ", ";
		std::cout << "idx.average=" << h.idx.average << ", ";
		std::cout << "idx.slice=" << h.idx.slice << ", ";
		std::cout << "idx.contrast=" << h.idx.contrast << ", ";
		std::cout << "idx.phase=" << h.idx.phase << ", ";
		std::cout << "idx.repetition=" << h.idx.repetition << ", ";
		std::cout << "idx.set=" << h.idx.set << ", ";
		std::cout << "idx.segment=" << h.idx.segment << ", ";
		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++)
			std::cout << "idx.user[" << i << "]=" << h.idx.user[i] << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++)
			std::cout << "user_int[" << i << "]=" << h.user_int[i] << ", ";
		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++)
			std::cout << "user_float[" << i << "]=" << h.user_float[i] << ", ";
	} else {
		// print differences
		if (h.version != last_header.version)
			std::cout << "version=" << h.version << ", ";
		if (h.flags != last_header.flags)
			std::cout << "flags=" << h.flags << ", ";
		if (h.measurement_uid != last_header.measurement_uid)
			std::cout << "measurement_uid=" << h.measurement_uid << ", ";
		if (h.scan_counter != last_header.scan_counter)
			std::cout << "scan_counter=" << h.scan_counter << ", ";
		if (h.acquisition_time_stamp != last_header.acquisition_time_stamp)
			std::cout << "acquisition_time_stamp=" << h.acquisition_time_stamp << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++)
			if (h.physiology_time_stamp[i] != last_header.physiology_time_stamp[i])
				std::cout << "physiology_time_stamp[" << i << "]=" << h.physiology_time_stamp[i] << ", ";

		if (h.number_of_samples != last_header.number_of_samples)
			std::cout << "number_of_samples=" << h.number_of_samples << ", ";
		if (h.available_channels != last_header.available_channels)
			std::cout << "available_channels=" << h.available_channels << ", ";
		if (h.active_channels != last_header.active_channels)
			std::cout << "active_channels=" << h.active_channels << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_CHANNEL_MASKS; i++)
			if (h.channel_mask[i] != last_header.channel_mask[i])
				std::cout << "channel_mask[" << i << "]=" << h.channel_mask[i] << ", ";

		if (h.discard_pre != last_header.discard_pre)
			std::cout << "discard_pre=" << h.discard_pre << ", ";
		if (h.discard_post != last_header.discard_post)
			std::cout << "discard_post=" << h.discard_post << ", ";
		if (h.center_sample != last_header.center_sample)
			std::cout << "center_sample=" << h.center_sample << ", ";
		if (h.encoding_space_ref != last_header.encoding_space_ref)
			std::cout << "encoding_space_ref=" << h.encoding_space_ref << ", ";
		if (h.trajectory_dimensions != last_header.trajectory_dimensions)
			std::cout << "trajectory_dimensions=" << h.trajectory_dimensions << ", ";
		if (h.sample_time_us != last_header.sample_time_us)
			std::cout << "sample_time_us=" << h.sample_time_us << ", ";

		for (int i = 0; i < 3; i++)
			if (h.position[i] != last_header.position[i])
				std::cout << "position[" << i << "]=" << h.position[i] << ", ";
		for (int i = 0; i < 3; i++)
			if (h.read_dir[i] != last_header.read_dir[i])
				std::cout << "read_dir[" << i << "]=" << h.read_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			if (h.phase_dir[i] != last_header.phase_dir[i])
				std::cout << "phase_dir[" << i << "]=" << h.phase_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			if (h.slice_dir[i] != last_header.slice_dir[i])
				std::cout << "slice_dir[" << i << "]=" << h.slice_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			if (h.patient_table_position[i] != last_header.patient_table_position[i])
				std::cout << "patient_table_position[" << i << "]=" << h.patient_table_position[i] << ", ";

		if (h.idx.kspace_encode_step_1 != last_header.idx.kspace_encode_step_1)
			std::cout << "idx.kspace_encode_step_1=" << h.idx.kspace_encode_step_1 << ", ";
		if (h.idx.kspace_encode_step_2 != last_header.idx.kspace_encode_step_2)
			std::cout << "idx.kspace_encode_step_2=" << h.idx.kspace_encode_step_2 << ", ";
		if (h.idx.average != last_header.idx.average)
			std::cout << "idx.average=" << h.idx.average << ", ";
		if (h.idx.slice != last_header.idx.slice)
			std::cout << "idx.slice=" << h.idx.slice << ", ";
		if (h.idx.contrast != last_header.idx.contrast)
			std::cout << "idx.contrast=" << h.idx.contrast << ", ";
		if (h.idx.phase != last_header.idx.phase)
			std::cout << "idx.phase=" << h.idx.phase << ", ";
		if (h.idx.repetition != last_header.idx.repetition)
			std::cout << "idx.repetition=" << h.idx.repetition << ", ";
		if (h.idx.set != last_header.idx.set)
			std::cout << "idx.set=" << h.idx.set << ", ";
		if (h.idx.segment != last_header.idx.segment)
			std::cout << "idx.segment=" << h.idx.segment << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++)
			if (h.idx.user[i] != last_header.idx.user[i])
				std::cout << "idx.user[" << i << "]=" << h.idx.user[i] << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++)
			if (h.user_int[i] != last_header.user_int[i])
				std::cout << "user_int[" << i << "]=" << h.user_int[i] << ", ";
		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++)
			if (h.user_float[i] != last_header.user_float[i])
				std::cout << "user_float[" << i << "]=" << h.user_float[i] << ", ";
	}
	// set header
	last_header = h;
	header_set = true;
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
	GDEBUG("Matrix size: %d, %d, %d\n", dimensionsIn_[0], dimensionsIn_[1], dimensionsIn_[2]);

	// get FOV
	field_of_view_[0] = r_space.fieldOfView_mm.x;
	field_of_view_[1] = e_space.fieldOfView_mm.y;
	field_of_view_[2] = e_space.fieldOfView_mm.z;
	GDEBUG("FOV: %f, %f, %f\n", field_of_view_[0], field_of_view_[1], field_of_view_[2]);

	// get echo line and echo partition
	iEchoLine_ = e_limits.kspace_encoding_step_1.get().center;
	iEchoPartition_ = e_limits.kspace_encoding_step_2.get().center;
	GDEBUG("echo line: %i, echo partition: %i\n", iEchoLine_, iEchoPartition_);

	// repetition time
	GlobalVar::instance()->fTR_ = h.sequenceParameters.get().TR.get().at(0);
	GDEBUG("TR: %f\n", GlobalVar::instance()->fTR_);
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
	GDEBUG("Matrix size: %d, %d, %d\n", dimensionsIn_[0], dimensionsIn_[1], dimensionsIn_[2]);

	// get FOV
	field_of_view_[0] = r_space.fieldOfView_mm().x();
	field_of_view_[1] = e_space.fieldOfView_mm().y();
	field_of_view_[2] = e_space.fieldOfView_mm().z();
	GDEBUG("FOV: %f, %f, %f\n", field_of_view_[0], field_of_view_[1], field_of_view_[2]);

	// get echo line and echo partition
	iEchoLine_ = e_limits.kspace_encoding_step_1().get().center();
	iEchoPartition_ = e_limits.kspace_encoding_step_2().get().center();
	GDEBUG("echo line: %i, echo partition: %i\n", iEchoLine_, iEchoPartition_);

	// repetition time
	GlobalVar::instance()->fTR_ = cfg->sequenceParameters().get().TR().at(0);
	GDEBUG("TR: %f\n", GlobalVar::instance()->fTR_);
#endif

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
			GWARN("No trajectory description present!");
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
	/*---------------------------------------------------*/
	/*----------- init buffer for k-space data ----------*/
	/*---------------------------------------------------*/
	if (!bufferkSpace_) {
		iNoChannels_ = m1->getObjectPtr()->active_channels;

		// initialize k-space buffer
		if (!(bufferkSpace_ = new hoNDArray<std::complex<float> >())) {
			GERROR("Failed to create k-space buffer\n");
			return GADGET_FAIL;
		}

		// get number of samples in acquisition (equals base resolution)
		iBaseRes_ = m1->getObjectPtr()->number_of_samples;

		GDEBUG("base res.: %d, no. scans: %lu, no. channel: %u\n", iBaseRes_, lNoScans_, iNoChannels_);

		// dimension vector of k-space array
		dimkSpace_.push_back(iBaseRes_);
		dimkSpace_.push_back(lNoScans_);
		dimkSpace_.push_back(iNoChannels_);

		// create buffer array for incoming k-space data (readout, time, channel)
		try {
			bufferkSpace_->create(&dimkSpace_);
		} catch (Gadgetron::bad_alloc) {
			print_not_enough_ram_msg<size_t>(dimkSpace_, sizeof(std::complex<float>));

			return GADGET_FAIL;
		} catch (std::runtime_error &err) {
			GEXCEPTION(err, "Failed to allocate k-space buffer array\n");
			return GADGET_FAIL;
		}

		// copy header information of first acquisition to global variable (create new header, set memory zero, copy header, push onto global vector)
		GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *tmp_m1 = new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();
		memset(tmp_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));
		fCopyAcqHeader(tmp_m1, m1);
		GlobalVar::instance()->AcqVec_.push_back(tmp_m1);
		GINFO("Receiving data..\n");
	}

	/*---------------------------------------------------*/
	/*--------- init buffer for navigator data ----------*/
	/*---------------------------------------------------*/
	if (!bufferNav_) {
		// initialize k-space buffer
		if (!(bufferNav_ = new hoNDArray<std::complex<float> >())) {
			GERROR("Failed to create navigator buffer\n");
			return GADGET_FAIL;
		}

		// dimension vector of navigator array
		dimNav_.push_back(iBaseRes_);
		dimNav_.push_back(lNoScans_);
		dimNav_.push_back(GlobalVar::instance()->iNavPERes_);
		dimNav_.push_back(iNoChannels_);
		iNoNav_ = 0;
		iNoNavLine_ = 0;

		GDEBUG("navigator dimensions: base res: %d, no. scans: %lu, PE resolution: %d, no. channels: %u\n", iBaseRes_, lNoScans_, GlobalVar::instance()->iNavPERes_, iNoChannels_);

		// create buffer array for incoming navigator data (readout, time, PE, channel)
		try {
			bufferNav_->create(&dimNav_);
		} catch (Gadgetron::bad_alloc) {
			print_not_enough_ram_msg<size_t>(dimNav_, sizeof(std::complex<float>));

			return GADGET_FAIL;
		} catch (std::runtime_error &err) {
			GEXCEPTION(err, "Failed to allocate navigator buffer array\n");
			return GADGET_FAIL;
		}

		GINFO("bufferNav_:\n");
		bufferNav_->print(std::cout);

		lCurrentScan_ = 0;
	}

	// protect Gadget from more inputs than expected
	if (lCurrentScan_ >= lNoScans_) {
		return GADGET_OK;
	}

	/*---------------------------------------------------*/
	/*----------- store incoming k-space data -----------*/
	/*---------------------------------------------------*/
	// navigator flag
	bool bNavigator = false;
	if (m1->getObjectPtr()->idx.set == 1) {		//m1->getObjectPtr()->user_int[1] & 0x1 != 0) {
		bNavigator = true;
	}

	// get current loop counters
	int samples		= m1->getObjectPtr()->number_of_samples;
	int line		= m1->getObjectPtr()->idx.kspace_encode_step_1;
	int partition	= m1->getObjectPtr()->idx.kspace_encode_step_2;

	// push current loop counters on according vector (temporal)
	GlobalVar::instance()->vPE_.push_back(line);
	GlobalVar::instance()->vPA_.push_back(partition);

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
		} else if (iNoNavLine_ == (GlobalVar::instance()->iNavPERes_/2)) {
			GlobalVar::instance()->vNavInd_.push_back(static_cast<float>(lCurrentScan_));
		}
	}

	lCurrentScan_++;

	/*---------------------------------------------------*/
	/*--------------- process sampled data --------------*/
	/*---------------------------------------------------*/
	if (lCurrentScan_ == lNoScans_) {
		GINFO("data received.. try to process data\n");

		// crop non-empty data from navigator array
		// check if last measurement is in between a navigator block
		if (iNoNavLine_ != 0) {
			iNoNav_--;
			GlobalVar::instance()->vNavInd_.pop_back();
		}

		GINFO("%i navigator data found..\n", iNoNav_);

		std::vector<size_t> vStart, vSize;

		vStart.push_back(0);
		vStart.push_back(0);
		vStart.push_back(0);
		vStart.push_back(0);

		vSize.push_back(bufferNav_->get_size(0));
		vSize.push_back(iNoNav_);
		vSize.push_back(bufferNav_->get_size(2));
		vSize.push_back(bufferNav_->get_size(3));

		bufferNav_->print(std::cout);
		get_subarray(*bufferNav_, vStart, vSize, *bufferNav_);
		bufferNav_->print(std::cout);

		GINFO("Flag - last in measurement detected..\n");
		GINFO("\n------------- navigator buffer ------------------\n");
		bufferNav_->print(std::cout);
		GINFO("\n------------- kSpace buffer -----------------\n");
		bufferkSpace_->print(std::cout);

		// create new ContainerMessages for header, navigator and kSpace data header
		GadgetContainerMessage<ISMRMRD::ImageHeader> *tmp_m1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();

		// initialize the image header
		memset(tmp_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));

		// initialize flags
		tmp_m1->getObjectPtr()->flags = 0;

		//tmp_m1->getObjectPtr()->user_int[0] = 7;
		tmp_m1->getObjectPtr()->user_int[0]			= iNPhases_;
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
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		tmp_m1->getObjectPtr()->data_type		= ISMRMRD::ISMRMRD_CXFLOAT;
#else
		tmp_m1->getObjectPtr()->image_data_type	= ISMRMRD::DATA_COMPLEX_FLOAT;
#endif
		tmp_m1->getObjectPtr()->image_index = static_cast<uint16_t>(++image_counter_);
		tmp_m1->getObjectPtr()->image_series_index = static_cast<uint16_t>(image_series_);

		// set user values, if Compressed Sensing is active
		//if(this->get_bool_value("CS_on") == true){
		tmp_m1->getObjectPtr()->user_float[0] = fCSAcc_;
		tmp_m1->getObjectPtr()->user_float[1] = fFullySa_/100;
		tmp_m1->getObjectPtr()->user_float[2] = fPartialFourierVal_;

		tmp_m1->getObjectPtr()->user_int[1] = iBodyRegion_;
		tmp_m1->getObjectPtr()->user_int[2] = iSamplingType_;
		tmp_m1->getObjectPtr()->user_int[3] = iVDMap_;
		tmp_m1->getObjectPtr()->user_int[4] = iESPReSSoDirection_;
		//}

		tmp_m1->getObjectPtr()->user_int[5] = iNoNav_;

		// navigator
		GadgetContainerMessage<hoNDArray<std::complex<float> > > *tmp_m2 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();

		// concatenate data with header
		tmp_m1->cont(tmp_m2);

		// create output
		try {
			tmp_m2->getObjectPtr()->create(bufferNav_->get_dimensions());
		} catch (std::runtime_error &err) {
			GEXCEPTION(err, "Unable to allocate new image array\n");

			tmp_m1->release();

			return -1;
		}

		// copy data
		memcpy(tmp_m2->getObjectPtr()->get_data_ptr(), bufferNav_->get_data_ptr(), sizeof(std::complex<float>)*bufferNav_->get_number_of_elements());

		// kSpace
		GadgetContainerMessage<hoNDArray<std::complex<float> > > *tmp_m3 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();

		// concatenate data
		tmp_m2->cont(tmp_m3);

		// create output
		try {
			tmp_m3->getObjectPtr()->create(bufferkSpace_->get_dimensions());
		} catch (std::runtime_error &err) {
			GEXCEPTION(err, "Unable to allocate new image array\n");

			tmp_m1->release();

			return -1;
		}

		// copy data
		memcpy(tmp_m3->getObjectPtr()->get_data_ptr(), bufferkSpace_->get_data_ptr(), sizeof(std::complex<float>)*bufferkSpace_->get_number_of_elements());

		// put on stream
		if (this->next()->putq(tmp_m1) < 0) {
			return GADGET_FAIL;
		}

		GDEBUG("global PE: %i, PA: %i\n", GlobalVar::instance()->vPE_.size(), GlobalVar::instance()->vPA_.size());

		return GADGET_OK;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_AccumulatorGadget)
