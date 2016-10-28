#include "CS_MRIImageToAcquisitionGadget.h"

using namespace Gadgetron;

// convert the image data set to individual acquisitions
int CS_MRIImageToAcquisitionGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
	// get dimensions of incoming data (x,y,z,t,c)
	vDims_ = *m2->getObjectPtr()->get_dimensions();
	
	// ---------------------------------------------------------------------------------
	// ------------------------------ conversion ---------------------------------------
	// ---------------------------------------------------------------------------------
	for (int iPhase = 0; iPhase < vDims_[3]; iPhase++)
			for (int iPartition = 0; iPartition < vDims_[2]; iPartition++)
				for (int iLine = 0; iLine < vDims_[1]; iLine++){

					// new AcquisitionHeader
					GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* GC_acq_hdr_m1 = new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();

					// init acquisition header
					memset(GC_acq_hdr_m1->getObjectPtr(), 0, sizeof(ISMRMRD::AcquisitionHeader));

					// copy data from global stored header
					fCopyHeader(GC_acq_hdr_m1, CS_GlobalVar::instance()->AcqVec_.at(0));

					// correct some header information
					fCorrectHeader(GC_acq_hdr_m1, iLine, iPartition, iPhase);

					// create empty 1D array (COL*CHA)
					GadgetContainerMessage< hoNDArray< std::complex<float> > >* hacfTmp = new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
					try{hacfTmp->getObjectPtr()->create(vDims_[0]*m1->getObjectPtr()->channels) ;}
					catch (std::runtime_error &err){
						GADGET_DEBUG_EXCEPTION(err,"Unable to allocate new image array\n");
						hacfTmp->release();
						return -1;
					}

					// fill 1D array
					std::complex<float> *cfPtrOut = hacfTmp->getObjectPtr()->get_data_ptr(), *cfPtrIn = m2->getObjectPtr()->get_data_ptr();
					for (int iC = 0; iC < m1->getObjectPtr()->channels; iC++){
						int iOffsetOut	= iC*vDims_[0];
						int iOffsetIn	= iLine*vDims_[0] + iPartition*vDims_[0]*vDims_[1] + iPhase*vDims_[0]*vDims_[1]*vDims_[2] + iC*vDims_[0]*vDims_[1]*vDims_[2]*vDims_[3];
						memcpy(cfPtrOut + iOffsetOut, cfPtrIn + iOffsetIn, sizeof(std::complex<float>)*vDims_[0]*m1->getObjectPtr()->channels);
					}

					// concat header to data
					GC_acq_hdr_m1->cont(hacfTmp);
			
					
					// concatenate
					if (this->next()->putq(GC_acq_hdr_m1) < 0)
    					return GADGET_FAIL;
				}

	// clear AcquisitionHeader vector
	CS_GlobalVar::instance()->AcqVec_.clear();

	return GADGET_OK;
}

// correct header information for filled full kSpace data
int CS_MRIImageToAcquisitionGadget::fCorrectHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_hdr_m1, int iLine, int iPartition, int iPhase){
	// fill temporal acquisition header with values from first scan
	fCopyHeader(GC_acq_hdr_m1, CS_GlobalVar::instance()->AcqVec_.at(0));

	// set loop counter
	GC_acq_hdr_m1->getObjectPtr()->idx.kspace_encode_step_1 = iLine;
	GC_acq_hdr_m1->getObjectPtr()->idx.kspace_encode_step_2 = iPartition;
	GC_acq_hdr_m1->getObjectPtr()->idx.phase				= iPhase;
	
	// correct scan counter
	GC_acq_hdr_m1->getObjectPtr()->scan_counter = iLine + iPartition*vDims_[1] + iPhase*vDims_[1]*vDims_[2];

	// clear false flags
	GC_acq_hdr_m1->getObjectPtr()->flags = 0;

	// ---------------------------------------------------------------------------------
	// -------------------------- set several flags ------------------------------------
	// ---------------------------------------------------------------------------------
	// first scan in partition encoding
	if (iPartition == 0)
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_FIRST_IN_ENCODE_STEP2-1);

	// last scan in partition encoding
	if (iPartition == vDims_[2]-1)
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_LAST_IN_ENCODE_STEP2-1);
	
	// first scan in phase encoding
	if (iLine == 0)
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_FIRST_IN_ENCODE_STEP1-1);

	// last scan in phase encoding
	if (iLine == vDims_[1]-1)
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_LAST_IN_ENCODE_STEP1-1);

	// last scan in measurement
	if (iLine == vDims_[1]-1 && iPartition == vDims_[2]-1 && iPhase == vDims_[3]-1)
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_LAST_IN_MEASUREMENT-1);

	return GADGET_OK;
}

// copy all AcquisitionHeader values
int CS_MRIImageToAcquisitionGadget::fCopyHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1, GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1_new){
	GC_acq_m1_new->getObjectPtr()->acquisition_time_stamp		= GC_acq_m1->getObjectPtr()->acquisition_time_stamp;
	GC_acq_m1_new->getObjectPtr()->active_channels				= GC_acq_m1->getObjectPtr()->active_channels;
	GC_acq_m1_new->getObjectPtr()->available_channels			= GC_acq_m1->getObjectPtr()->available_channels;
	GC_acq_m1_new->getObjectPtr()->center_sample				= GC_acq_m1->getObjectPtr()->center_sample;
	GC_acq_m1_new->getObjectPtr()->channel_mask[0]				= GC_acq_m1->getObjectPtr()->channel_mask[0];
	GC_acq_m1_new->getObjectPtr()->channel_mask[1]				= GC_acq_m1->getObjectPtr()->channel_mask[1];
	GC_acq_m1_new->getObjectPtr()->channel_mask[2]				= GC_acq_m1->getObjectPtr()->channel_mask[2];
	GC_acq_m1_new->getObjectPtr()->channel_mask[3]				= GC_acq_m1->getObjectPtr()->channel_mask[3];
	GC_acq_m1_new->getObjectPtr()->channel_mask[4]				= GC_acq_m1->getObjectPtr()->channel_mask[4];
	GC_acq_m1_new->getObjectPtr()->channel_mask[5]				= GC_acq_m1->getObjectPtr()->channel_mask[5];
	GC_acq_m1_new->getObjectPtr()->channel_mask[6]				= GC_acq_m1->getObjectPtr()->channel_mask[6];
	GC_acq_m1_new->getObjectPtr()->channel_mask[7]				= GC_acq_m1->getObjectPtr()->channel_mask[7];
	GC_acq_m1_new->getObjectPtr()->channel_mask[8]				= GC_acq_m1->getObjectPtr()->channel_mask[8];
	GC_acq_m1_new->getObjectPtr()->channel_mask[9]				= GC_acq_m1->getObjectPtr()->channel_mask[9];
	GC_acq_m1_new->getObjectPtr()->channel_mask[10]				= GC_acq_m1->getObjectPtr()->channel_mask[10];
	GC_acq_m1_new->getObjectPtr()->channel_mask[11]				= GC_acq_m1->getObjectPtr()->channel_mask[11];
	GC_acq_m1_new->getObjectPtr()->channel_mask[12]				= GC_acq_m1->getObjectPtr()->channel_mask[12];
	GC_acq_m1_new->getObjectPtr()->channel_mask[13]				= GC_acq_m1->getObjectPtr()->channel_mask[13];
	GC_acq_m1_new->getObjectPtr()->channel_mask[14]				= GC_acq_m1->getObjectPtr()->channel_mask[14];
	GC_acq_m1_new->getObjectPtr()->channel_mask[15]				= GC_acq_m1->getObjectPtr()->channel_mask[15];
	GC_acq_m1_new->getObjectPtr()->discard_post					= GC_acq_m1->getObjectPtr()->discard_post;
	GC_acq_m1_new->getObjectPtr()->discard_pre					= GC_acq_m1->getObjectPtr()->discard_pre;
	GC_acq_m1_new->getObjectPtr()->encoding_space_ref			= GC_acq_m1->getObjectPtr()->encoding_space_ref;
	GC_acq_m1_new->getObjectPtr()->flags						= GC_acq_m1->getObjectPtr()->flags;
	GC_acq_m1_new->getObjectPtr()->idx.average					= GC_acq_m1->getObjectPtr()->idx.average;
	GC_acq_m1_new->getObjectPtr()->idx.contrast					= GC_acq_m1->getObjectPtr()->idx.contrast;
	GC_acq_m1_new->getObjectPtr()->idx.kspace_encode_step_1		= GC_acq_m1->getObjectPtr()->idx.kspace_encode_step_1;
	GC_acq_m1_new->getObjectPtr()->idx.kspace_encode_step_2		= GC_acq_m1->getObjectPtr()->idx.kspace_encode_step_2;
	GC_acq_m1_new->getObjectPtr()->idx.phase					= GC_acq_m1->getObjectPtr()->idx.phase;
	GC_acq_m1_new->getObjectPtr()->idx.repetition				= GC_acq_m1->getObjectPtr()->idx.repetition;
	GC_acq_m1_new->getObjectPtr()->idx.segment					= GC_acq_m1->getObjectPtr()->idx.segment;
	GC_acq_m1_new->getObjectPtr()->idx.set						= GC_acq_m1->getObjectPtr()->idx.set;
	GC_acq_m1_new->getObjectPtr()->idx.slice					= GC_acq_m1->getObjectPtr()->idx.slice;
	GC_acq_m1_new->getObjectPtr()->idx.user[0]					= GC_acq_m1->getObjectPtr()->idx.user[0];
	GC_acq_m1_new->getObjectPtr()->idx.user[1]					= GC_acq_m1->getObjectPtr()->idx.user[1];
	GC_acq_m1_new->getObjectPtr()->idx.user[2]					= GC_acq_m1->getObjectPtr()->idx.user[2];
	GC_acq_m1_new->getObjectPtr()->idx.user[3]					= GC_acq_m1->getObjectPtr()->idx.user[3];
	GC_acq_m1_new->getObjectPtr()->idx.user[4]					= GC_acq_m1->getObjectPtr()->idx.user[4];
	GC_acq_m1_new->getObjectPtr()->idx.user[5]					= GC_acq_m1->getObjectPtr()->idx.user[5];
	GC_acq_m1_new->getObjectPtr()->idx.user[6]					= GC_acq_m1->getObjectPtr()->idx.user[6];
	GC_acq_m1_new->getObjectPtr()->idx.user[7]					= GC_acq_m1->getObjectPtr()->idx.user[7];
	GC_acq_m1_new->getObjectPtr()->measurement_uid				= GC_acq_m1->getObjectPtr()->measurement_uid;
	GC_acq_m1_new->getObjectPtr()->number_of_samples			= GC_acq_m1->getObjectPtr()->number_of_samples;
	GC_acq_m1_new->getObjectPtr()->patient_table_position[0]	= GC_acq_m1->getObjectPtr()->patient_table_position[0];
	GC_acq_m1_new->getObjectPtr()->patient_table_position[1]	= GC_acq_m1->getObjectPtr()->patient_table_position[1];
	GC_acq_m1_new->getObjectPtr()->patient_table_position[2]	= GC_acq_m1->getObjectPtr()->patient_table_position[2];
	GC_acq_m1_new->getObjectPtr()->phase_dir[0]					= GC_acq_m1->getObjectPtr()->phase_dir[0];
	GC_acq_m1_new->getObjectPtr()->phase_dir[1]					= GC_acq_m1->getObjectPtr()->phase_dir[1];
	GC_acq_m1_new->getObjectPtr()->phase_dir[2]					= GC_acq_m1->getObjectPtr()->phase_dir[2];
	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[0]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[0];
	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[1]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[1];
	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[2]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[2];
	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[3]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[3];
	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[4]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[4];
	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[5]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[5];
	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[6]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[6];
	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[7]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[7];
	GC_acq_m1_new->getObjectPtr()->position[0]					= GC_acq_m1->getObjectPtr()->position[0];
	GC_acq_m1_new->getObjectPtr()->position[1]					= GC_acq_m1->getObjectPtr()->position[1];
	GC_acq_m1_new->getObjectPtr()->position[2]					= GC_acq_m1->getObjectPtr()->position[2];
	GC_acq_m1_new->getObjectPtr()->read_dir[0]					= GC_acq_m1->getObjectPtr()->read_dir[0];
	GC_acq_m1_new->getObjectPtr()->read_dir[1]					= GC_acq_m1->getObjectPtr()->read_dir[1];
	GC_acq_m1_new->getObjectPtr()->read_dir[2]					= GC_acq_m1->getObjectPtr()->read_dir[2];
	GC_acq_m1_new->getObjectPtr()->sample_time_us				= GC_acq_m1->getObjectPtr()->sample_time_us;
	GC_acq_m1_new->getObjectPtr()->scan_counter					= GC_acq_m1->getObjectPtr()->scan_counter;
	GC_acq_m1_new->getObjectPtr()->slice_dir[0]					= GC_acq_m1->getObjectPtr()->slice_dir[0];
	GC_acq_m1_new->getObjectPtr()->slice_dir[1]					= GC_acq_m1->getObjectPtr()->slice_dir[1];
	GC_acq_m1_new->getObjectPtr()->slice_dir[2]					= GC_acq_m1->getObjectPtr()->slice_dir[2];
	GC_acq_m1_new->getObjectPtr()->trajectory_dimensions		= GC_acq_m1->getObjectPtr()->trajectory_dimensions;
	GC_acq_m1_new->getObjectPtr()->user_float[0]				= GC_acq_m1->getObjectPtr()->user_float[0];
	GC_acq_m1_new->getObjectPtr()->user_float[1]				= GC_acq_m1->getObjectPtr()->user_float[1];
	GC_acq_m1_new->getObjectPtr()->user_float[2]				= GC_acq_m1->getObjectPtr()->user_float[2];
	GC_acq_m1_new->getObjectPtr()->user_float[3]				= GC_acq_m1->getObjectPtr()->user_float[3];
	GC_acq_m1_new->getObjectPtr()->user_float[4]				= GC_acq_m1->getObjectPtr()->user_float[4];
	GC_acq_m1_new->getObjectPtr()->user_float[5]				= GC_acq_m1->getObjectPtr()->user_float[5];
	GC_acq_m1_new->getObjectPtr()->user_float[6]				= GC_acq_m1->getObjectPtr()->user_float[6];
	GC_acq_m1_new->getObjectPtr()->user_float[7]				= GC_acq_m1->getObjectPtr()->user_float[7];
	GC_acq_m1_new->getObjectPtr()->user_int[0]					= GC_acq_m1->getObjectPtr()->user_int[0];
	GC_acq_m1_new->getObjectPtr()->user_int[1]					= GC_acq_m1->getObjectPtr()->user_int[1];
	GC_acq_m1_new->getObjectPtr()->user_int[2]					= GC_acq_m1->getObjectPtr()->user_int[2];
	GC_acq_m1_new->getObjectPtr()->user_int[3]					= GC_acq_m1->getObjectPtr()->user_int[3];
	GC_acq_m1_new->getObjectPtr()->user_int[4]					= GC_acq_m1->getObjectPtr()->user_int[4];
	GC_acq_m1_new->getObjectPtr()->user_int[5]					= GC_acq_m1->getObjectPtr()->user_int[5];
	GC_acq_m1_new->getObjectPtr()->user_int[6]					= GC_acq_m1->getObjectPtr()->user_int[6];
	GC_acq_m1_new->getObjectPtr()->user_int[7]					= GC_acq_m1->getObjectPtr()->user_int[7];
	GC_acq_m1_new->getObjectPtr()->version						= GC_acq_m1->getObjectPtr()->version;
	
	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_MRIImageToAcquisitionGadget)
