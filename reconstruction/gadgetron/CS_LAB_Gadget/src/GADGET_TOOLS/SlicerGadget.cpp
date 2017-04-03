/*	
file name	: 	SlicerGadget.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	implementation of the class "SlicerGadget.h"

references	:	-
*/

#include "SlicerGadget.h"

namespace Gadgetron{

// class constructor
SlicerGadget::SlicerGadget() : image_counter_(0), image_series_(0), bIs2D_(false), bIs3D_(false), bIs4D_(false){}
 
// class destructor 
SlicerGadget::~SlicerGadget(){}

// read ACE message block - flexible data header - (nothing to read)
int SlicerGadget::process_config(ACE_Message_Block* mb){return GADGET_OK;}

// copy the header information from header "m1" to header "tmp_m1"
int copy_header(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *tmp_m1, GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1){
	tmp_m1->getObjectPtr()->acquisition_time_stamp		= m1->getObjectPtr()->acquisition_time_stamp;
	tmp_m1->getObjectPtr()->active_channels				= m1->getObjectPtr()->active_channels;
	tmp_m1->getObjectPtr()->available_channels			= m1->getObjectPtr()->available_channels;
	tmp_m1->getObjectPtr()->center_sample				= m1->getObjectPtr()->center_sample;
	tmp_m1->getObjectPtr()->channel_mask[0]				= m1->getObjectPtr()->channel_mask[0];
	tmp_m1->getObjectPtr()->channel_mask[1]				= m1->getObjectPtr()->channel_mask[1];
	tmp_m1->getObjectPtr()->channel_mask[2]				= m1->getObjectPtr()->channel_mask[2];
	tmp_m1->getObjectPtr()->channel_mask[3]				= m1->getObjectPtr()->channel_mask[3];
	tmp_m1->getObjectPtr()->channel_mask[4]				= m1->getObjectPtr()->channel_mask[4];
	tmp_m1->getObjectPtr()->channel_mask[5]				= m1->getObjectPtr()->channel_mask[5];
	tmp_m1->getObjectPtr()->channel_mask[6]				= m1->getObjectPtr()->channel_mask[6];
	tmp_m1->getObjectPtr()->channel_mask[7]				= m1->getObjectPtr()->channel_mask[7];
	tmp_m1->getObjectPtr()->channel_mask[8]				= m1->getObjectPtr()->channel_mask[8];
	tmp_m1->getObjectPtr()->channel_mask[9]				= m1->getObjectPtr()->channel_mask[9];
	tmp_m1->getObjectPtr()->channel_mask[10]			= m1->getObjectPtr()->channel_mask[10];
	tmp_m1->getObjectPtr()->channel_mask[11]			= m1->getObjectPtr()->channel_mask[11];
	tmp_m1->getObjectPtr()->channel_mask[12]			= m1->getObjectPtr()->channel_mask[12];
	tmp_m1->getObjectPtr()->channel_mask[13]			= m1->getObjectPtr()->channel_mask[13];
	tmp_m1->getObjectPtr()->channel_mask[14]			= m1->getObjectPtr()->channel_mask[14];
	tmp_m1->getObjectPtr()->channel_mask[15]			= m1->getObjectPtr()->channel_mask[15];
	tmp_m1->getObjectPtr()->discard_post				= m1->getObjectPtr()->discard_post;
	tmp_m1->getObjectPtr()->discard_pre					= m1->getObjectPtr()->discard_pre;
	tmp_m1->getObjectPtr()->encoding_space_ref			= m1->getObjectPtr()->encoding_space_ref;
	tmp_m1->getObjectPtr()->flags						= m1->getObjectPtr()->flags;
	tmp_m1->getObjectPtr()->idx.average					= m1->getObjectPtr()->idx.average;
	tmp_m1->getObjectPtr()->idx.contrast				= m1->getObjectPtr()->idx.contrast;
	tmp_m1->getObjectPtr()->idx.kspace_encode_step_1	= m1->getObjectPtr()->idx.kspace_encode_step_1;
	tmp_m1->getObjectPtr()->idx.kspace_encode_step_2	= m1->getObjectPtr()->idx.kspace_encode_step_2;
	tmp_m1->getObjectPtr()->idx.phase					= m1->getObjectPtr()->idx.phase;
	tmp_m1->getObjectPtr()->idx.repetition				= m1->getObjectPtr()->idx.repetition;
	tmp_m1->getObjectPtr()->idx.segment					= m1->getObjectPtr()->idx.segment;
	tmp_m1->getObjectPtr()->idx.set						= m1->getObjectPtr()->idx.set;
	tmp_m1->getObjectPtr()->idx.slice					= m1->getObjectPtr()->idx.slice;
	tmp_m1->getObjectPtr()->idx.user[0]					= m1->getObjectPtr()->idx.user[0];
	tmp_m1->getObjectPtr()->idx.user[1]					= m1->getObjectPtr()->idx.user[1];
	tmp_m1->getObjectPtr()->idx.user[2]					= m1->getObjectPtr()->idx.user[2];
	tmp_m1->getObjectPtr()->idx.user[3]					= m1->getObjectPtr()->idx.user[3];
	tmp_m1->getObjectPtr()->idx.user[4]					= m1->getObjectPtr()->idx.user[4];
	tmp_m1->getObjectPtr()->idx.user[5]					= m1->getObjectPtr()->idx.user[5];
	tmp_m1->getObjectPtr()->idx.user[6]					= m1->getObjectPtr()->idx.user[6];
	tmp_m1->getObjectPtr()->idx.user[7]					= m1->getObjectPtr()->idx.user[7];
	tmp_m1->getObjectPtr()->measurement_uid				= m1->getObjectPtr()->measurement_uid;
	tmp_m1->getObjectPtr()->number_of_samples			= m1->getObjectPtr()->number_of_samples;
	tmp_m1->getObjectPtr()->patient_table_position[0]	= m1->getObjectPtr()->patient_table_position[0];
	tmp_m1->getObjectPtr()->patient_table_position[1]	= m1->getObjectPtr()->patient_table_position[1];
	tmp_m1->getObjectPtr()->patient_table_position[2]	= m1->getObjectPtr()->patient_table_position[2];
	tmp_m1->getObjectPtr()->phase_dir[0]				= m1->getObjectPtr()->phase_dir[0];
	tmp_m1->getObjectPtr()->phase_dir[1]				= m1->getObjectPtr()->phase_dir[1];
	tmp_m1->getObjectPtr()->phase_dir[2]				= m1->getObjectPtr()->phase_dir[2];
	tmp_m1->getObjectPtr()->physiology_time_stamp[0]	= m1->getObjectPtr()->physiology_time_stamp[0];
	tmp_m1->getObjectPtr()->physiology_time_stamp[1]	= m1->getObjectPtr()->physiology_time_stamp[1];
	tmp_m1->getObjectPtr()->physiology_time_stamp[2]	= m1->getObjectPtr()->physiology_time_stamp[2];
	tmp_m1->getObjectPtr()->physiology_time_stamp[3]	= m1->getObjectPtr()->physiology_time_stamp[3];
	tmp_m1->getObjectPtr()->physiology_time_stamp[4]	= m1->getObjectPtr()->physiology_time_stamp[4];
	tmp_m1->getObjectPtr()->physiology_time_stamp[5]	= m1->getObjectPtr()->physiology_time_stamp[5];
	tmp_m1->getObjectPtr()->physiology_time_stamp[6]	= m1->getObjectPtr()->physiology_time_stamp[6];
	tmp_m1->getObjectPtr()->physiology_time_stamp[7]	= m1->getObjectPtr()->physiology_time_stamp[7];
	tmp_m1->getObjectPtr()->position[0]					= m1->getObjectPtr()->position[0];
	tmp_m1->getObjectPtr()->position[1]					= m1->getObjectPtr()->position[1];
	tmp_m1->getObjectPtr()->position[2]					= m1->getObjectPtr()->position[2];
	tmp_m1->getObjectPtr()->read_dir[0]					= m1->getObjectPtr()->read_dir[0];
	tmp_m1->getObjectPtr()->read_dir[1]					= m1->getObjectPtr()->read_dir[1];
	tmp_m1->getObjectPtr()->read_dir[2]					= m1->getObjectPtr()->read_dir[2];
	tmp_m1->getObjectPtr()->sample_time_us				= m1->getObjectPtr()->sample_time_us;
	tmp_m1->getObjectPtr()->scan_counter				= m1->getObjectPtr()->scan_counter;
	tmp_m1->getObjectPtr()->slice_dir[0]				= m1->getObjectPtr()->slice_dir[0];
	tmp_m1->getObjectPtr()->slice_dir[1]				= m1->getObjectPtr()->slice_dir[1];
	tmp_m1->getObjectPtr()->slice_dir[2]				= m1->getObjectPtr()->slice_dir[2];
	tmp_m1->getObjectPtr()->trajectory_dimensions		= m1->getObjectPtr()->trajectory_dimensions;
	tmp_m1->getObjectPtr()->user_float[0]				= m1->getObjectPtr()->user_float[0];
	tmp_m1->getObjectPtr()->user_float[1]				= m1->getObjectPtr()->user_float[1];
	tmp_m1->getObjectPtr()->user_float[2]				= m1->getObjectPtr()->user_float[2];
	tmp_m1->getObjectPtr()->user_float[3]				= m1->getObjectPtr()->user_float[3];
	tmp_m1->getObjectPtr()->user_float[4]				= m1->getObjectPtr()->user_float[4];
	tmp_m1->getObjectPtr()->user_float[5]				= m1->getObjectPtr()->user_float[5];
	tmp_m1->getObjectPtr()->user_float[6]				= m1->getObjectPtr()->user_float[6];
	tmp_m1->getObjectPtr()->user_float[7]				= m1->getObjectPtr()->user_float[7];
	tmp_m1->getObjectPtr()->user_int[0]					= m1->getObjectPtr()->user_int[0];
	tmp_m1->getObjectPtr()->user_int[1]					= m1->getObjectPtr()->user_int[1];
	tmp_m1->getObjectPtr()->user_int[2]					= m1->getObjectPtr()->user_int[2];
	tmp_m1->getObjectPtr()->user_int[3]					= m1->getObjectPtr()->user_int[3];
	tmp_m1->getObjectPtr()->user_int[4]					= m1->getObjectPtr()->user_int[4];
	tmp_m1->getObjectPtr()->user_int[5]					= m1->getObjectPtr()->user_int[5];
	tmp_m1->getObjectPtr()->user_int[6]					= m1->getObjectPtr()->user_int[6];
	tmp_m1->getObjectPtr()->user_int[7]					= m1->getObjectPtr()->user_int[7];
	tmp_m1->getObjectPtr()->version						= m1->getObjectPtr()->version;

	return GADGET_OK;
}

// process(...): data processing - cuts the data into 2D data sets
int SlicerGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,GadgetContainerMessage< hoNDArray< float> >* m2)
{
	// check partitions/repetitions - more than one -> 3D/4D data - do "slicing"!!
	std::vector<size_t> dimension = *m2->getObjectPtr()->get_dimensions();
	size_t num_dims = m2->getObjectPtr()->get_number_of_dimensions();
	size_t num_rep = 0, num_par = 0;
	
	// get dimensions flag
	/*if (num_dims == 3){
		if ((dimension.at(0) > 1)&&(dimension.at(1) > 1)&&(dimension.at(2) == 1 )&&(m1->getObjectPtr()->user_int[0] == 1 || m1->getObjectPtr()->user_int[0] == 3))
			bIs2D_ = true;
		else if((dimension.at(0) > 1)&&(dimension.at(1) > 1)&&(dimension.at(2) > 1)&&(m1->getObjectPtr()->user_int[0] == 2 || m1->getObjectPtr()->user_int[0] == 5)){
			bIs3D_ = true;
			num_rep = 1;
			num_par = dimension.at(2);
		}
	}
	else if(num_dims == 4 &&(m1->getObjectPtr()->user_int[0] == 6)){
		bIs4D_ = true;
		num_rep = dimension.at(3);
		num_par = dimension.at(2);
	}*/
	if (num_dims == 2){
		bIs2D_ = true;
	}
	else if (num_dims == 3){
		bIs3D_ = true;
	}
	else if (num_dims == 4){
		bIs4D_ = true;
	}	

	if (bIs2D_){
		// data is 2D - do nothing
		if (this->next()->putq(m1) < 0) {
			return GADGET_FAIL;
		}
		return GADGET_OK;
	}
	else if(bIs3D_ || bIs4D_){
		// handle 3D/4D data - separate slices in own container and header

		// loop over repetitions
		for (int rep = 0; rep < num_rep; rep++){
		
			// loop over partitions
			for (int par = 0; par < num_par; par++){

				// new AcquisitionHeader
				GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* tmp_m1 = new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();

				// initialize the image header
				memset(tmp_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));

				// copy acquisition header (global variable) - needed for the correct slice position,..
				copy_header(tmp_m1, GlobalVar::instance()->AcqVec_.at(par));

				// create empty 2D array
				GadgetContainerMessage< hoNDArray< float> >* sec_buffer_ = new GadgetContainerMessage<hoNDArray< float > >();
				try{sec_buffer_->getObjectPtr()->create(dimension[0], dimension[1], 1, 1, 1) ;}
				catch (std::runtime_error &err){
#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				  GDEBUG("Unable to allocate new image array\n");
#else
				  GADGET_DEBUG_EXCEPTION(err,"Unable to allocate new image array\n");
#endif
				  sec_buffer_->release();
				  return -1;
				}

				// create header for 2D data			
				GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1_sec_buffer_ = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
				
				// initialize the image header
				memset(cm1_sec_buffer_->getObjectPtr(),0,sizeof(ISMRMRD::ImageHeader));
				cm1_sec_buffer_->getObjectPtr()->flags = 0;

				// connect header to data
				cm1_sec_buffer_->cont(sec_buffer_);

				// get the data pointers
				float* old_ptr = m2->getObjectPtr()->get_data_ptr();
				float* new_ptr = sec_buffer_->getObjectPtr()->get_data_ptr();

				// data offset - partition*Nx*Ny + rep*Nx*Ny*Nz
				size_t offset = par*dimension[0]*dimension[1] + rep*dimension[0]*dimension[1]*dimension[2];

				//copy the data for one slice
				memcpy(new_ptr, old_ptr+offset, sizeof(float)*(dimension[0]*dimension[1]));
				
				// set several header data (dimensions, FOV, image_index, number of active channels, slice and repetition index, data type, position, read-, phase-, slice-direction)
				cm1_sec_buffer_->getObjectPtr()->matrix_size[0] = (uint16_t) dimension[0];
				cm1_sec_buffer_->getObjectPtr()->matrix_size[1] = (uint16_t) dimension[1];
				cm1_sec_buffer_->getObjectPtr()->matrix_size[2] = 1;
				float* FoV = m1->getObjectPtr()->field_of_view;			
				cm1_sec_buffer_->getObjectPtr()->field_of_view[0]   =  FoV[0];
				cm1_sec_buffer_->getObjectPtr()->field_of_view[1]   =  FoV[1];
				cm1_sec_buffer_->getObjectPtr()->field_of_view[2]   =  FoV[2];
				cm1_sec_buffer_->getObjectPtr()->image_index = (uint16_t)(++image_counter_);
				cm1_sec_buffer_->getObjectPtr()->channels = 1;
				cm1_sec_buffer_->getObjectPtr()->slice = tmp_m1->getObjectPtr()->idx.slice;
				cm1_sec_buffer_->getObjectPtr()->repetition = tmp_m1->getObjectPtr()->idx.repetition;
#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				cm1_sec_buffer_->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;
				cm1_sec_buffer_->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
#else
				cm1_sec_buffer_->getObjectPtr()->image_data_type = ISMRMRD::DATA_FLOAT;
				cm1_sec_buffer_->getObjectPtr()->image_type = ISMRMRD::TYPE_MAGNITUDE;
#endif
				cm1_sec_buffer_->getObjectPtr()->image_series_index = (uint16_t)image_series_;
				cm1_sec_buffer_->getObjectPtr()->user_int[0] = tmp_m1->getObjectPtr()->idx.kspace_encode_step_2;
				
				// correct position vector (dependency to the slice orientation)
				if(tmp_m1->getObjectPtr()->slice_dir[0] != 0){
					cm1_sec_buffer_->getObjectPtr()->position[0] = tmp_m1->getObjectPtr()->position[0]-FoV[2]/2*(1-1/(float)num_par) + par*FoV[2]/num_par;
					cm1_sec_buffer_->getObjectPtr()->position[1] = tmp_m1->getObjectPtr()->position[1];
					cm1_sec_buffer_->getObjectPtr()->position[2] = tmp_m1->getObjectPtr()->position[2];
				}
				else if(tmp_m1->getObjectPtr()->slice_dir[1] != 0){
					cm1_sec_buffer_->getObjectPtr()->position[0] = tmp_m1->getObjectPtr()->position[0];
					cm1_sec_buffer_->getObjectPtr()->position[1] = tmp_m1->getObjectPtr()->position[1]-FoV[2]/2*(1-1/(float)num_par) + par*FoV[2]/num_par;
					cm1_sec_buffer_->getObjectPtr()->position[2] = tmp_m1->getObjectPtr()->position[2];
				}
				else if(tmp_m1->getObjectPtr()->slice_dir[2] != 0){
					cm1_sec_buffer_->getObjectPtr()->position[0] = tmp_m1->getObjectPtr()->position[0];
					cm1_sec_buffer_->getObjectPtr()->position[1] = tmp_m1->getObjectPtr()->position[1];
					cm1_sec_buffer_->getObjectPtr()->position[2] = tmp_m1->getObjectPtr()->position[2]-FoV[2]/2*(1-1/(float)num_par) + par*FoV[2]/num_par;
				}
				
				memcpy(cm1_sec_buffer_->getObjectPtr()->read_dir,tmp_m1->getObjectPtr()->read_dir,sizeof(float)*3);
				memcpy(cm1_sec_buffer_->getObjectPtr()->phase_dir,tmp_m1->getObjectPtr()->phase_dir,sizeof(float)*3);
				memcpy(cm1_sec_buffer_->getObjectPtr()->slice_dir,tmp_m1->getObjectPtr()->slice_dir, sizeof(float)*3);
				memcpy(cm1_sec_buffer_->getObjectPtr()->patient_table_position,tmp_m1->getObjectPtr()->patient_table_position, sizeof(float)*3);

				// concatenate
				if (this->next()->putq(cm1_sec_buffer_) < 0) {
    				return GADGET_FAIL;
				}
			} 
		}

		// clear AcquisitionHeader vector
		GlobalVar::instance()->AcqVec_.clear();

		m1->release();
		return GADGET_OK;
	}
}

GADGET_FACTORY_DECLARE(SlicerGadget)
}

