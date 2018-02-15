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
				fCopyAcqHeader(tmp_m1, GlobalVar::instance()->AcqVec_.at(par));

				// create empty 2D array
				GadgetContainerMessage< hoNDArray< float> >* sec_buffer_ = new GadgetContainerMessage<hoNDArray< float > >();
				try{sec_buffer_->getObjectPtr()->create(dimension[0], dimension[1], 1, 1, 1) ;}
				catch (std::runtime_error &err){
				  GEXCEPTION(err, "Unable to allocate new image array\n");

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

