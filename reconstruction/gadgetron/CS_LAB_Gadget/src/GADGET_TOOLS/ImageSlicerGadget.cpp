#include "ImageSlicerGadget.h"

namespace Gadgetron{

// class constructor
ImageSlicerGadget::ImageSlicerGadget(){}
 
// class destructor 
ImageSlicerGadget::~ImageSlicerGadget(){}

// read ACE message block - flexible data header - (nothing to read)
int ImageSlicerGadget::process_config(ACE_Message_Block* mb){return GADGET_OK;}

// process(...): data processing - cuts the data into 2D data sets
int ImageSlicerGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,GadgetContainerMessage< hoNDArray< float> >* m2)
{	
	// loop over partitions
	for (int iPar = 0; iPar < (*m2->getObjectPtr()->get_dimensions())[2]; iPar++){
		
		// loop over phases
		for (int iPhs = 0; iPhs < (*m2->getObjectPtr()->get_dimensions())[3]; iPhs){

			GDEBUG("SlicerGadget: par: %i phs: %i\n", iPar, iPhs);
			// new image header
			GadgetContainerMessage<ISMRMRD::ImageHeader>* new_m1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();

			// init header
			memset(new_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));
		
			// copy header information from global stacked vector
			fCopyImageHeader(new_m1, GlobalVar::instance()->ImgHeadVec_.at(iPar));
			new_m1->getObjectPtr()->matrix_size[2] = 1;

			// create 2D array and fill 
			GadgetContainerMessage< hoNDArray<float> >* hafSecBuf = new GadgetContainerMessage< hoNDArray<float> >();
			try{hafSecBuf->getObjectPtr()->create((*m2->getObjectPtr()->get_dimensions())[0], (*m2->getObjectPtr()->get_dimensions())[1]) ;}
			catch (std::runtime_error &err){
				GEXCEPTION(err,"Unable to allocate new image array\n");
				hafSecBuf->release();
				return -1;
			}

			// concatenate
			new_m1->cont(hafSecBuf);

			// get current partition offset
			size_t tOffset = iPar * (*m2->getObjectPtr()->get_dimensions())[0]*(*m2->getObjectPtr()->get_dimensions())[1]+iPhs*(*m2->getObjectPtr()->get_dimensions())[0]*(*m2->getObjectPtr()->get_dimensions())[1]*(*m2->getObjectPtr()->get_dimensions())[2];
			float* old_ptr = m2->getObjectPtr()->get_data_ptr();
			float* new_ptr = hafSecBuf->getObjectPtr()->get_data_ptr();

			// copy data
			memcpy(new_ptr, old_ptr + tOffset, sizeof(float)*(*m2->getObjectPtr()->get_dimensions())[0]*(*m2->getObjectPtr()->get_dimensions())[1]);

			// put on Q
			if (this->next()->putq(new_m1) < 0) {
				return GADGET_FAIL;
			}
		}
	}

	// clear global variable
	GlobalVar::instance()->ImgHeadVec_.clear();

	m1->release();
	return GADGET_OK;
}

// copy header data
//int ImageSlicerGadget::fCopyHeader(GadgetContainerMessage<ISMRMRD::ImageHeader> *tmp_m1, GadgetContainerMessage<ISMRMRD::ImageHeader>* m1){
//tmp_m1->getObjectPtr()->average						= m1->getObjectPtr()->average;
//	tmp_m1->getObjectPtr()->channels					= m1->getObjectPtr()->channels;
//	tmp_m1->getObjectPtr()->contrast					= m1->getObjectPtr()->contrast;
//	tmp_m1->getObjectPtr()->field_of_view[0]			= m1->getObjectPtr()->field_of_view[0];
//	tmp_m1->getObjectPtr()->field_of_view[1]			= m1->getObjectPtr()->field_of_view[1];
//	tmp_m1->getObjectPtr()->field_of_view[2]			= m1->getObjectPtr()->field_of_view[2];
//	tmp_m1->getObjectPtr()->image_data_type				= m1->getObjectPtr()->image_data_type;
//	tmp_m1->getObjectPtr()->image_series_index			= m1->getObjectPtr()->image_series_index;
//	tmp_m1->getObjectPtr()->image_index					= m1->getObjectPtr()->image_index;
//	tmp_m1->getObjectPtr()->image_type					= m1->getObjectPtr()->image_type;
//	tmp_m1->getObjectPtr()->matrix_size[0]				= m1->getObjectPtr()->matrix_size[0];
//	tmp_m1->getObjectPtr()->matrix_size[1]				= m1->getObjectPtr()->matrix_size[1];
//	tmp_m1->getObjectPtr()->matrix_size[2]				= m1->getObjectPtr()->matrix_size[2];
//	tmp_m1->getObjectPtr()->phase						= m1->getObjectPtr()->phase;
//	tmp_m1->getObjectPtr()->repetition					= m1->getObjectPtr()->repetition;
//	tmp_m1->getObjectPtr()->set							= m1->getObjectPtr()->set;
//	tmp_m1->getObjectPtr()->slice						= m1->getObjectPtr()->slice;
//	tmp_m1->getObjectPtr()->acquisition_time_stamp		= m1->getObjectPtr()->acquisition_time_stamp;
//	tmp_m1->getObjectPtr()->flags						= m1->getObjectPtr()->flags;
//	tmp_m1->getObjectPtr()->measurement_uid				= m1->getObjectPtr()->measurement_uid;
//	tmp_m1->getObjectPtr()->patient_table_position[0]	= m1->getObjectPtr()->patient_table_position[0];
//	tmp_m1->getObjectPtr()->patient_table_position[1]	= m1->getObjectPtr()->patient_table_position[1];
//	tmp_m1->getObjectPtr()->patient_table_position[2]	= m1->getObjectPtr()->patient_table_position[2];
//	tmp_m1->getObjectPtr()->phase_dir[0]				= m1->getObjectPtr()->phase_dir[0];
//	tmp_m1->getObjectPtr()->phase_dir[1]				= m1->getObjectPtr()->phase_dir[1];
//	tmp_m1->getObjectPtr()->phase_dir[2]				= m1->getObjectPtr()->phase_dir[2];
//	tmp_m1->getObjectPtr()->physiology_time_stamp[0]	= m1->getObjectPtr()->physiology_time_stamp[0];
//	tmp_m1->getObjectPtr()->physiology_time_stamp[1]	= m1->getObjectPtr()->physiology_time_stamp[1];
//	tmp_m1->getObjectPtr()->physiology_time_stamp[2]	= m1->getObjectPtr()->physiology_time_stamp[2];
//	tmp_m1->getObjectPtr()->physiology_time_stamp[3]	= m1->getObjectPtr()->physiology_time_stamp[3];
//	tmp_m1->getObjectPtr()->physiology_time_stamp[4]	= m1->getObjectPtr()->physiology_time_stamp[4];
//	tmp_m1->getObjectPtr()->physiology_time_stamp[5]	= m1->getObjectPtr()->physiology_time_stamp[5];
//	tmp_m1->getObjectPtr()->physiology_time_stamp[6]	= m1->getObjectPtr()->physiology_time_stamp[6];
//	tmp_m1->getObjectPtr()->physiology_time_stamp[7]	= m1->getObjectPtr()->physiology_time_stamp[7];
//	tmp_m1->getObjectPtr()->position[0]					= m1->getObjectPtr()->position[0];
//	tmp_m1->getObjectPtr()->position[1]					= m1->getObjectPtr()->position[1];
//	tmp_m1->getObjectPtr()->position[2]					= m1->getObjectPtr()->position[2];
//	tmp_m1->getObjectPtr()->read_dir[0]					= m1->getObjectPtr()->read_dir[0];
//	tmp_m1->getObjectPtr()->read_dir[1]					= m1->getObjectPtr()->read_dir[1];
//	tmp_m1->getObjectPtr()->read_dir[2]					= m1->getObjectPtr()->read_dir[2];
//	tmp_m1->getObjectPtr()->slice_dir[0]				= m1->getObjectPtr()->slice_dir[0];
//	tmp_m1->getObjectPtr()->slice_dir[1]				= m1->getObjectPtr()->slice_dir[1];
//	tmp_m1->getObjectPtr()->slice_dir[2]				= m1->getObjectPtr()->slice_dir[2];
//	tmp_m1->getObjectPtr()->user_float[0]				= m1->getObjectPtr()->user_float[0];
//	tmp_m1->getObjectPtr()->user_float[1]				= m1->getObjectPtr()->user_float[1];
//	tmp_m1->getObjectPtr()->user_float[2]				= m1->getObjectPtr()->user_float[2];
//	tmp_m1->getObjectPtr()->user_float[3]				= m1->getObjectPtr()->user_float[3];
//	tmp_m1->getObjectPtr()->user_float[4]				= m1->getObjectPtr()->user_float[4];
//	tmp_m1->getObjectPtr()->user_float[5]				= m1->getObjectPtr()->user_float[5];
//	tmp_m1->getObjectPtr()->user_float[6]				= m1->getObjectPtr()->user_float[6];
//	tmp_m1->getObjectPtr()->user_float[7]				= m1->getObjectPtr()->user_float[7];
//	tmp_m1->getObjectPtr()->user_int[0]					= m1->getObjectPtr()->user_int[0];
//	tmp_m1->getObjectPtr()->user_int[1]					= m1->getObjectPtr()->user_int[1];
//	tmp_m1->getObjectPtr()->user_int[2]					= m1->getObjectPtr()->user_int[2];
//	tmp_m1->getObjectPtr()->user_int[3]					= m1->getObjectPtr()->user_int[3];
//	tmp_m1->getObjectPtr()->user_int[4]					= m1->getObjectPtr()->user_int[4];
//	tmp_m1->getObjectPtr()->user_int[5]					= m1->getObjectPtr()->user_int[5];
//	tmp_m1->getObjectPtr()->user_int[6]					= m1->getObjectPtr()->user_int[6];
//	tmp_m1->getObjectPtr()->user_int[7]					= m1->getObjectPtr()->user_int[7];
//	tmp_m1->getObjectPtr()->version						= m1->getObjectPtr()->version;
//	
//	return GADGET_OK;
//}

GADGET_FACTORY_DECLARE(ImageSlicerGadget)
}

