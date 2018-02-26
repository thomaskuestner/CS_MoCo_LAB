#include "ImageAccumulatorGadget.h"

using namespace Gadgetron;

ImageAccumulatorGadget::ImageAccumulatorGadget()
{
}

ImageAccumulatorGadget::~ImageAccumulatorGadget()
{
	if (hafBuffer_) {
		delete hafBuffer_;
	}
}

// read flexible data header
int ImageAccumulatorGadget::process_config(ACE_Message_Block *mb)
{
  return GADGET_OK;
}

// buffers the incoming data to a fully buffered k-space array
int ImageAccumulatorGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2)
{
	// create temporal buffer if not already exists
	if (!hafBuffer_) {
		// get dimensions of slice
		vtDimensions_ = *m2->getObjectPtr()->get_dimensions();

		// get number of partitions
		vtDimensions_.push_back(m1->getObjectPtr()->matrix_size[2]);

		// extend dimensions by number of phases
		vtDimensions_.push_back(GlobalVar::instance()->iNPhases_);
		for (size_t iI = 0; iI < vtDimensions_.size(); iI++) {
			GDEBUG("image size - %i: %i\n", iI, vtDimensions_.at(iI));
		}
		
		// initialize buffer array for incoming data
		if (!(hafBuffer_ = new hoNDArray<float>())) {
			GERROR("Failed create buffer\n");

			return GADGET_FAIL;
		}

		// create buffer array for incoming data
		try {
			hafBuffer_->create(&vtDimensions_);
		} catch (std::runtime_error &err) {
			GEXCEPTION(err,"Failed allocate buffer array\n");
			return GADGET_FAIL;
		}

		GINFO("receiving data...\n");
	}
	
	// get pointers to the temporal buffer and the incoming data
	float* b = hafBuffer_->get_data_ptr();
	float* d = m2->getObjectPtr()->get_data_ptr();
	iPhs_ = static_cast<size_t>(m1->getObjectPtr()->phase);

	// get partition
	iPartition_++;
	iImageLoopCounter_++;

	GDEBUG("stack images: %i/%i - %i/%i - %i..\n", iPartition_+1, vtDimensions_[2], iPhs_+1, vtDimensions_[3], iImageLoopCounter_);

	// calculate the offset and copy the data into the temporal buffer
	size_t tOffset = iPartition_*vtDimensions_[1]*vtDimensions_[0] + iPhs_*vtDimensions_[1]*vtDimensions_[0]*vtDimensions_[2];
	memcpy(b+tOffset, d, sizeof(float)*vtDimensions_[1]*vtDimensions_[0]);

	// save header information of each slice (don't mess up position..)
	GadgetContainerMessage<ISMRMRD::ImageHeader> *tmp_m1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();

	// On some platforms, it is necessary to initialize the image header
	memset(tmp_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));
	
	// copy header data
	fCopyImageHeader(tmp_m1, m1);

	// push header to global header vector - for kSpace
	GlobalVar::instance()->ImgHeadVec_.push_back(tmp_m1);

	// reset partition number if end is reached
	if (iPartition_++ == vtDimensions_[2]) {
		iPartition_ = 0;
	}

	// end of gadget reached when all images received
	if (iImageLoopCounter_ == vtDimensions_[2]*vtDimensions_[3]) {
		GINFO("data received and stacked..\n");

		// copy data to new container
		GadgetContainerMessage<hoNDArray<float> > *cm2 = new GadgetContainerMessage<hoNDArray<float> >();

		// concatenate data with header
		tmp_m1->cont(cm2);

		try {
			cm2->getObjectPtr()->create(hafBuffer_->get_dimensions());
		} catch (std::runtime_error &err) {
			GEXCEPTION(err, "ImageAccumulatorGadget - Unable to allocate new image array\n");
			tmp_m1->release();
			return -1;
		}

		memcpy(cm2->getObjectPtr()->get_data_ptr(),b, sizeof(float)*hafBuffer_->get_number_of_elements());

		m1->release();

		// put on stream
		if (this->next()->putq(tmp_m1) < 0) {
    		return GADGET_FAIL;
		}

		return GADGET_OK;
	}

	return GADGET_OK;
}

//int ImageAccumulatorGadget::fCopyHeader(GadgetContainerMessage<ISMRMRD::ImageHeader> *tmp_m1, GadgetContainerMessage<ISMRMRD::ImageHeader>* m1) {
//	tmp_m1->getObjectPtr()->average						= m1->getObjectPtr()->average;
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

GADGET_FACTORY_DECLARE(ImageAccumulatorGadget)
