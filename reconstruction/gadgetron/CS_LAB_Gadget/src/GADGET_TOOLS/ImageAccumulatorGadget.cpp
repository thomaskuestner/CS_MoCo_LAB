#include "ImageAccumulatorGadget.h"

using namespace Gadgetron;

ImageAccumulatorGadget::ImageAccumulatorGadget()
{
}

ImageAccumulatorGadget::~ImageAccumulatorGadget()
{
	// free hafBuffer_
	delete hafBuffer_;
	hafBuffer_ = NULL;
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
		vtDimensions_.push_back(get_number_of_gates(m1->getObjectPtr()->user_int[0], 0));
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
	fCopyImageHeader(tmp_m1, m1->getObjectPtr());

	// push header to global header vector - for kSpace
	GlobalVar::instance()->ImgHeadVec_.push_back(tmp_m1->getObjectPtr());

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

GADGET_FACTORY_DECLARE(ImageAccumulatorGadget)
