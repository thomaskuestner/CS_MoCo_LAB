#include "ImageSlicerGadget.h"

using namespace Gadgetron;

// class constructor
ImageSlicerGadget::ImageSlicerGadget()
{
}

// class destructor 
ImageSlicerGadget::~ImageSlicerGadget()
{
}

int ImageSlicerGadget::process_config(ACE_Message_Block *mb)
{
	return GADGET_OK;
}

// data processing - cuts the data into 2D data sets
int ImageSlicerGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1,GadgetContainerMessage<hoNDArray<float> > *m2)
{
	// loop over partitions
	for (size_t iPar = 0; iPar < (*m2->getObjectPtr()->get_dimensions())[2]; iPar++) {
		// loop over phases
		for (size_t iPhs = 0; iPhs < (*m2->getObjectPtr()->get_dimensions())[3]; iPhs++) {
			GDEBUG("SlicerGadget: par: %i phs: %i\n", iPar, iPhs);

			// new image header
			GadgetContainerMessage<ISMRMRD::ImageHeader>* new_m1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();

			// init header
			memset(new_m1->getObjectPtr(), 0, sizeof(ISMRMRD::ImageHeader));
		
			// copy header information from global stacked vector
			fCopyImageHeader(new_m1, GlobalVar::instance()->ImgHeadVec_.at(iPar));
			new_m1->getObjectPtr()->matrix_size[2] = 1;

			// create 2D array and fill
			GadgetContainerMessage<hoNDArray<float> > *hafSecBuf = new GadgetContainerMessage<hoNDArray<float> >();
			try {
				hafSecBuf->getObjectPtr()->create((*m2->getObjectPtr()->get_dimensions())[0], (*m2->getObjectPtr()->get_dimensions())[1]);
			} catch (std::runtime_error &err) {
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

			// put on Queue
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

GADGET_FACTORY_DECLARE(ImageSlicerGadget)
