#include "ImageSaveHDFGadgetFLOAT.h"

namespace Gadgetron{

// save the incoming acquisition data to an hdf5 file
int ImageSaveHDFGadgetFLOAT::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray<float> >* m2)
{
	std::vector<size_t> dim = *m2->getObjectPtr()->get_dimensions();

	try{
		ismrmrd_dataset_->appendImageHeader(*m1->getObjectPtr(),"image.head");

		hoNDArray< float >* buffer_ = new hoNDArray< float >(dim, m2->getObjectPtr()->get_data_ptr(),false);
		
		std::vector<unsigned int> dims(dim.size());

		size_t i;
		for (i = 0; i< dim.size(); i++){
			dims[i] = dim[i];
		}

		if (ismrmrd_dataset_->appendArray(dims, m2->getObjectPtr()->get_data_ptr(), "image_0.img") < 0){
			GADGET_DEBUG1("Failed to write image data\n");
			return GADGET_FAIL;
		}
	}
	catch (...) {
        GADGET_DEBUG1("Error attempting to append images to HDF5 file\n");
        return GADGET_FAIL;
    }

	if (this->next()->putq(m1) == -1){
		m1->release();
		ACE_ERROR_RETURN( (LM_ERROR, ACE_TEXT("%p\n"), ACE_TEXT("ImageSaveHDFGadget::process, passing on data on to next gadget")), -1);
	}
	
	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(ImageSaveHDFGadgetFLOAT)
}
