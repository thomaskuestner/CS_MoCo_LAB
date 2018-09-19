#include "CS_Retro_ImageFormatForSavingGadget.h"

#include <mri_core_data.h>

using namespace Gadgetron;

// class constructor
CS_Retro_ImageFormatForSavingGadget::CS_Retro_ImageFormatForSavingGadget()
{
}

// class destructor
CS_Retro_ImageFormatForSavingGadget::~CS_Retro_ImageFormatForSavingGadget()
{
}

int CS_Retro_ImageFormatForSavingGadget::process_config(ACE_Message_Block *mb)
{
	return GADGET_OK;
}

int CS_Retro_ImageFormatForSavingGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2)
{
	GDEBUG("Performing image reformatting and saving.\n");

	// set shortcut
	const hoNDArray<float> &data = *m2->getObjectPtr();

	// incoming image is [x y z resp_gates card_gates]
	// handle each card_gate separately
	for (size_t i = 0; i < data.get_size(4); i++) {
		// create new header
		GadgetContainerMessage<ISMRMRD::ImageHeader> *cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		fCopyImageHeader(cm1, m1->getObjectPtr());

		// reset channel value
		cm1->getObjectPtr()->channels = data.get_size(3);

		// now set image indices
		cm1->getObjectPtr()->image_index = i;
		cm1->getObjectPtr()->image_series_index = receive_counter_;

		// create data element
		GadgetContainerMessage<hoNDArray<float> > *cm2 = new GadgetContainerMessage<hoNDArray<float> >();
		cm2->getObjectPtr()->create(data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3));
		memcpy(cm2->getObjectPtr()->get_data_ptr(), data.get_data_ptr(), cm2->getObjectPtr()->get_number_of_bytes());

		// concatenate data
		cm1->cont(cm2);

		// put data on pipeline
		if (this->next()->putq(cm1) < 0) {
			return GADGET_FAIL;
		}
	}

	// increase receive_counter_
	receive_counter_++;

	// free memory
	m1->release();

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_ImageFormatForSavingGadget)
