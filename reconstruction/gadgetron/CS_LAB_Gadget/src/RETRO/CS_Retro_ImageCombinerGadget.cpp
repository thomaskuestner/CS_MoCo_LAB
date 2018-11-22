#include "CS_Retro_ImageCombinerGadget.h"

#include <mri_core_data.h>

using namespace Gadgetron;

// class constructor
CS_Retro_ImageCombinerGadget::CS_Retro_ImageCombinerGadget()
{
}

// class destructor
CS_Retro_ImageCombinerGadget::~CS_Retro_ImageCombinerGadget()
{
	delete data;
}

int CS_Retro_ImageCombinerGadget::process_config(ACE_Message_Block *mb)
{
	return GADGET_OK;
}

int CS_Retro_ImageCombinerGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	// receive [x y z resp_phases card_phases]
	hoNDArray<std::complex<float> > &received_data = *m2->getObjectPtr();

	// handle first initialization
	if (data == NULL) {
		number_of_respiratory_phases_	= get_number_of_gates(m1->getObjectPtr()->user_int[0], 0);
		number_of_cardiac_phases_		= get_number_of_gates(m1->getObjectPtr()->user_int[0], 1);

		// just pass if whole data is processed at once
		if (received_data.get_size(3) == number_of_respiratory_phases_ &&  received_data.get_size(4) == number_of_cardiac_phases_) {
			GINFO("Whole data was processed at once, so nothing has to be combined. Gadget is bypassed.\n");

			// put data on pipeline
			if (this->next()->putq(m1) < 0) {
				return GADGET_FAIL;
			}

			return GADGET_OK;
		}

		// order of array: [x y z resp_phases card_phases]
		data = new hoNDArray<std::complex<float> >(received_data.get_size(0), received_data.get_size(1), received_data.get_size(2), number_of_respiratory_phases_, number_of_cardiac_phases_);
	}

	// get current image position (which phase?)
	const unsigned int current_resp_phase = m1->getObjectPtr()->image_index;
	const unsigned int current_card_phase = m1->getObjectPtr()->image_series_index;

	// copy data to position
	size_t offset = current_resp_phase * data->get_size(0) * data->get_size(1) * data->get_size(2)
		+ current_card_phase * data->get_size(0) * data->get_size(1) * data->get_size(2) * data->get_size(3);
	memcpy(data->get_data_ptr()+offset, received_data.get_data_ptr(), received_data.get_number_of_bytes());

	// increase receive counter
	receive_counter_ += received_data.get_size(3)*received_data.get_size(4);

	// continue sending when finished
	if (receive_counter_ >= number_of_respiratory_phases_*number_of_cardiac_phases_) {
		// create new header
		GadgetContainerMessage<ISMRMRD::ImageHeader> *cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		fCopyImageHeader(cm1, m1->getObjectPtr());

		// create data element
		GadgetContainerMessage<hoNDArray<std::complex<float> > > *cm2 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();
		cm2->getObjectPtr()->create(data->get_dimensions());
		memcpy(cm2->getObjectPtr()->get_data_ptr(), data->get_data_ptr(), cm2->getObjectPtr()->get_number_of_bytes());

		// concatenate data
		cm1->cont(cm2);

		// put data on pipeline
		if (this->next()->putq(cm1) < 0) {
			return GADGET_FAIL;
		}
	}

	// free memory
	m1->release();

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_ImageCombinerGadget)
