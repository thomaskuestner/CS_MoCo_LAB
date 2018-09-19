#include "CS_Retro_ImageSplitterGadget.h"

#include <mri_core_data.h>

using namespace Gadgetron;

// class constructor
CS_Retro_ImageSplitterGadget::CS_Retro_ImageSplitterGadget()
{
}

// class destructor
CS_Retro_ImageSplitterGadget::~CS_Retro_ImageSplitterGadget()
{
}

int CS_Retro_ImageSplitterGadget::process_config(ACE_Message_Block *mb)
{
	return GADGET_OK;
}

int CS_Retro_ImageSplitterGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	// input: kx-ky-kz-g1-g2-c
	hoNDArray<std::complex<float> > &data = *m2->getObjectPtr();

	// permute kSpace: kx-ky-kz-g1-g2-c -> kx-ky-kz-c-g1-g2
	std::vector<size_t> vtDimOrder;
	vtDimOrder.push_back(0);
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(5);
	vtDimOrder.push_back(3);
	vtDimOrder.push_back(4);
	data = *permute(&data, &vtDimOrder, false);

	for (size_t card_phase = 0; card_phase < data.get_size(5); card_phase++) {
		for (size_t resp_phase = 0; resp_phase < data.get_size(4); resp_phase++) {
			// create image header
			GadgetContainerMessage<ISMRMRD::ImageHeader> *cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
			fCopyImageHeader(cm1, m1->getObjectPtr());

			// now set image indeces
			cm1->getObjectPtr()->image_index = resp_phase;
			cm1->getObjectPtr()->image_series_index = card_phase;

			// create data element [kx ky kz c 1 1]
			GadgetContainerMessage<hoNDArray<std::complex<float> > > *cm2 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();
			cm2->getObjectPtr()->create(data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3), 1, 1);
			memcpy(cm2->getObjectPtr()->get_data_ptr(), data.get_data_ptr(), cm2->getObjectPtr()->get_number_of_bytes());

			// re-permute data to [kx ky kz 1 1 c]
			std::vector<size_t> repermute_vector;
			repermute_vector.push_back(0);
			repermute_vector.push_back(1);
			repermute_vector.push_back(2);
			repermute_vector.push_back(4);
			repermute_vector.push_back(5);
			repermute_vector.push_back(3);
			data = *permute(&data, &repermute_vector, false);

			// send data
			cm1->cont(cm2);
			if (this->next()->putq(cm1) < 0) {
				return GADGET_FAIL;
			}
		}
	}

	// free memory
	m1->release();

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_ImageSplitterGadget)
