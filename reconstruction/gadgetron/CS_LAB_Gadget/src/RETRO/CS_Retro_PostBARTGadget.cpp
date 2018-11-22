#include "CS_Retro_PostBARTGadget.h"

#include <mri_core_data.h>

using namespace Gadgetron;

// class constructor
CS_Retro_PostBARTGadget::CS_Retro_PostBARTGadget()
{
}

// class destructor
CS_Retro_PostBARTGadget::~CS_Retro_PostBARTGadget()
{
}

int CS_Retro_PostBARTGadget::process_config(ACE_Message_Block *mb)
{
	return GADGET_OK;
}

int CS_Retro_PostBARTGadget::process(GadgetContainerMessage<IsmrmrdImageArray> *m1)
{
	// get data
	hoNDArray<std::complex<float> > data = m1->getObjectPtr()->data_;

	// save image indices for later usage
	const unsigned int resp_phase = m1->getObjectPtr()->headers_.at(0).user_int[0];
	const unsigned int card_phase = m1->getObjectPtr()->headers_.at(0).user_int[1];

	// free memory
	m1->release();

	// permute array: kx-ky-kz-c-t-s-slc -> kx-ky-kz-t-s-c-slc (with c=slc=1 - will be cropped later on)
	std::vector<size_t> vtDimOrder;
	vtDimOrder.push_back(0);
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(5);
	vtDimOrder.push_back(6);
	vtDimOrder.push_back(4);
	vtDimOrder.push_back(3);
	data = *permute(&data, &vtDimOrder, false);
	
	// crop array and convert to GadgetContainerMessage
	GadgetContainerMessage<hoNDArray<std::complex<float> > > *cm2 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();
	cm2->getObjectPtr()->create(data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3), data.get_size(4));
	memcpy(cm2->getObjectPtr()->get_data_ptr(), data.get_data_ptr(), cm2->getObjectPtr()->get_number_of_bytes());

	// restore image header from pre BART state
	GadgetContainerMessage<ISMRMRD::ImageHeader> *cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	fCopyImageHeader(cm1, GlobalVar::instance()->ImgHeadVec_.at(0));

	// free memory in GlobalVar
	delete GlobalVar::instance()->ImgHeadVec_.at(0);
	GlobalVar::instance()->ImgHeadVec_.erase(GlobalVar::instance()->ImgHeadVec_.begin(), GlobalVar::instance()->ImgHeadVec_.begin()+1);

	// restore image indices
	cm1->getObjectPtr()->image_index = resp_phase;
	cm1->getObjectPtr()->image_series_index = card_phase;

	// concatenate data
	cm1->cont(cm2);

	// put data on pipeline
	if (this->next()->putq(cm1) < 0) {
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_PostBARTGadget)
