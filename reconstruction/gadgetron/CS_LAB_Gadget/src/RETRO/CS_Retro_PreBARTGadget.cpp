#include "CS_Retro_PreBARTGadget.h"

#include <mri_core_data.h>

using namespace Gadgetron;

// class constructor
CS_Retro_PreBARTGadget::CS_Retro_PreBARTGadget()
{
}

// class destructor
CS_Retro_PreBARTGadget::~CS_Retro_PreBARTGadget()
{
}

int CS_Retro_PreBARTGadget::process_config(ACE_Message_Block *mb)
{
	return GADGET_OK;
}

int CS_Retro_PreBARTGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	// save image header
	GlobalVar::instance()->ImgHeadVec_.clear();

	// make copy of image header to be able to release m1 later
	ISMRMRD::ImageHeader *m1_cpy = new ISMRMRD::ImageHeader();
	memcpy(m1_cpy, m1->getObjectPtr(), sizeof(ISMRMRD::ImageHeader));

	// save copy of image header
	GlobalVar::instance()->ImgHeadVec_.push_back(m1_cpy);

	// get the pipeline content
	ISMRMRD::AcquisitionHeader header = *GlobalVar::instance()->AcqVec_.at(0);
	hoNDArray<std::complex<float> > data = *m2->getObjectPtr();

	// permute kSpace: kx-ky-kz-g1-g2-c -> kx-ky-kz-c-g1-g2
	std::vector<size_t> vtDimOrder;
	vtDimOrder.push_back(0);
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(5);
	vtDimOrder.push_back(3);
	vtDimOrder.push_back(4);
	data = *permute(&data, &vtDimOrder, false);

	// create hoNDArray with header
	hoNDArray<ISMRMRD::AcquisitionHeader> header_array;
	header_array.create(data.get_size(1), data.get_size(2), data.get_size(4), data.get_size(5));

	#pragma omp parallel for
	for (size_t i = 0; i < header_array.get_number_of_elements(); i++) {
		header_array.at(i) = header;
	}

	// pack content together
	Gadgetron::IsmrmrdDataBuffered buffered_data;
	buffered_data.data_ = data;
	buffered_data.headers_ = header_array;

	IsmrmrdReconBit *recon_bit = new IsmrmrdReconBit;
	recon_bit->data_ = buffered_data;

	// put new recon data on queue
	GadgetContainerMessage<IsmrmrdReconData> *recon_mc = new GadgetContainerMessage<IsmrmrdReconData>();
	recon_mc->getObjectPtr()->rbit_.push_back(*recon_bit);

	if (this->next()->putq(recon_mc) < 0) {
		return GADGET_FAIL;
	}

	// free memory
	m1->release();

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_PreBARTGadget)
