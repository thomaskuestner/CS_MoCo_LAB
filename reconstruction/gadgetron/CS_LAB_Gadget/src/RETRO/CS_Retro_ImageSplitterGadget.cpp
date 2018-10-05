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
	// get gadget property
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	simultaneous_cardiac_phases_		= SimultaneousCardiacPhases.value();
	simultaneous_respiratory_phases_	= SimultaneousRespiratoryPhases.value();
#else
	simultaneous_cardiac_phases_		= *(get_int_value("SimultaneousCardiacPhases").get());
	simultaneous_respiratory_phases_	= *(get_int_value("SimultaneousRespiratoryPhases").get());
#endif

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

	// set parameters to all dimensions if zero
	if (simultaneous_cardiac_phases_ == 0) {
		simultaneous_cardiac_phases_ = data.get_size(5);
	}

	if (simultaneous_respiratory_phases_ == 0) {
		simultaneous_respiratory_phases_ = data.get_size(4);
	}

	// get loop counts (add +1 additionally if block won't match exactly (--> + X != 0))
	const size_t card_loop_amount = data.get_size(5)/simultaneous_cardiac_phases_ + (data.get_size(5) % simultaneous_cardiac_phases_ != 0);
	const size_t resp_loop_amount = data.get_size(4)/simultaneous_respiratory_phases_ + (data.get_size(4) % simultaneous_respiratory_phases_ != 0);

	for (size_t card_phase = 0; card_phase < card_loop_amount; card_phase++) {
		for (size_t resp_phase = 0; resp_phase < resp_loop_amount; resp_phase++) {
			// create image header
			GadgetContainerMessage<ISMRMRD::ImageHeader> *cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
			fCopyImageHeader(cm1, m1->getObjectPtr());

			// now set image indeces
			cm1->getObjectPtr()->image_index = resp_phase*simultaneous_respiratory_phases_;
			cm1->getObjectPtr()->image_series_index = card_phase*simultaneous_cardiac_phases_;

			// determine block size (if last iteration is reached: set block size to rest if there is a rest, otherwise keep standard length)
			const size_t card_block_size = card_phase == (card_loop_amount-1) ? (data.get_size(5) % simultaneous_cardiac_phases_ == 0 ? simultaneous_cardiac_phases_ : data.get_size(5) % simultaneous_cardiac_phases_) : simultaneous_cardiac_phases_;
			const size_t resp_block_size = resp_phase == (resp_loop_amount-1) ? (data.get_size(4) % simultaneous_respiratory_phases_ == 0 ? simultaneous_respiratory_phases_ : data.get_size(4) % simultaneous_respiratory_phases_) : simultaneous_respiratory_phases_;

			// create data element [kx ky kz c resp_phases card_phases]
			hoNDArray<std::complex<float> > data_slice;
			data_slice.create(data.get_size(0), data.get_size(1), data.get_size(2), data.get_size(3), resp_block_size, card_block_size);

			// copy elements
			for (size_t card_block_cnt = 0; card_block_cnt < card_block_size; card_block_cnt++) {
				const size_t origin_offset = card_phase*simultaneous_cardiac_phases_*data.get_size(0)*data.get_size(1)*data.get_size(2)*data.get_size(3)*data.get_size(4)
					+ resp_phase*simultaneous_respiratory_phases_*data.get_size(0)*data.get_size(1)*data.get_size(2)*data.get_size(3);
				const size_t copy_offset = card_block_cnt*data.get_size(0)*data.get_size(1)*data.get_size(2)*data.get_size(3)*data.get_size(4);
				const size_t target_offset = card_block_cnt*data_slice.get_size(0)*data_slice.get_size(1)*data_slice.get_size(2)*data_slice.get_size(3)*data_slice.get_size(4);
				memcpy(data_slice.get_data_ptr()+target_offset, data.get_data_ptr()+origin_offset+copy_offset, sizeof(std::complex<float>)*data_slice.get_size(0)*data_slice.get_size(1)*data_slice.get_size(2)*data_slice.get_size(3)*resp_block_size);
			}

			// re-permute data to [kx ky kz resp_phases card_phases c]
			std::vector<size_t> repermute_vector;
			repermute_vector.push_back(0);
			repermute_vector.push_back(1);
			repermute_vector.push_back(2);
			repermute_vector.push_back(4);
			repermute_vector.push_back(5);
			repermute_vector.push_back(3);
			data_slice = *permute(&data_slice, &repermute_vector, false);

			// create data message
			GadgetContainerMessage<hoNDArray<std::complex<float> > > *cm2 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();
			cm2->getObjectPtr()->create(data_slice.get_dimensions());
			memcpy(cm2->getObjectPtr()->get_data_ptr(), data.get_data_ptr(), cm2->getObjectPtr()->get_number_of_bytes());

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
