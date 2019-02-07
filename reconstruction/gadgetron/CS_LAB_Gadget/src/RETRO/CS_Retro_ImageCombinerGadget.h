#ifndef CS_RETRO_IMAGECOMBINERGADGET_H
#define CS_RETRO_IMAGECOMBINERGADGET_H

#pragma once
#include "Gadget.h"
#include "hoNDArray.h"
#include "CS_LAB_export.h"
#include <ismrmrd.h>
#include <complex>
#include "GadgetIsmrmrdReadWrite.h"
#include "GlobalVar.h"
#include "SomeFunctions.h"

namespace Gadgetron
{
	class EXPORTCSLAB CS_Retro_ImageCombinerGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray<std::complex<float> > >
	{
	public:
		CS_Retro_ImageCombinerGadget();
		~CS_Retro_ImageCombinerGadget();

		GADGET_DECLARE(CS_Retro_ImageCombinerGadget);

		int close(unsigned long flags) override;

	protected:
		int process_config(ACE_Message_Block *mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2);

	private:
		hoNDArray<std::complex<float> > *data_ = NULL;
		GadgetContainerMessage<ISMRMRD::ImageHeader> *return_message_ = NULL;
		unsigned int receive_counter_ = 0;
		unsigned int number_of_cardiac_phases_, number_of_respiratory_phases_;
	};
}

#endif //CS_RETRO_IMAGECOMBINERGADGET_H
