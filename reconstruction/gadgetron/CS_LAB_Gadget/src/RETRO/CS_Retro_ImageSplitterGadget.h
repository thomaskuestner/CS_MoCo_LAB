#ifndef CS_RETRO_IMAGESPLITTERGADGET_H
#define CS_RETRO_IMAGESPLITTERGADGET_H

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
	class EXPORTCSLAB CS_Retro_ImageSplitterGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray<std::complex<float> > >
	{
	public:
		CS_Retro_ImageSplitterGadget();
		~CS_Retro_ImageSplitterGadget();

		GADGET_DECLARE(CS_Retro_ImageSplitterGadget);

	protected:
		int process_config(ACE_Message_Block *mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2);
	};
}

#endif //CS_RETRO_IMAGESPLITTERGADGET_H
