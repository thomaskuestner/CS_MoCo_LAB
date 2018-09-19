#ifndef CS_RETRO_IMAGEFORMATFORSAVINGGADGET_H
#define CS_RETRO_IMAGEFORMATFORSAVINGGADGET_H

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
	class EXPORTCSLAB CS_Retro_ImageFormatForSavingGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray<float> >
	{
	public:
		CS_Retro_ImageFormatForSavingGadget();
		~CS_Retro_ImageFormatForSavingGadget();

		GADGET_DECLARE(CS_Retro_ImageFormatForSavingGadget);

	protected:
		int process_config(ACE_Message_Block *mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2);

	private:
		unsigned int receive_counter_ = 0;
	};
}

#endif //CS_RETRO_IMAGEFORMATFORSAVINGGADGET_H
