#ifndef CS_RETRO_POSTBARTGADGET_H
#define CS_RETRO_POSTBARTGADGET_H

#pragma once
#include "Gadget.h"
#include "hoNDArray.h"
#include "CS_LAB_export.h"
#include <ismrmrd.h>
#include <complex>
#include "GadgetIsmrmrdReadWrite.h"
#include "GlobalVar.h"
#include "SomeFunctions.h"
#include <mri_core_data.h>

namespace Gadgetron
{
	class EXPORTCSLAB CS_Retro_PostBARTGadget : public Gadget1<IsmrmrdImageArray>
	{
	public:
		CS_Retro_PostBARTGadget();
		~CS_Retro_PostBARTGadget();

		GADGET_DECLARE(CS_Retro_PostBARTGadget);

	protected:
		int process_config(ACE_Message_Block *mb);
		int process(GadgetContainerMessage<IsmrmrdImageArray> *m1);
	};
}

#endif //CS_RETRO_POSTBARTGADGET_H
