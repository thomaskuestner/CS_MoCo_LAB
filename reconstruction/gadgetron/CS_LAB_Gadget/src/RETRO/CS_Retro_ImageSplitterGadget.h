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

	private:
		unsigned int simultaneous_cardiac_phases_, simultaneous_respiratory_phases_;

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	public:
		GADGET_PROPERTY(SimultaneousCardiacPhases, int, "SimultaneousCardiacPhases", 0);
		GADGET_PROPERTY(SimultaneousRespiratoryPhases, int, "SimultaneousRespiratoryPhases", 0);
#endif
	};
}

#endif //CS_RETRO_IMAGESPLITTERGADGET_H
