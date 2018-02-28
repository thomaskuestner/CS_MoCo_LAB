#ifndef IMAGESLICERGADGET_H
#define IMAGESLICERGADGET_H

#pragma once
#include "CS_LAB_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include <ismrmrd.h>
#include <complex>

#include "GlobalVar.h"
#include "SomeFunctions.h"

#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron
{
	class EXPORTCSLAB ImageSlicerGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray<float> >
	{
	public:
		ImageSlicerGadget();
		~ImageSlicerGadget();
		int process_config(ACE_Message_Block* mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2);
		GADGET_DECLARE(ImageSlicerGadget);
	};
}

#endif // IMAGESLICERGADGET_H
