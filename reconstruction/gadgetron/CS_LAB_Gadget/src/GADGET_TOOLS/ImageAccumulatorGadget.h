#ifndef IMAGEACCUMULATORGADGET_H
#define IMAGEACCUMULATORGADGET_H

#pragma once
#include "CS_LAB_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include <ismrmrd.h>
#include <complex>

#include "GlobalVar.h"
#include "SomeFunctions.h"

#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron{
	
	class EXPORTCSLAB ImageAccumulatorGadget : public Gadget2< ISMRMRD::ImageHeader, hoNDArray< float > >
    {
    public:      
      ImageAccumulatorGadget();
      ~ImageAccumulatorGadget();
      int process_config(ACE_Message_Block* mb);
	  int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage<hoNDArray< float >>* m2);
	  GADGET_DECLARE(ImageAccumulatorGadget);
	  
	  int fCopyHeader(GadgetContainerMessage<ISMRMRD::ImageHeader>* tmp_m1, GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);

    protected:
		std::vector<size_t> vtDimensions_;
		int iPartition_;
		int iPhs_;
		int iImageLoopCounter_;
		hoNDArray< float >* hafBuffer_;
    };
}
#endif //IMAGEACCUMULATORGADGET_H
