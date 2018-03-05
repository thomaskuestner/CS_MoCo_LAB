#ifndef ACQUISITIONFINISHGADGET_H
#define ACQUISITIONFINISHGADGET_H

#include "CS_LAB_export.h"
#include "Gadget.h"
#include "NDArray.h"
#include "GadgetMRIHeaders.h"
#include "GadgetMessageInterface.h"
#include "GadgetStreamController.h"
#include <ismrmrd.h>
#include <complex>

#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron{

  class EXPORTCSLAB AcquisitionFinishGadget : public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(AcquisitionFinishGadget);
      
    protected:
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
    };
}

#endif //ACQUISITIONFINISHGADGET_H
