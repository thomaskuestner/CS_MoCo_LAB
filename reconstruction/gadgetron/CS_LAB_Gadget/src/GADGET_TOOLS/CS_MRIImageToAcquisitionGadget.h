#ifndef CS_MRIIMAGETOACQUISITIONGADGET_H
#define CS_MRIIMAGETOACQUISITIONGADGET_H

#pragma once
#include "Gadget.h"
#include "NDArray.h"
#include "GadgetMRIHeaders.h"
#include "gadgetron_mricore_export.h"
#include <ismrmrd.h>
#include <complex>
#include "CS_GlobalVar.h"
#include "GadgetMessageInterface.h"
#include "GadgetStreamController.h"
#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron{
	
	class CS_MRIImageToAcquisitionGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:
		GADGET_DECLARE(CS_MRIImageToAcquisitionGadget);
      
    protected:
		virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

		// correct header information (used for Compressed Sensing data - incomplete information)
		int fCorrectHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_hdr_m1, int iLine, int iPartition, int iPhase);
		int fCopyHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1, GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1_new); 
      
		// array dimension
		std::vector<size_t> vDims_;
	};
}

#endif //CS_MRIIMAGETOACQUISITIONGADGET_H
