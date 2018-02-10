#ifndef CS_ACCUMULATORGADGET_H
#define CS_ACCUMULATORGADGET_H

#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
	//#define GET_MACRO(_1,_2,_3,NAME,...) NAME
	//#define GDEBUG(...) GET_MACRO(__VA_ARGS__, GADGET_DEBUG1, GADGET_DEBUG2)(__VA_ARGS__)
	#define GADGET_DEBUG1(...) GDEBUG(__VA_ARGS__)
	#define GADGET_DEBUG2(x, ...) GDEBUG(x, ##__VA_ARGS__)
	#define GADGET_DEBUG_EXCEPTION(x,y) GEXCEPTION(x,y)
#endif

#pragma once
#if __GADGETRON_VERSION_HIGHER_3_6__ == 0
#include "ismrmrd.hxx"
#endif
#include "Gadget.h"
#include "hoNDArray.h"
#include "CS_LAB_export.h"
#include <ismrmrd.h>
#include <complex>
#include "SomeFunctions.h"
#include "GlobalVar.h"
#include "GadgetIsmrmrdReadWrite.h"

#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
	#include "xml.h"
#else	
	#include "ismrmrd/xml.h"
#endif

namespace Gadgetron{
	
	class EXPORTCSLAB CS_AccumulatorGadget : public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:      
      CS_AccumulatorGadget();
      ~CS_AccumulatorGadget();
      int process_config(ACE_Message_Block* mb);
	  int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
	  GADGET_DECLARE(CS_AccumulatorGadget);
	  
    protected:
	  int fCopyData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* GC_acq_m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* GC_img_m2, std::complex<float>* pcfBuffer);
	  //int fCopyHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* GC_acq_m1, GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* GC_acq_m1_new); 
      hoNDArray< std::complex<float> >* hacfBuffer_;
	  std::vector<size_t> vFOV_;
      std::vector<size_t> vDim_;
      long long image_counter_;
      long long image_series_;
	  int iBodyRegion_, iSamplingType_;
	  float fFullySa_, fLESPReSSo_, fLQ_;
    };
}
#endif //CS_ACCUMULATORGADGET_H
