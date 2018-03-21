#ifndef CS_ACCUMULATORGADGET_H
#define CS_ACCUMULATORGADGET_H

#include "gadgetron_messages.h"

#pragma once
#ifndef __GADGETRON_VERSION_HIGHER_3_6__
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

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	#include "xml.h"
#else
	#include "ismrmrd/xml.h"
#endif

namespace Gadgetron
{
	class EXPORTCSLAB CS_AccumulatorGadget : public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float> > >
	{
	public:
		CS_AccumulatorGadget();
		~CS_AccumulatorGadget();
		int process_config(ACE_Message_Block* mb);
		int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader > *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2);
		GADGET_DECLARE(CS_AccumulatorGadget);

	protected:
		int fCopyData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *GC_img_m2, std::complex<float> *pcfBuffer);
		hoNDArray<std::complex<float> > *hacfBuffer_;
		std::vector<size_t> vFOV_;
		std::vector<size_t> vDim_;
		long long image_counter_;
		long long image_series_;
		int iBodyRegion_, iSamplingType_;
		float fFullySa_, fLESPReSSo_, fLQ_;
		int iESPReSSoDirection_ = 0;	// ESPReSSo direction (y: 1, z: 2)
		int iVDMap_ = 0;				// density map
		int iNPhases_ = 0;
		float fCSAcc_ = 0;
		float fPartialFourierVal_ = 0;	// Partial Fourier value (4/8, 5/8, 6/8, 7/8)
	};
}

#endif // CS_ACCUMULATORGADGET_H
