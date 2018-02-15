#ifndef CS_RETRO_ACCUMULATORGADGET_H
#define CS_RETRO_ACCUMULATORGADGET_H

#include "gadgetron_messages.h"

#pragma once
#include "Gadget.h"
#include "hoNDArray.h"
#include "CS_LAB_export.h"
#include <ismrmrd.h>
#include <complex>
#include "GlobalVar.h"
#include "SomeFunctions.h"
#include "hoNDArray_blas.h"
#include "hoNDArray_math_util.h"
#include "hoNDImage_util.h"
#include "hoNDBSpline.h"
#include <cmath>

#include "GadgetIsmrmrdReadWrite.h"

#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
	#include "xml.h"
#else
	#include "ismrmrd/xml.h"
#endif

namespace Gadgetron {
	class EXPORTCSLAB CS_Retro_AccumulatorGadget : public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
	{
	private:
		hoNDArray< std::complex<float> >* bufferkSpace_;
		hoNDArray< std::complex<float> >* bufferNav_;

		std::vector<size_t> dimkSpace_;
		std::vector<size_t> dimNav_;

		long lCurrentScan_;
		std::vector<size_t> dimensionsIn_;
		std::vector<float> field_of_view_;
		long long image_counter_;
		long long image_series_;

		unsigned long lNoScans_;
		unsigned int iNoChannels_;
		unsigned int iNoNav_;
		unsigned int iNoNavLine_;
		int iEchoLine_;
		int iEchoPartition_;

		// CS_Retro variables
		int iBaseRes_;

		// number of phases
		int iNPhases_;

		// Compressed Sensing variables
		int iESPReSSoDirection_;
		float fPartialFourierVal_;
		int iBodyRegion_;
		int iVDMap_;
		int iSamplingType_;
		float fCSAcc_;
		float fFullySa_;

	public:
		CS_Retro_AccumulatorGadget();
		~CS_Retro_AccumulatorGadget();

		GADGET_DECLARE(CS_Retro_AccumulatorGadget);

	protected:
		int process_config(ACE_Message_Block* mb);
		int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>*m1, GadgetContainerMessage<hoNDArray<std::complex<float>>>*m2);

	public:
#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		GADGET_PROPERTY(NavPeriod, int, "NavPeriod", 0);
		GADGET_PROPERTY(NavPERes, int, "NavPERes", 0);
		GADGET_PROPERTY(MeasurementTime, int, "MeasurementTime", 0);
		GADGET_PROPERTY(Phases, int, "Phases", 0);
#endif
	};
} // close namespace Gadgetron
#endif //CS_RETRO_ACCUMULATORGADGET_H
