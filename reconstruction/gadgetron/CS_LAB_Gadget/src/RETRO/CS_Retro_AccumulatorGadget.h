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

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	#include "xml.h"
#else
	#include "ismrmrd/xml.h"
#endif

namespace Gadgetron {
	class EXPORTCSLAB CS_Retro_AccumulatorGadget : public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float> > >
	{
	private:
		std::vector<GadgetContainerMessage<hoNDArray<std::complex<float> > >*> buffer_kspace_;
		std::vector<hoNDArray<std::complex<float> >*> buffer_nav_;

		GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *acq_header_ = NULL;

		size_t dimensionsIn_[3];
		float field_of_view_[3];

		unsigned int iNoNav_ = 0;
		unsigned int iNoNavLine_ = 0;
		int iEchoLine_ = 0;
		int iEchoPartition_ = 0;

		// number of phases
		int iNPhases_ = 0;

		// Compressed Sensing variables
		int iESPReSSoDirection_ = 0;
		float fPartialFourierVal_ = 0;			// Partial Fourier value (4/8, 5/8, 6/8, 7/8)
		int iBodyRegion_ = 0;
		int iVDMap_ = 0;
		int iSamplingType_ = 0;
		float fCSAcc_ = 0;
		float fFullySa_ = 0.065;

	public:
		CS_Retro_AccumulatorGadget();
		~CS_Retro_AccumulatorGadget();

		GADGET_DECLARE(CS_Retro_AccumulatorGadget);

		int close(unsigned long flags) override;

	protected:
		int process_config(ACE_Message_Block *mb);
		int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2);

	private:
		bool is_image_dataset(ISMRMRD::AcquisitionHeader &header)
		{
			return header.idx.set == 0;
		}

		bool is_navigator_dataset(ISMRMRD::AcquisitionHeader &header)
		{
			return header.idx.set == 1;
		}

		bool process_data(void);

	public:
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		GADGET_PROPERTY(NavPeriod, int, "NavPeriod", 0);
		GADGET_PROPERTY(NavPERes, int, "NavPERes", 0);
		GADGET_PROPERTY(Phases, int, "Phases", 1);
#endif
	};
} // close namespace Gadgetron

#endif // CS_RETRO_ACCUMULATORGADGET_H
