#ifndef CS_RETRO_NAVIGATORGADGET_H
#define CS_RETRO_NAVIGATORGADGET_H

#pragma once
#include "Gadget.h"
#include "hoNDArray.h"
#include "CS_LAB_export.h"
#include <ismrmrd.h>
#include <complex>
#include "GlobalVar.h"
#include "hoNDFFT_CS.h"
#include "SomeFunctions.h"
#include "hoNDArray_blas.h"
#include "hoNDArray_math_util.h"
#include "hoNDImage_util.h"
#include "hoNDBSpline.h"
#include "hoNDKLT.h"
#include "hoNDFFT.h"
#include <cmath>

#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron {
	class EXPORTCSLAB CS_Retro_NavigatorGadget : public Gadget3<ISMRMRD::ImageHeader, hoNDArray<std::complex<float> >, hoNDArray<std::complex<float> > >
	{
	private:
		// navigator signal interpolated to TRs
		std::vector<float> navigator_resp_interpolated_;

		// navigator signal
		std::vector<float> vNav_;

		// field of view
		float field_of_view_[3];

		// number of acquisition
		long lNoScans_;

		// number of channels
		int iNoChannels_;

		// number of navigator acquisitions
		int iNoNav_;

		// navigation method (set by gadget property)
		int iNavMethod_;

		// min/max frequencies of resp/card
		float min_card_freq_, max_card_freq_, min_resp_freq_, max_resp_freq_;

		// Principal Component search range (index count MATLAB-like: counting from 1)
		const size_t search_range_min_ = 1;
		const size_t search_range_max_ = 15;

	public:
		CS_Retro_NavigatorGadget();
		~CS_Retro_NavigatorGadget();

		GADGET_DECLARE(CS_Retro_NavigatorGadget);

	protected:
		int process_config(ACE_Message_Block* mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m3);

	private:
		void getNav2D(hoNDArray<std::complex<float> > &aNav);
		void getNav2DPCA(hoNDArray<std::complex<float> > &aNav);
		void butterworth_filtering(const double fl, const double fh, std::vector<float> &signal);

	public:
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		GADGET_PROPERTY(NavigationMethod, int, "NavigationMethod", 0);
		GADGET_PROPERTY(MinRespFreq, float, "MinRespFreq", 7.5);
		GADGET_PROPERTY(MaxRespFreq, float, "MaxRespFreq", 40.0);
		GADGET_PROPERTY(MinCardFreq, float, "MinCardFreq", 40.0);
		GADGET_PROPERTY(MaxCardFreq, float, "MaxCardFreq", 150.0);
#endif
	};
} // close namespace Gadgetron
#endif // CS_RETRO_NAVIGATORGADGET_H
