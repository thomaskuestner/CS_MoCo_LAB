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
	class EXPORTCSLAB CS_Retro_NavigatorGadget : public Gadget3< ISMRMRD::ImageHeader, hoNDArray< std::complex<float> >, hoNDArray<std::complex<float>> >
	{
	private:
		// navigator signal interpolated to TRs
		std::vector<float> vNavInt_;

		// navigator signal
		std::vector<float> vNav_;

		// field of view
		std::vector<float> field_of_view_;

		// number of acquisition
		long lNoScans_;

		// number of channels
		int iNoChannels_;

		// number of navigator acquisitions
		int iNoNav_;

		// matlab debug output
		bool bMatlab_;
		
		// navigation method (set by gadget property)
		int iNavMethod_;

	public:
		CS_Retro_NavigatorGadget();
		~CS_Retro_NavigatorGadget();

		GADGET_DECLARE(CS_Retro_NavigatorGadget);

	protected:
		int process_config(ACE_Message_Block* mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader>*m1, GadgetContainerMessage<hoNDArray<std::complex<float>>>*m2, GadgetContainerMessage<hoNDArray<std::complex<float>>>* m3);

	private:
		void getNav2D(hoNDArray<std::complex<float>> &aNav);
		void getNav2DPCA(hoNDArray<std::complex<float>> &aNav);

	public:
#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		GADGET_PROPERTY(NavigationMethod, int, "NavigationMethod", 0);
#endif
	};
} // close namespace Gadgetron
#endif //CS_RETRO_NAVIGATORGADGET_H
