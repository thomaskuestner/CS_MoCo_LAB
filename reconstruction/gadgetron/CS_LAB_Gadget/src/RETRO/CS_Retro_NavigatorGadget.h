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
#include <cmath>
#include <math.h>

#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron{
	
	class EXPORTCSLAB CS_Retro_NavigatorGadget : public Gadget3< ISMRMRD::ImageHeader, hoNDArray< std::complex<float> >, hoNDArray<std::complex<float>> >
    {
    public:      
		CS_Retro_NavigatorGadget();
		~CS_Retro_NavigatorGadget();
		int process_config(ACE_Message_Block* mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader>*m1, GadgetContainerMessage<hoNDArray<std::complex<float>>>*m2, GadgetContainerMessage<hoNDArray<std::complex<float>>>* m3);
		GADGET_DECLARE(CS_Retro_NavigatorGadget);
	  
    
		bool CS_Retro_NavigatorGadget::getNav2D(hoNDArray<std::complex<float>> &aNav);

		// navigator index vector
		std::vector<float> vNavInd_;

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

		// navigator period
		int iNavPeriod_;

		// matlab debug output
		bool bMatlab_;

    };
}
#endif //CS_RETRO_NAVIGATORGADGET_H
