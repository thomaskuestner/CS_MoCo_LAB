#ifndef CS_RETRO_POPULATIONGADGET_H
#define CS_RETRO_POPULATIONGADGET_H

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

#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron{
	
	class EXPORTCSLAB CS_Retro_PopulationGadget : public Gadget3< ISMRMRD::ImageHeader, hoNDArray<float>, hoNDArray<std::complex<float> > >
    {
    public:      
		CS_Retro_PopulationGadget();
		~CS_Retro_PopulationGadget();
		int process_config(ACE_Message_Block* mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray<float> >* m2, GadgetContainerMessage< hoNDArray<std::complex<float> > >* m3);
		GADGET_DECLARE(CS_Retro_PopulationGadget);
	  
    //private:
		bool fDiscard();
		bool fCalcCentroids(int iNoGates);
		bool fPopulatekSpace(int iNoGates);

		// populated/reordered kspace
		hoNDArray<std::complex<float>> hacfKSpace_reordered_;

		// unordered k-space
		hoNDArray<std::complex<float>> hacfKSpace_unordered_;

		// array dimensions of unordered k-space
		std::vector<size_t> vtDims_unordered_;

		// dimensions of incoming array
		std::vector<size_t> dimensionsIn_;
		
		// gating mode
		int iGatingMode_;

		// population mode
		int iPopulationMode_;

		// number of channels
		int iNoChannels_;

		// echo in k-space line
		int iEchoLine_;

		// echo partition in k-space
		int iEchoPartition_;

		// vector containing the tolerance values
		std::vector<float> vTolerance_;

		// tolerance/blending factor
		float fTolerance_;

		// navigator signal
		std::vector<float> vNavInt_;

		// centroids of gates
		std::vector<float> vfCentroids_;
		
		// TR of sequence
		float fTR_;

		bool bMatlab_;
		std::vector<int> vPE_;
		std::vector<int> vPA_;
    };
}
#endif //CS_RETRO_POPULATIONGADGET_H
