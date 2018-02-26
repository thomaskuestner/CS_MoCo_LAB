/*
file name	: 	CS_Retro_PopulationGadget.h

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
				Thomas Kuestner (thomas.kuestner@med.uni-tuebingen.de)

version		: 	1.1

date		: 	13.10.2015: initial file
				15.01.2018: include Gaussian-weighted gating

description	: 	k-space population/gating
*/

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

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	#include "xml.h"
#else
	#include "ismrmrd/xml.h"
#endif

namespace Gadgetron {
	class EXPORTCSLAB CS_Retro_PopulationGadget : public Gadget3<ISMRMRD::ImageHeader, hoNDArray<float>, hoNDArray<std::complex<float> > >
	{
	private:
		// populated/reordered kspace
		hoNDArray<std::complex<float> > hacfKSpace_reordered_;

		// unordered k-space
		hoNDArray<std::complex<float> > hacfKSpace_unordered_;

		// array dimensions of unordered k-space
		std::vector<size_t> vtDims_unordered_;

		// dimensions of incoming array
		size_t dimensionsIn_[3];

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

		// number of gates
		int iNPhases_;

		// centroids of gates
		std::vector<float> vfCentroids_;

	public:
		CS_Retro_PopulationGadget();
		~CS_Retro_PopulationGadget();

		GADGET_DECLARE(CS_Retro_PopulationGadget);

	protected:
		int process_config(ACE_Message_Block *mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m3);

	private:
		bool fDiscard();
		bool fCalcCentroids(int iNoGates);
		bool fPopulatekSpace(int iNoGates);

	public:
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		// declare gadget properties
		GADGET_PROPERTY(Gates, int, "Gates", 4);
		GADGET_PROPERTY(PopulationMode, int, "PopulationMode", 0);
		GADGET_PROPERTY(GatingMode, int, "GatingMode", 0);
#endif
	};
} // close namespace Gadgetron
#endif //CS_RETRO_POPULATIONGADGET_H
