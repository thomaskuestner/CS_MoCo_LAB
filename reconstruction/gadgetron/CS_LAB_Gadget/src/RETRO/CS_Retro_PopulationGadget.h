/*
file name	: 	CS_Retro_PopulationGadget.h
author		: 	Martin Schwartz (martin.schwartz@med.uni-tuebingen.de)
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
		std::vector<float> navigator_card_interpolated_, navigator_resp_interpolated_;

		// centroids of gates
		std::vector<float> respiratory_centroids_, cardiac_centroids_;

	public:
		CS_Retro_PopulationGadget();
		~CS_Retro_PopulationGadget();

		GADGET_DECLARE(CS_Retro_PopulationGadget);

	protected:
		int process_config(ACE_Message_Block *mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m3);

	private:
		bool fDiscard();
		bool get_cardiac_gates(int cardiac_gate_count);
		bool get_respiratory_gates(int respiratory_gate_count);
		bool fPopulatekSpace(int iNoGates);
		void calculate_weights(std::vector<float> &weights, const int population_mode, const int phase);

		template <typename T>
		void discard_empty_elements_from_back(std::vector<T> &v) {
			while (v.size() > 0 && v.at(v.size()-1) == 0) {
				v.pop_back();
			}
		}

		/**
		* @brief returns populated data as hoNDArray via first argument. Dimensions: [RX Channels] (e.g. [256 10])
		*/
		void get_populated_data(hoNDArray<std::complex<float> > &populated_data, const int population_mode, const hoNDArray<std::complex<float> > &unordered, const std::vector<size_t> &indices, const std::vector<float> &centroid_distances);

	public:
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		// declare gadget properties
		GADGET_PROPERTY(PopulationMode, int, "PopulationMode", 0);
		GADGET_PROPERTY(GatingMode, int, "GatingMode", 0);
		GADGET_PROPERTY(Tolerance, float, "Tolerance", 1.0);
#endif
	};
} // close namespace Gadgetron
#endif // CS_RETRO_POPULATIONGADGET_H
