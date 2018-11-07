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

		// center mask
		hoNDArray<bool> mask_center_;

		// samples to fill center
		hoNDArray<std::complex<float> > center_samples_;

		// parameters for center mask
		float low_res_vs_;			//% low-resolution view-sharing [0,1] (0: none, 1: use complete fully sampled center)
		bool omit_center_vs_;		//% omit DC component (1D navigator) to be shared amongst all motion states

		// number of channels
		int iNoChannels_;

		// echo in k-space line
		int iEchoLine_;

		// echo partition in k-space
		int iEchoPartition_;

		// vector containing the tolerance values
		std::vector<float> respiratory_tolerance_vector_;

		// tolerance/blending factor
		float cardiac_tolerance_parameter_, respiratory_tolerance_parameter_;

		// navigator signal
		std::vector<float> navigator_card_interpolated_, navigator_resp_interpolated_;

		// centroids of gates
		std::vector<float> respiratory_centroids_;
		std::vector<unsigned int> cardiac_gates_;

	public:
		CS_Retro_PopulationGadget();
		~CS_Retro_PopulationGadget();

		GADGET_DECLARE(CS_Retro_PopulationGadget);

	protected:
		int process_config(ACE_Message_Block *mb);
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m3);

	private:
		bool fDiscard();
		bool get_cardiac_gates(const unsigned int cardiac_gate_count, const float f_s);
		bool get_respiratory_gates(const unsigned int respiratory_gate_count);
		bool fPopulatekSpace(const unsigned int cardiac_gate_count, const unsigned int respiratory_gate_count);
		std::vector<float> calculate_weights(const unsigned int phase);

		template <typename T>
		void discard_empty_elements_from_back(std::vector<T> &v) {
			while (v.size() > 0 && v.at(v.size()-1) == 0) {
				v.pop_back();
			}
		}

		/**
		* @brief returns populated data as hoNDArray via first argument. Dimensions: [RX Channels] (e.g. [256 10])
		* @param indices The indices of the interesting data
		* @param centroid_distances Distances of respiratory phase centroid
		* @param cardiac_gate_count Number of cardiac gates
		* @param line Ky in which the data should be copied (important for view sharing)
		* @param partition Kz in which the data should be copied (important for view sharing)
		* @param respiratory_phase The respiratory phase the data is populated for
		*/
		hoNDArray<std::complex<float> > get_populated_data(const std::vector<size_t> &indices, const std::vector<float> &centroid_distances, const unsigned int cardiac_gate_count, const unsigned int line, const unsigned int partition, const unsigned int respiratory_phase);

		template <typename T> void remove_high_peaks(std::vector<T> &signal);
		template <typename T> void search_peaks(const std::vector<T> &signal, std::vector<size_t> &x_pos);

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	public:
		// declare gadget properties
		GADGET_PROPERTY(PopulationMode, int, "PopulationMode", 0);
		GADGET_PROPERTY(CardiacGatingMode, int, "CardiacGatingMode", 0);
		GADGET_PROPERTY(RespiratoryGatingMode, int, "RespiratoryGatingMode", 0);
		GADGET_PROPERTY(CardiacTolerance, float, "CardiacTolerance", 1.0);
		GADGET_PROPERTY(RespiratoryTolerance, float, "RespiratoryTolerance", 1.0);
		GADGET_PROPERTY(LowResVS, float, "LowResVS", 0.0);
		GADGET_PROPERTY(OmitCenterVS, int, "OmitCenterVS", 0);
#endif
	};
} // close namespace Gadgetron
#endif // CS_RETRO_POPULATIONGADGET_H
