/*	
file name	: 	GlobalVar.h
author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
version		: 	1.0
date		: 	03.01.2015
description	: 	singleton class implementation for global variables				
variables	:	bESPRActiveCS_		:	ESPReSSo constraint is active for purely CS data set without Partial Fourier sampling (true: active)
				cfLambdaESPReSSo_	:	Lagrangian multiplier for the ESPReSSo constraint
				cfLambda_			:	Lagrangian multiplier for the FOCUSS constraint
				fFullySampled_		:	fully sampled region (in percent)
				fPartialFourierVal_	:	Partial Fourier value (4/8, 5/8, 6/8, 7/8)
references	:	http://de.wikibooks.org/wiki/C%2B%2B-Programmierung:_Entwurfsmuster:_Singleton
*/

#ifndef GlobalVar_H
#define GlobalVar_H

#pragma once
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"
#include "hoNDArray.h"
#include "hoNDKLT_CS.h"

namespace Gadgetron{

	class GlobalVar{
		public:
			static GlobalVar* instance(){
				static CGuard g;
				if (!_instance)
					_instance = new GlobalVar();
				return _instance;
			}

		//std::vector< hoNDArray< std::complex< float > > * > vfPrincipleComponents_;
		std::vector< bool > vbStatPrinc_;
		std::vector< hoNDKLT_CS< std::complex< float > > * > KLTVec_;

		// acquisition header
		std::vector< GadgetContainerMessage<ISMRMRD::AcquisitionHeader> * > AcqVec_;

		// image header
		std::vector< GadgetContainerMessage<ISMRMRD::ImageHeader> * > ImgHeadVec_;

		// nav indices
		std::vector< float > vNavInd_;

		// loop counter vectors
		std::vector< int > vPE_;
		std::vector< int > vPA_;

		// TR of sequence
		float fTR_;

		// FOCUSS parameters
		int iNOuter_;				//k-t FOCUSS loops
		int iNInner_;				// CG loops
		int iVDMap_;				// density map
		int iESPReSSoDirection_;	// ESPReSSo direction (y: 1, z: 2)
		int iDimFFT_;
		int iDimDCTSparse_;
		int iDimPCASparse_;
		int iDimKernelFFT_;
		int iScrambleDim_;
		int iTransformFFTBA_;
		int ikSpaceOut_;
		
		float fPartialFourierVal_;
		float fCSAcc_;
		float fFullySampled_;		// sequence parameters

		std::complex< float > cfLambdaESPReSSo_;	// lambda for ESPReSSo conjugate similarity
		std::complex< float > cfLambda_;			//FOCUSS stability in noisy environment (default:5, max:75)

		bool bESPRActiveCS_;		// using ESPReSSo constraint for non-ESPReSSo acquisitions

		// Retro Vars
		int iNoGates_;
		int iGatingMode_;			// Gating mode (0: percentile, 1: kMeans)
		int iMeasurementTime_;
		int iNavPERes_;
		int iNavPeriod_;
		int iNPhases_;
		int iPopulationMode_;		// mode for k-space population (0: closes, 1: average, 2: collect)

		private:
			static GlobalVar* _instance;
			GlobalVar() {}
			GlobalVar(const GlobalVar&);
			~GlobalVar(){}

			class CGuard{
				public:
					~CGuard()
					{
						if( NULL != GlobalVar::_instance)
						{
							delete GlobalVar::_instance;
							GlobalVar::_instance = NULL;
						}
					}
			};
			friend class CGuard;
	};	
}

#endif  //GlobalVar_H