/*	
file name	: 	GlobalVar.h
author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
version		: 	1.0
date		: 	03.01.2015
description	: 	singleton class implementation for global variables				
variables	:		
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

		// FOCUSS parameters
		int iNOuter_;
		int iNInner_;		
		int iVDMap_;		
		int iESPReSSoDirection_;		
		int iDimFFT_;
		int iDimDCTSparse_;
		int iDimPCASparse_;
		int iDimKernelFFT_;
		int iScrambleDim_;
		int iTransformFFTBA_;
		int ikSpaceOut_;
		
		float fPartialFourierVal_;
		float fCSAcc_;
		float fFullySampled_;

		std::complex< float > cfLambdaESPReSSo_;
		std::complex< float > cfLambda_;

		bool bESPRActiveCS_;

		// Retro Vars
		int iNoGates_;
		int iGatingMode_;
		int iMeasurementTime_;
		int iNavPERes_;
		int iNavPeriod_;
		int iNPhases_;
		int iPopulationMode_;

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