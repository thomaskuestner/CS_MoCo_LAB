/*	
file name	: 	GlobalVar_FOCUSS.h

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	singleton class implementation for global variables
				
variables	:	vPrincipleComponents_:	principle components for the PCA transformation
				iStatPrinc_	:	boolean vector - true, if principle components for this dimensions already exist
				
references	:	http://de.wikibooks.org/wiki/C%2B%2B-Programmierung:_Entwurfsmuster:_Singleton
*/

#ifndef GlobalVar_FOCUSS_H
#define GlobalVar_FOCUSS_H

#pragma once
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"
#include "hoNDArray.h"
#include "CS_LAB_export.h"

namespace Gadgetron{

	class EXPORTCSLAB GlobalVar_FOCUSS{
		public:
			static GlobalVar_FOCUSS* instance(){
				static CGuard g;
				if (!_instance)
					_instance = new GlobalVar_FOCUSS();
				return _instance;
			}
			
		// wisdom string vector
		std::vector<hoNDArray<std::complex<float> > *> vfPrincipleComponents_;
		std::vector<bool> vbStatPrinc_;

		// FOCUSS parameters
		int iNOuter_;
		int iNInner_;
		int iCGResidual_;
		bool bESPRActiveCS_;
		int iVDMap_;
		float fFullySampled_;
		std::complex< float > cfLambdaESPReSSo_;
		std::complex< float > cfLambda_;
		int iESPReSSoDirection_;
		float fPartialFourierVal_;
		float fCSAcc_;
		int iDimFFT_;
		int iDimDCTSparse_;
		int iDimPCASparse_;
		int iDimKernelFFT_;
		int iTransformFFTBA_;
		int ikSpaceOut_;

		private:
			static GlobalVar_FOCUSS* _instance;
			GlobalVar_FOCUSS() {}
			GlobalVar_FOCUSS(const GlobalVar_FOCUSS&);
			~GlobalVar_FOCUSS(){}

			class CGuard{
				public:
					~CGuard()
					{
						if( NULL != GlobalVar_FOCUSS::_instance)
						{
							delete GlobalVar_FOCUSS::_instance;
							GlobalVar_FOCUSS::_instance = NULL;
						}
					}
			};
			friend class CGuard;
	};	
}

#endif  //GlobalVar_FOCUSS_H
