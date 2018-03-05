/*	
file name	: 	CS_LAB.h
author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
version		: 	1.0
date		: 	16.12.2015
description	: 	evaluates the object dimensions of the incoming data set and creates the respective class object.
input		:	m1					: 	appended acquisition header information
				m2					: 	data set				
output		:	m1					: 	untouched header information
				m2					:	reconstructed image (FOCUSS algorithm)				
functions	:	process(...)		:	data processing - create FOCUSS object and reconstruct k-space data
				GADGET_DECLARE(...)	:	Gadget declaration (for the Gadgetron pipeline)
				fGradESPReSSo)...)	: 	empty function definition (CS_FOCUSS is abstract class)
				fInitESPReSSo(...)	: 	empty function definition (CS_FOCUSS is abstract class)
				fWindowing(...)		:	empty function definition (CS_FOCUSS is abstract class)				
variables	:	opCS_				:	pointer to CS_FOCUSS object
				iAlgorithm_			:	algorithm selection
				iDataset_			:	dataset type				
references	:	-
notes		:	This class have to be called in the Gadgetron XML configuration file for automatically FOCUSS algorithm selection
*/

#ifndef CS_LAB_H
#define CS_LAB_H

#pragma once
#include "CS_LAB_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include <complex>
#include <ctime>
#include <ismrmrd.h>
#include "CS_FOCUSS.h"

#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron
{
	class EXPORTCSLAB CS_LAB : public CS_FOCUSS
	{
		GADGET_DECLARE(CS_LAB)
	public:
		int fRecon(hoNDArray<std::complex<float> >  &hacfInput, hoNDArray<std::complex<float> >  &hacfRecon){ return true; };

		void fExternalControl();

		int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

		void fGradESPReSSo(hoNDArray<std::complex<float> > & hacfRho, hoNDArray<std::complex<float> > &hacfFullMask, hoNDArray<std::complex<float> > &hacfKSpace, hoNDArray<std::complex<float> > &hacfW, hoNDArray<std::complex<float> > &hacfQ){};

		void fInitESPReSSo(hoNDArray<bool>& habFullMask){};

		void fWindowing(hoNDArray<std::complex<float> > & hacfWWindowed){};

		// pointer to CS_FOCUSS class object
		CS_FOCUSS *opCS_;

		// algorithm selection
		int iAlgorithm_;

		// dataset?
		int iDataset_;

		// Matlab transformation parameters..
		int iFFT_Sparse_, iDCT_Sparse_, iPCA_Sparse_, iKernel_FFT_dim_, iScrambleDim_, iFFTBA_, kSpaceOut_, iTime_;
};

}
#endif //CS_LAB_H
