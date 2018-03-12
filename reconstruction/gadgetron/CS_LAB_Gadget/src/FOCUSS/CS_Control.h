/*	
file name	: 	CS_Control.h
author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
version		: 	1.0
date		: 	03.01.2015
description	: 	evaluates the object dimensions of the incoming data set and creates the respective class object. This class is derived from "CS_FOCUSS" and is based on the Composite design pattern. For a detailed description it is suggested to read p. 55 of the thesis.
input		:	m1					: 	appended acquisition header information
				m2					: 	data set				
output		:	m1					: 	untouched header information
				m2					:	reconstructed image (FOCUSS algorithm)
functions	:	process(...)		:	data processing - create FOCUSS object and reconstruct k-space data
				GADGET_DECLARE(...)	:	Gadget declaration (for the Gadgetron pipeline)
				gradTV(...)			:	empty function definition (CS_FOCUSS is abstract class)
				gradESPReSSo)...)	: 	empty function definition (CS_FOCUSS is abstract class)
				initESPReSSo(...)	: 	empty function definition (CS_FOCUSS is abstract class)
				windowing(...)		:	empty function definition (CS_FOCUSS is abstract class)				
variables	:	pCS					:	pointer to CS_FOCUSS object				
references	:	-
notes		:	This class have to be called in the Gadgetron XML configuration file for automatically FOCUSS algorithm selection
*/

#ifndef CS_CONTROL_H
#define CS_CONTROL_H

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
	class EXPORTCSLAB CS_CONTROL : public CS_FOCUSS
	{
	public:
		GADGET_DECLARE(CS_CONTROL)

	protected:
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray< std::complex<float> > > *m2);

		void fGradESPReSSo(hoNDArray<std::complex<float> > &hacfRho, hoNDArray<std::complex<float> > &hacfFullMask, hoNDArray<std::complex<float> > &hacfKSpace, hoNDArray<std::complex<float> > &hacfW, hoNDArray<std::complex<float> > &hacfQ) {}

		void fInitESPReSSo(hoNDArray<bool> &habFullMask) {}

		void fWindowing(hoNDArray<std::complex<float> > &hacfWWindowed) {}

		int fRecon(hoNDArray<std::complex<float> > &hacfInput, hoNDArray<std::complex<float> > &hacfRecon) { return GADGET_OK; };

		// pointer to CS_FOCUSS class object
		CS_FOCUSS *pCS;
	};
}

#endif //CS_CONTROL_H
