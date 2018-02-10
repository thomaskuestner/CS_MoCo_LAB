/*
file name	: 	CS_Control.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	implementation of the class "CS_Control.h"

references	:	-
*/

#include "CS_Control.h"

using namespace Gadgetron;

int CS_CONTROL::process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2){

	// get dimension of the incoming data object
	std::vector<size_t> vDims = *m2->getObjectPtr()->get_dimensions();

	// copy GadgetContainer and init with m2 data
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* tmp_m2 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
	tmp_m2->getObjectPtr()->create(vDims);
	memcpy(tmp_m2->getObjectPtr()->get_data_ptr(), m2->getObjectPtr()->get_data_ptr(), m2->getObjectPtr()->get_number_of_elements()*sizeof(std::complex< float >));
	
	// evaluate dimension and create suitable class object
	if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) == 1 && vDims.at(3) == 1){
		pCS = new CS_FOCUSS_2D();
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Incoming data is 2D - starting 2D FOCUSS reconstruction\n");
		#else
			GADGET_DEBUG1("Incoming data is 2D - starting 2D FOCUSS reconstruction\n");
		#endif
	}
	else if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) == 1 && vDims.at(3) > 1){
		//pCS = new CS_FOCUSS_2Dt();
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Incoming data is 2Dt - starting 2Dt FOCUSS reconstruction\n");
		#else
			GADGET_DEBUG1("Incoming data is 2Dt - starting 2Dt FOCUSS reconstruction\n");
		#endif
	}
	else if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) > 1 && vDims.at(3) == 1){
		// squeeze array due to x,y,z,c dimension of 3D FOCUSS class
		sum_dim(*tmp_m2->getObjectPtr(), 3, *tmp_m2->getObjectPtr());
		
		//tmp_m2->getObjectPtr()->print(std::cout);

		pCS = new CS_FOCUSS_3D();
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Incoming data is 3D - starting 3D FOCUSS reconstruction\n");
		#else
			GADGET_DEBUG1("Incoming data is 3D - starting 3D FOCUSS reconstruction\n");
		#endif
	}
	else if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) > 1 && vDims.at(3) > 1){
		
		pCS = new CS_FOCUSS_4D();

		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Incoming data is 4D - starting 4D FOCUSS reconstruction\n");
		#else
			GADGET_DEBUG1("Incoming data is 4D - starting 4D FOCUSS reconstruction\n");
		#endif
	}

	// set parameters of the FOCUSS class - required, because the xml config file is read in by CS_CONTROL class and not by FOCUSS class
	//pCS->iCGResidual_				= iCGResidual_;
	pCS->iNChannels_				= iNChannels_;
	pCS->fP_						= fP_;
	pCS->fEpsilon_					= fEpsilon_;
	pCS->fCSAccel_					= fCSAccel_;
	pCS->hacfFilter_				= hacfFilter_;
	pCS->Transform_KernelTransform_ = Transform_KernelTransform_;
	pCS->Transform_fftBA_			= Transform_fftBA_;
	pCS->Transform_fftAA_			= Transform_fftAA_;

	// disable standalone Gadget behaviour
	pCS->bControl_	= true;
	pCS->bDebug_	= true;
	pCS->bMatlab_	= false;

	// process data in class member function
	pCS->process(m1, tmp_m2);

	//Now pass on image
	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}
	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_CONTROL)
