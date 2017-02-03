/*
file name	: 	CS_Control.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	implementation of the class "CS_Control.h"

references	:	-
*/

#include "CS_Control.h"

#include "GadgetIsmrmrdReadWrite.h"

using namespace Gadgetron;
/*int CS_CONTROL::process_config(ACE_Message_Block* mb){

	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		bXMLControl_ = bXMLControl.value();
	#else
		bXMLControl_ = this->get_bool_value("XMLControl");
	#endif	

	if (bXMLControl_) {

		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1		
			iNOuter_ = OuterIterations.value();
		  	iCGResidual_ = residual.value();
			iNInner_ = InnerIterations.value();
			bESPRActiveCS_ = CSESPReSSo.value();
			cfLambda_ = lambda.value();
		#else
		  	// how to calculate the beta value
			iCGResidual_ = this->get_int_value("residual");	

			// maximum number of FOCUSS iterations	
		  	iNOuter_ = this->get_int_value("OuterIterations");
		  	
			// maximum number of CG iterations	
		  	iNInner_ = this->get_int_value("InnerIterations");
		  	
			// use ESPReSSo-constraint for pure CS data
			bESPRActiveCS_ = this->get_bool_value("CSESPReSSo");

			// FOCUSS lambda
			cfLambda_ = this->get_double_value("lambda");
		#endif
	}
	else{
		iNOuter_ = GlobalVar::instance()->iNOuter_;
		iNInner_ = GlobalVar::instance()->iNInner_;
		iCGResiudal_ = GlobalVar::instance()->iCGResiudal_;
		bESPRActiveCS_ = GlobalVar::instance()->bESPRActiveCS_;
		cfLambda_ = GlobalVar::instance()->lambda_;
	
	}

	if (iNInner_ <= 0) iNInner_ = 20;
	if (iNOuter_ <= 0) iNOuter_ = 2;	

	// p-value for the lp-norm
	fP_ = .5;

	// convergence boundary
	fEpsilon_ = (float)1e-6;

	// setup of the transformation parameters - sparsity dim, fft dim, ..
	fSetupTransformation();

	return GADGET_OK;
};*/

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
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("not implemented in this version\n");
		#else
			GADGET_DEBUG1("not implemented in this version\n");
		#endif
	}

	// set parameters of the FOCUSS class - required, because the xml config file is read in by CS_CONTROL class and not by FOCUSS class
	pCS->iCGResidual_				= iCGResidual_;
	pCS->iNChannels_				= iNChannels_;
	pCS->iNOuter_					= iNOuter_;
	pCS->iNInner_					= iNInner_;
	pCS->fP_						= fP_;
	pCS->cfLambda_					= cfLambda_;
	pCS->cfLambdaESPReSSo_			= cfLambdaESPReSSo_;
	pCS->fEpsilon_					= fEpsilon_;
	pCS->fCSAccel_					= fCSAccel_;
	pCS->iESPReSSoDirection_		= iESPReSSoDirection_;
	pCS->fPartialFourierVal_		= fPartialFourierVal_;
	pCS->bESPRActiveCS_				= bESPRActiveCS_;
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
