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
int CS_CONTROL::process_config(ACE_Message_Block* mb){

	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		GDEBUG("process config..\n");
	#else
		GADGET_DEBUG1("process config..\n");
	#endif	
	bXMLControl_ = true;
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		bXMLControl_ = bXMLControl.value();
		GDEBUG("XML Control enabled..\n");
	#else
		bXMLControl_ = this->get_int_value("bXMLControl");
		GADGET_DEBUG1("XML Control enabled..\n");
	#endif	

	if (bXMLControl_) {

		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1		
			iNOuter_ = OuterIterations.value();
		  	iCGResidual_ = iCGResidual.value();
			iNInner_ = InnerIterations.value();
			bESPRActiveCS_ = CSESPReSSo.value();
			cfLambda_ = lambda.value();			
			cfLambdaESPReSSo_ = lambdaESPReSSo.value();
			int iDimFFT = fftSparseDim.value();
			int iDimDCTSparse = dctSparseDim.value();
			int iDimPCASparse = pcaSparseDim.value();
			int iDimKernelFFT = kernelFftDim.value();
			int iTransformFFTBA = transformFftBaDim.value();
			int ikSpaceOut = kSpaceOutDim.value();
		#else
		  	// how to calculate the beta value
			iCGResidual_ = this->get_int_value("iCGResidual");
		  	iNOuter_ = this->get_int_value("OuterIterations");		  	
		  	iNInner_ = this->get_int_value("InnerIterations");		  
			bESPRActiveCS_ = this->get_int_value("CSESPReSSo");
			cfLambda_ = this->get_double_value("lambda");
			cfLambdaESPReSSo_ = this->get_double_value("lambdaESPReSSo");
			int iDimFFT = this->get_int_value("fftSparseDim");
			int iDimDCTSparse = this->get_int_value("dctSparseDim");
			int iDimPCASparse = this->get_int_value("pcaSparseDim");
			int iDimKernelFFT = this->get_int_value("kernelFftDim");
			int iTransformFFTBA = this->get_int_value("transformFftBaDim");
			int ikSpaceOut = this->get_int_value("kSpaceOutDim");
		#endif

		// update global parameters
		GlobalVar_FOCUSS::instance()->iNOuter_ = iNOuter_;
		GlobalVar_FOCUSS::instance()->iNInner_ = iNInner_;
		GlobalVar_FOCUSS::instance()->iCGResidual_ = iCGResidual_;
		GlobalVar_FOCUSS::instance()->bESPRActiveCS_ = bESPRActiveCS_;
		GlobalVar_FOCUSS::instance()->cfLambda_ = cfLambda_;	
		GlobalVar_FOCUSS::instance()->cfLambdaESPReSSo_ = cfLambdaESPReSSo_;
		GlobalVar_FOCUSS::instance()->iDimFFT_ = iDimFFT;
		GlobalVar_FOCUSS::instance()->iDimDCTSparse_ = iDimDCTSparse;
		GlobalVar_FOCUSS::instance()->iDimPCASparse_ = iDimPCASparse;
		GlobalVar_FOCUSS::instance()->iDimKernelFFT_ = iDimKernelFFT;
		GlobalVar_FOCUSS::instance()->iTransformFFTBA_ = iTransformFFTBA;
		GlobalVar_FOCUSS::instance()->ikSpaceOut_ = ikSpaceOut;
	}
	else{
		iNOuter_ = GlobalVar_FOCUSS::instance()->iNOuter_;
		iNInner_ = GlobalVar_FOCUSS::instance()->iNInner_;
		iCGResidual_ = GlobalVar_FOCUSS::instance()->iCGResidual_;
		bESPRActiveCS_ = GlobalVar_FOCUSS::instance()->bESPRActiveCS_;
		iVDMap_ = GlobalVar_FOCUSS::instance()->iVDMap_;
		fFullySampled_ = GlobalVar_FOCUSS::instance()->fFullySampled_;
		cfLambdaESPReSSo_ = GlobalVar_FOCUSS::instance()->cfLambdaESPReSSo_;
		cfLambda_ = GlobalVar_FOCUSS::instance()->cfLambda_;
		iESPReSSoDirection_ = GlobalVar_FOCUSS::instance()->iESPReSSoDirection_;
		fPartialFourierVal_ = GlobalVar_FOCUSS::instance()->fPartialFourierVal_;
	
	}
	
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		GDEBUG("lambda is %f \n", GlobalVar_FOCUSS::instance()->cfLambda_);
		GDEBUG("Lambda ESPReSSo is %f \n", GlobalVar_FOCUSS::instance()->cfLambdaESPReSSo_);
		GDEBUG("Fully Sampled is %f \n", GlobalVar_FOCUSS::instance()->fFullySampled_);
		GDEBUG("CS Acceleration is %f \n", GlobalVar_FOCUSS::instance()->fCSAcc_);
		GDEBUG("bESPRActiveCS is %i \n", GlobalVar_FOCUSS::instance()->bESPRActiveCS_);
		GDEBUG("kSpaceOutDim is %i \n", GlobalVar_FOCUSS::instance()->ikSpaceOut_);
		GDEBUG("transformFftBaDim is %i \n", GlobalVar_FOCUSS::instance()->iTransformFFTBA_);
		GDEBUG("kernelFftDim is %i \n", GlobalVar_FOCUSS::instance()->iDimKernelFFT_);
		GDEBUG("pcaSparseDim is %i \n", GlobalVar_FOCUSS::instance()->iDimPCASparse_);
		GDEBUG("dctSparseDim is %i \n", GlobalVar_FOCUSS::instance()->iDimDCTSparse_);
		GDEBUG("fftSparseDim is %i  \n", GlobalVar_FOCUSS::instance()->iDimFFT_);
		GDEBUG("InnerIterations is %i \n", GlobalVar_FOCUSS::instance()->iNInner_);
		GDEBUG("OuterIterations is %i \n", GlobalVar_FOCUSS::instance()->iNOuter_);
		GDEBUG("CG Residual is %i \n", GlobalVar_FOCUSS::instance()->iCGResidual_);
		GDEBUG("VDMap is %i \n", GlobalVar_FOCUSS::instance()->iVDMap_);
	#else
		GADGET_DEBUG2("lambda is %f \n", GlobalVar_FOCUSS::instance()->cfLambda_);
		GADGET_DEBUG2("Lambda ESPReSSo is %f \n", GlobalVar_FOCUSS::instance()->cfLambdaESPReSSo_);
		GADGET_DEBUG2("Fully Sampled is %f \n", GlobalVar_FOCUSS::instance()->fFullySampled_);
		GADGET_DEBUG2("CS Acceleration is %f \n", GlobalVar_FOCUSS::instance()->fCSAcc_);
		GADGET_DEBUG2("bESPRActiveCS is %i \n", GlobalVar_FOCUSS::instance()->bESPRActiveCS_);
		GADGET_DEBUG2("kSpaceOutDim is %i \n", GlobalVar_FOCUSS::instance()->ikSpaceOut_);
		GADGET_DEBUG2("transformFftBaDim is %i \n", GlobalVar_FOCUSS::instance()->iTransformFFTBA_);
		GADGET_DEBUG2("kernelFftDim is %i \n", GlobalVar_FOCUSS::instance()->iDimKernelFFT_);
		GADGET_DEBUG2("pcaSparseDim is %i \n", GlobalVar_FOCUSS::instance()->iDimPCASparse_);
		GADGET_DEBUG2("dctSparseDim is %i \n", GlobalVar_FOCUSS::instance()->iDimDCTSparse_);
		GADGET_DEBUG2("fftSparseDim is %i  \n", GlobalVar_FOCUSS::instance()->iDimFFT_);
		GADGET_DEBUG2("InnerIterations is %i \n", GlobalVar_FOCUSS::instance()->iNInner_);
		GADGET_DEBUG2("OuterIterations is %i \n", GlobalVar_FOCUSS::instance()->iNOuter_);
		GADGET_DEBUG2("CG Residual is %i \n", GlobalVar_FOCUSS::instance()->iCGResidual_);
		GADGET_DEBUG2("VDMap is %i \n", GlobalVar_FOCUSS::instance()->iVDMap_);
	#endif

	if (iNInner_ <= 0) iNInner_ = 20;
	if (iNOuter_ <= 0) iNOuter_ = 2;	

	// p-value for the lp-norm
	fP_ = .5;

	// convergence boundary
	fEpsilon_ = (float)1e-6;

	// setup of the transformation parameters - sparsity dim, fft dim, ..
	fSetupTransformation();

	return GADGET_OK;
};


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
