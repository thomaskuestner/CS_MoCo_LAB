#include "CS_LAB.h"

#include "GadgetIsmrmrdReadWrite.h"

using namespace Gadgetron;

void CS_LAB::fExternalControl(){

	// algorithm: 0 - FOCUSS
	if (iAlgorithm_ == 0) {

		// evaluate dimension and create suitable class object
		if (iDataset_ == 2 && iTime_ == 0){
			opCS_ = new CS_FOCUSS_2D();
			// if (bMatlab_)
					// mexPrintf("2D dataset detected - init 2D FOCUSS..\n");mexEvalString("drawnow;");
		}
		else if (iDataset_ == 2 && iTime_){
			//opCS_ = new CS_FOCUSS_2Dt();
			// if (bMatlab_)
					// mexPrintf("2Dt dataset detected - init 2Dt FOCUSS..\n");mexEvalString("drawnow;");
		}
		else if (iDataset_ == 3){
			opCS_ = new CS_FOCUSS_3D();
			// if (bMatlab_)
					// mexPrintf("3D dataset detected - init 3D FOCUSS..\n");mexEvalString("drawnow;");
		}
	}
	else {
		// not implemented in this version
	}
	
	//opCS_->iCGResidual_					= iCGResidual_;
	opCS_->iNChannels_					= iNChannels_;
	opCS_->iNOuter_						= iNOuter_;
	opCS_->iNInner_						= iNInner_;
	opCS_->fP_							= fP_;
	opCS_->cfLambda_					= cfLambda_;
	opCS_->cfLambdaESPReSSo_			= cfLambdaESPReSSo_;
	opCS_->fEpsilon_					= fEpsilon_;		
	opCS_->fCSAccel_					= fCSAccel_;
	opCS_->iESPReSSoDirection_			= iESPReSSoDirection_;
	opCS_->fPartialFourierVal_			= fPartialFourierVal_;
	opCS_->fFullySampled_				= fFullySampled_;
	opCS_->bESPRActiveCS_				= bESPRActiveCS_;
	opCS_->hacfFilter_					= hacfFilter_;
	opCS_->Transform_KernelTransform_	= Transform_KernelTransform_;
	opCS_->Transform_fftBA_				= Transform_fftBA_;
	opCS_->Transform_fftAA_				= Transform_fftAA_;

	// disable standalone Gadget behaviour
	opCS_->bControl_		= false;
	opCS_->bDebug_			= true;
	opCS_->bMatlab_			= true;
	
	// instantiate transformation objects
	opCS_->Transform_KernelTransform_	= new Transform();
	opCS_->Transform_fftAA_				= new Transform();
	opCS_->Transform_fftBA_				= new Transform();

	/*-------------------------------------------------------------------------
	---------------------------- KernelTransform ------------------------------
	--------------------------------------------------------------------------*/
	// KernelTransform - active
	opCS_->Transform_KernelTransform_->set_active();

	// configure KernelTransformation - sparsifying transform
	// check FFT entry
	if (iFFT_Sparse_ != 0){
		for(int i = 0; i < 7; i++){
			int bit = (iFFT_Sparse_ & (1 << i)) >> i;
			if (bit == 1) {
				opCS_->Transform_KernelTransform_->set_transformation_sparsity(0,i);
				// if (bMatlab_)
					// mexPrintf("KernelTransform - FFT sparse - dim: %i \n", i); mexEvalString("drawnow;");
			}
		}
	}
	// check DCT entry
	if (iDCT_Sparse_ != 0){
		for(int i = 0; i < 7; i++){
			int bit = (iDCT_Sparse_ & (1 << i)) >> i;
			if (bit == 1) {
				opCS_->Transform_KernelTransform_->set_transformation_sparsity(1,i);
				// if (bMatlab_)
					// mexPrintf("KernelTransform - DCT sparse - dim: %i \n", i); mexEvalString("drawnow;");
			}
		}
	}

	// configure KernelTransformation - FFT
	if (iKernel_FFT_dim_ != 0){
		for(int i = 0; i < 7; i++){
			int bit = (iKernel_FFT_dim_ & (1 << i)) >> i;
			if (bit == 1){
				opCS_->Transform_KernelTransform_->set_transformation_fft(i);
				// if (bMatlab_)
					// mexPrintf("KernelTransform - FFT - dim: %i \n", i); mexEvalString("drawnow;");
			}
		}
	}

	/*-------------------------------------------------------------------------
	---------------------------------- fftBA ----------------------------------
	--------------------------------------------------------------------------*/
	// configure fftBA - transform dimension before start FOCUSS
	if (iFFTBA_!= 0){
		for(int i = 0; i < 7; i++){
			int bit = (iFFTBA_ & (1 << i)) >> i;
			if (bit == 1){
				opCS_->Transform_fftBA_->set_transformation_fft(i);
				// if (bMatlab_)
					// mexPrintf("KernelTransform - fftBA - dim: %i \n", i); mexEvalString("drawnow;");
			}
		}
		opCS_->Transform_fftBA_->set_active();
	}

	/*-------------------------------------------------------------------------
	---------------------------------- fftAA ----------------------------------
	--------------------------------------------------------------------------*/
	if (kSpaceOut_ == true){
		// if (bMatlab_)
					// mexPrintf("output kSpace data: true\n"); mexEvalString("drawnow;");
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,0);
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,1);
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,2);
		opCS_->Transform_fftAA_->set_transformation_fft(0);
		opCS_->Transform_fftAA_->set_transformation_fft(1);
		opCS_->Transform_fftAA_->set_transformation_fft(2);
		opCS_->Transform_fftAA_->set_active();
	}
	else{
		// if (bMatlab_)
					// mexPrintf("output kSpace data: false\n"); mexEvalString("drawnow;");
	}
}

int CS_LAB::process_config(ACE_Message_Block* mb){
#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		GDEBUG("process config..\n");
	#else
		GADGET_DEBUG1("process config..\n");
	#endif	
	//bXMLControl_ = true;
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		bXMLControl_ = bXMLControl.value();		
	#else
		bXMLControl_ = this->get_int_value("bXMLControl");		
	#endif	
	
	if (bXMLControl_) {

		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("XML Control enabled..\n");		
			iNOuter_ = OuterIterations.value();
		  	iNInner_ = InnerIterations.value();
			bESPRActiveCS_ = CSESPReSSo.value();
			cfLambda_ = std::complex< float > (lambda.value(), 0.0);			
			cfLambdaESPReSSo_ = lambdaESPReSSo.value();
			int iDimFFT = fftSparseDim.value();
			int iDimDCTSparse = dctSparseDim.value();
			int iDimPCASparse = pcaSparseDim.value();
			int iDimKernelFFT = kernelFftDim.value();
			int iTransformFFTBA = transformFftBaDim.value();
			int ikSpaceOut = kSpaceOutDim.value();
		#else
			GADGET_DEBUG1("XML Control enabled..\n");
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
		bESPRActiveCS_ = GlobalVar_FOCUSS::instance()->bESPRActiveCS_;
		iVDMap_ = GlobalVar_FOCUSS::instance()->iVDMap_;
		fFullySampled_ = GlobalVar_FOCUSS::instance()->fFullySampled_;
		cfLambdaESPReSSo_ = GlobalVar_FOCUSS::instance()->cfLambdaESPReSSo_;
		cfLambda_ = GlobalVar_FOCUSS::instance()->cfLambda_;
		iESPReSSoDirection_ = GlobalVar_FOCUSS::instance()->iESPReSSoDirection_;
		fPartialFourierVal_ = GlobalVar_FOCUSS::instance()->fPartialFourierVal_;
	
	}
	
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		GDEBUG("lambda is %f \n", GlobalVar_FOCUSS::instance()->cfLambda_.real());
		GDEBUG("Lambda ESPReSSo is %f \n", GlobalVar_FOCUSS::instance()->cfLambdaESPReSSo_.real());
		GDEBUG("Fully Sampled is %f \n", GlobalVar_FOCUSS::instance()->fFullySampled_);
		GDEBUG("bESPRActiveCS is %i \n", GlobalVar_FOCUSS::instance()->bESPRActiveCS_);
		GDEBUG("kSpaceOutDim is %i \n", GlobalVar_FOCUSS::instance()->ikSpaceOut_);
		GDEBUG("transformFftBaDim is %i \n", GlobalVar_FOCUSS::instance()->iTransformFFTBA_);
		GDEBUG("kernelFftDim is %i \n", GlobalVar_FOCUSS::instance()->iDimKernelFFT_);
		GDEBUG("pcaSparseDim is %i \n", GlobalVar_FOCUSS::instance()->iDimPCASparse_);
		GDEBUG("dctSparseDim is %i \n", GlobalVar_FOCUSS::instance()->iDimDCTSparse_);
		GDEBUG("fftSparseDim is %i  \n", GlobalVar_FOCUSS::instance()->iDimFFT_);
		GDEBUG("InnerIterations is %i \n", GlobalVar_FOCUSS::instance()->iNInner_);
		GDEBUG("OuterIterations is %i \n", GlobalVar_FOCUSS::instance()->iNOuter_);
	#else
		GADGET_DEBUG2("lambda is %f \n", GlobalVar_FOCUSS::instance()->cfLambda_);
		GADGET_DEBUG2("Lambda ESPReSSo is %f \n", GlobalVar_FOCUSS::instance()->cfLambdaESPReSSo_);
		GADGET_DEBUG2("Fully Sampled is %f \n", GlobalVar_FOCUSS::instance()->fFullySampled_);
		GADGET_DEBUG2("bESPRActiveCS is %i \n", GlobalVar_FOCUSS::instance()->bESPRActiveCS_);
		GADGET_DEBUG2("kSpaceOutDim is %i \n", GlobalVar_FOCUSS::instance()->ikSpaceOut_);
		GADGET_DEBUG2("transformFftBaDim is %i \n", GlobalVar_FOCUSS::instance()->iTransformFFTBA_);
		GADGET_DEBUG2("kernelFftDim is %i \n", GlobalVar_FOCUSS::instance()->iDimKernelFFT_);
		GADGET_DEBUG2("pcaSparseDim is %i \n", GlobalVar_FOCUSS::instance()->iDimPCASparse_);
		GADGET_DEBUG2("dctSparseDim is %i \n", GlobalVar_FOCUSS::instance()->iDimDCTSparse_);
		GADGET_DEBUG2("fftSparseDim is %i  \n", GlobalVar_FOCUSS::instance()->iDimFFT_);
		GADGET_DEBUG2("InnerIterations is %i \n", GlobalVar_FOCUSS::instance()->iNInner_);
		GADGET_DEBUG2("OuterIterations is %i \n", GlobalVar_FOCUSS::instance()->iNOuter_);
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


int CS_LAB::process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2){
	
	// get dimension of the incoming data object
	std::vector<size_t> vDims = *m2->getObjectPtr()->get_dimensions();

	// copy GadgetContainer and init with m2 data
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* tmp_m2 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
	tmp_m2->getObjectPtr()->create(vDims);
	memcpy(tmp_m2->getObjectPtr()->get_data_ptr(), m2->getObjectPtr()->get_data_ptr(), m2->getObjectPtr()->get_number_of_elements()*sizeof(std::complex< float >));
	
	// evaluate dimension and create suitable class object
	if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) == 1 && vDims.at(3) == 1){
		opCS_ = new CS_FOCUSS_2D();
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Incoming data is 2D - starting 2D FOCUSS reconstruction\n");
		#else
			GADGET_DEBUG1("Incoming data is 2D - starting 2D FOCUSS reconstruction\n");
		#endif
	}
	else if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) == 1 && vDims.at(3) > 1){
		//opCS_ = new CS_FOCUSS_2Dt();
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Incoming data is 2Dt - starting 2Dt FOCUSS reconstruction\n");
		#else
			GADGET_DEBUG1("Incoming data is 2Dt - starting 2Dt FOCUSS reconstruction\n");
		#endif
	}
	else if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) > 1 && vDims.at(3) == 1){
		// squeeze array due to x,y,z,c dimension of 3D FOCUSS class
		sum_dim(*tmp_m2->getObjectPtr(), 3, *tmp_m2->getObjectPtr());
		
		opCS_ = new CS_FOCUSS_3D();
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
	//opCS_->iCGResidual_					= iCGResidual_;
	opCS_->iNChannels_					= iNChannels_;
	opCS_->iNOuter_						= iNOuter_;	
	opCS_->iNInner_						= iNInner_;
	opCS_->fP_							= fP_;
	opCS_->cfLambda_					= cfLambda_;
	opCS_->cfLambdaESPReSSo_			= cfLambdaESPReSSo_;
	opCS_->fEpsilon_					= fEpsilon_;
	opCS_->fCSAccel_					= fCSAccel_;
	opCS_->iESPReSSoDirection_			= iESPReSSoDirection_;
	opCS_->fPartialFourierVal_			= fPartialFourierVal_;
	opCS_->bESPRActiveCS_				= bESPRActiveCS_;
	opCS_->hacfFilter_					= hacfFilter_;
	opCS_->Transform_KernelTransform_	= Transform_KernelTransform_;
	opCS_->Transform_fftBA_				= Transform_fftBA_;
	opCS_->Transform_fftAA_				= Transform_fftAA_;

	// disable standalone Gadget behaviour
	opCS_->bControl_	= true;
	opCS_->bDebug_		= true;
	opCS_->bMatlab_		= false;

	// process data in class member function
	opCS_->process(m1, tmp_m2);

	//Now pass on image
	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}
	return GADGET_OK;
}
GADGET_FACTORY_DECLARE(CS_LAB)
