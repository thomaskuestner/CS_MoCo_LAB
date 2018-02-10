#include "CS_LAB.h"

using namespace Gadgetron;

void CS_LAB::fExternalControl(){

	// algorithm: 0 - FOCUSS
	if (iAlgorithm_ == 0) {

		// evaluate dimension and create suitable class object
		if (iDataset_ == 2 && iTime_ == 0){
			opCS_ = new CS_FOCUSS_2D();
		}
		else if (iDataset_ == 2 && iTime_){
		}
		else if (iDataset_ == 3){
			opCS_ = new CS_FOCUSS_3D();
		}
		else if (iDataset_ == 4){
			opCS_ = new CS_FOCUSS_4D();
		}
	}
	else {
		// not implemented in this version
	}
	
	//opCS_->iCGResidual_					= iCGResidual_;
	opCS_->iNChannels_					= iNChannels_;
	opCS_->fP_							= fP_;
	opCS_->fEpsilon_					= fEpsilon_;		
	opCS_->fCSAccel_					= fCSAccel_;
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
			}
		}
	}
	// check DCT entry
	if (iDCT_Sparse_ != 0){
		for(int i = 0; i < 7; i++){
			int bit = (iDCT_Sparse_ & (1 << i)) >> i;
			if (bit == 1) {
				opCS_->Transform_KernelTransform_->set_transformation_sparsity(1,i);
			}
		}
	}

	// configure KernelTransformation - FFT
	if (iKernel_FFT_dim_ != 0){
		for(int i = 0; i < 7; i++){
			int bit = (iKernel_FFT_dim_ & (1 << i)) >> i;
			if (bit == 1){
				int bit2 = (GlobalVar::instance()->iScrambleDim_ & (1 << i)) >> i;
				if (bit2 == 1){
					opCS_->Transform_KernelTransform_->set_transformation_fft(i, true);
				}
				else{
					opCS_->Transform_KernelTransform_->set_transformation_fft(i, false);
				}			
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
				opCS_->Transform_fftBA_->set_transformation_fft(i, true);
			}
		}
		opCS_->Transform_fftBA_->set_active();
	}

	/*-------------------------------------------------------------------------
	---------------------------------- fftAA ----------------------------------
	--------------------------------------------------------------------------*/
	if (kSpaceOut_ == true){
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,0);
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,1);
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,2);
		opCS_->Transform_fftAA_->set_transformation_fft(0, true);
		opCS_->Transform_fftAA_->set_transformation_fft(1, true);
		opCS_->Transform_fftAA_->set_transformation_fft(2, true);
		opCS_->Transform_fftAA_->set_active();
	}
}

//int CS_LAB::process_config(ACE_Message_Block* mb){
//#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
//		GDEBUG("process config..\n");
//	#else
//		GADGET_DEBUG1("process config..\n");
//	#endif
//	//bXMLControl_ = true;
//	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
//		bXMLControl_ = bXMLControl.value();
//	#else
//		bXMLControl_ = this->get_int_value("bXMLControl");
//	#endif
//
//	if (bXMLControl_) {
//
//		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
//			GDEBUG("XML Control enabled..\n");
//			iNOuter_ = OuterIterations.value();
//		  	iNInner_ = InnerIterations.value();
//			bESPRActiveCS_ = CSESPReSSo.value();
//			cfLambda_ = lambda.value();
//			cfLambdaESPReSSo_ = lambdaESPReSSo.value();
//			int iDimFFT = fftSparseDim.value();
//			int iDimDCTSparse = dctSparseDim.value();
//			int iDimPCASparse = pcaSparseDim.value();
//			int iDimKernelFFT = kernelFftDim.value();
//			int iScrambleDim = scrambleDim.value();
//			int iTransformFFTBA = transformFftBaDim.value();
//			int ikSpaceOut = kSpaceOutDim.value();
//		#else
//			GADGET_DEBUG1("XML Control enabled..\n");
//		  	iNOuter_ = this->get_int_value("OuterIterations");
//		  	iNInner_ = this->get_int_value("InnerIterations");
//			bESPRActiveCS_ = this->get_int_value("CSESPReSSo");
//			cfLambda_ = std::complex<float>(this->get_double_value("lambda"),0.0);
//			cfLambdaESPReSSo_ =  std::complex<float>(this->get_double_value("lambdaESPReSSo"),0.0);
//			int iDimFFT = this->get_int_value("fftSparseDim");
//			int iDimDCTSparse = this->get_int_value("dctSparseDim");
//			int iDimPCASparse = this->get_int_value("pcaSparseDim");
//			int iDimKernelFFT = this->get_int_value("kernelFftDim");
//			int iScrambleDim = this->get_int_value("scrambleDim");
//			int iTransformFFTBA = this->get_int_value("transformFftBaDim");
//			int ikSpaceOut = this->get_int_value("kSpaceOutDim");
//		#endif
//
//		// update global parameters
//		GlobalVar::instance()->iNOuter_ = iNOuter_;
//		GlobalVar::instance()->iNInner_ = iNInner_;
//		GlobalVar::instance()->bESPRActiveCS_ = bESPRActiveCS_;
//		GlobalVar::instance()->cfLambda_ = cfLambda_;
//		GlobalVar::instance()->cfLambdaESPReSSo_ = cfLambdaESPReSSo_;
//		GlobalVar::instance()->iDimFFT_ = iDimFFT;
//		GlobalVar::instance()->iDimDCTSparse_ = iDimDCTSparse;
//		GlobalVar::instance()->iDimPCASparse_ = iDimPCASparse;
//		GlobalVar::instance()->iDimKernelFFT_ = iDimKernelFFT;
//		GlobalVar::instance()->iScrambleDim_ = iScrambleDim;
//		GlobalVar::instance()->iTransformFFTBA_ = iTransformFFTBA;
//		GlobalVar::instance()->ikSpaceOut_ = ikSpaceOut;
//	}
//	else{
//		iNOuter_ = GlobalVar::instance()->iNOuter_;
//		iNInner_ = GlobalVar::instance()->iNInner_;
//		bESPRActiveCS_ = GlobalVar::instance()->bESPRActiveCS_;
//		iVDMap_ = GlobalVar::instance()->iVDMap_;
//		fFullySampled_ = GlobalVar::instance()->fFullySampled_;
//		cfLambdaESPReSSo_ = GlobalVar::instance()->cfLambdaESPReSSo_;
//		cfLambda_ = GlobalVar::instance()->cfLambda_;
//		iESPReSSoDirection_ = GlobalVar::instance()->iESPReSSoDirection_;
//		fPartialFourierVal_ = GlobalVar::instance()->fPartialFourierVal_;
//	}
//
//	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
//		GDEBUG("lambda is %f \n", GlobalVar::instance()->cfLambda_.real());
//		GDEBUG("Lambda ESPReSSo is %f \n", GlobalVar::instance()->cfLambdaESPReSSo_.real());
//		GDEBUG("Fully Sampled is %f \n", GlobalVar::instance()->fFullySampled_);
//		GDEBUG("bESPRActiveCS is %i \n", GlobalVar::instance()->bESPRActiveCS_);
//		GDEBUG("kSpaceOutDim is %i \n", GlobalVar::instance()->ikSpaceOut_);
//		GDEBUG("transformFftBaDim is %i \n", GlobalVar::instance()->iTransformFFTBA_);
//		GDEBUG("kernelFftDim is %i \n", GlobalVar::instance()->iDimKernelFFT_);
//		GDEBUG("pcaSparseDim is %i \n", GlobalVar::instance()->iDimPCASparse_);
//		GDEBUG("dctSparseDim is %i \n", GlobalVar::instance()->iDimDCTSparse_);
//		GDEBUG("fftSparseDim is %i  \n", GlobalVar::instance()->iDimFFT_);
//		GDEBUG("scrambleDim is %i \n", GlobalVar::instance()->iScrambleDim_);
//		GDEBUG("InnerIterations is %i \n", GlobalVar::instance()->iNInner_);
//		GDEBUG("OuterIterations is %i \n", GlobalVar::instance()->iNOuter_);
//	#else
//		GADGET_DEBUG2("lambda is %f \n", GlobalVar::instance()->cfLambda_.real());
//		GADGET_DEBUG2("Lambda ESPReSSo is %f \n", GlobalVar::instance()->cfLambdaESPReSSo_.real());
//		GADGET_DEBUG2("Fully Sampled is %f \n", GlobalVar::instance()->fFullySampled_);
//		GADGET_DEBUG2("bESPRActiveCS is %i \n", GlobalVar::instance()->bESPRActiveCS_);
//		GADGET_DEBUG2("kSpaceOutDim is %i \n", GlobalVar::instance()->ikSpaceOut_);
//		GADGET_DEBUG2("transformFftBaDim is %i \n", GlobalVar::instance()->iTransformFFTBA_);
//		GADGET_DEBUG2("kernelFftDim is %i \n", GlobalVar::instance()->iDimKernelFFT_);
//		GADGET_DEBUG2("pcaSparseDim is %i \n", GlobalVar::instance()->iDimPCASparse_);
//		GADGET_DEBUG2("dctSparseDim is %i \n", GlobalVar::instance()->iDimDCTSparse_);
//		GADGET_DEBUG2("fftSparseDim is %i  \n", GlobalVar::instance()->iDimFFT_);
//		GADGET_DEBUG2("scrambleDim is %i \n", GlobalVar::instance()->iScrambleDim_);
//		GADGET_DEBUG2("InnerIterations is %i \n", GlobalVar::instance()->iNInner_);
//		GADGET_DEBUG2("OuterIterations is %i \n", GlobalVar::instance()->iNOuter_);
//	#endif
//
//	if (GlobalVar::instance()->iNInner_ <= 0) GlobalVar::instance()->iNInner_ = 20;
//	if (GlobalVar::instance()->iNOuter_ <= 0) GlobalVar::instance()->iNOuter_ = 2;
//
//	// p-value for the lp-norm
//	fP_ = .5;
//
//	// convergence boundary
//	fEpsilon_ = (float)1e-6;
//
//	// setup of the transformation parameters - sparsity dim, fft dim, ..
//	fSetupTransformation();
//
//	return GADGET_OK;
//};


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
		
		opCS_ = new CS_FOCUSS_4D();

		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Incoming data is 4D - starting 4D FOCUSS reconstruction\n");
		#else
			GADGET_DEBUG1("Incoming data is 4D - starting 4D FOCUSS reconstruction\n");
		#endif
	}

	// set parameters of the FOCUSS class - required, because the xml config file is read in by CS_CONTROL class and not by FOCUSS class
	//opCS_->iCGResidual_					= iCGResidual_;
	opCS_->iNChannels_					= iNChannels_;
	opCS_->fP_							= fP_;
	opCS_->fEpsilon_					= fEpsilon_;
	opCS_->fCSAccel_					= fCSAccel_;
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
