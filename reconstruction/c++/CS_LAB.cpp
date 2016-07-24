#include "CS_LAB.h"

#include "GadgetIsmrmrdReadWrite.h"
#include "mex.h"

using namespace Gadgetron;

void CS_LAB::fExternalControl(){

	// algorithm: 0 - FOCUSS
	if (iAlgorithm_ == 0) {

		// evaluate dimension and create suitable class object
		if (iDataset_ == 2 && iTime_ == 0){
			opCS_ = new CS_FOCUSS_2D();
			if (bMatlab_)
					mexPrintf("2D dataset detected - init 2D FOCUSS..\n");mexEvalString("drawnow;");
		}
		else if (iDataset_ == 2 && iTime_){
			//opCS_ = new CS_FOCUSS_2Dt();
			if (bMatlab_)
					mexPrintf("2Dt dataset detected - init 2Dt FOCUSS..\n");mexEvalString("drawnow;");
		}
		else if (iDataset_ == 3){
			opCS_ = new CS_FOCUSS_3D();
			if (bMatlab_)
					mexPrintf("3D dataset detected - init 3D FOCUSS..\n");mexEvalString("drawnow;");
		}
	}
	else {
		// not implemented in this version
	}
	
	opCS_->iCGResidual_					= iCGResidual_;
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
				if (bMatlab_)
					mexPrintf("KernelTransform - FFT sparse - dim: %i \n", i); mexEvalString("drawnow;");
			}
		}
	}
	// check DCT entry
	if (iDCT_Sparse_ != 0){
		for(int i = 0; i < 7; i++){
			int bit = (iDCT_Sparse_ & (1 << i)) >> i;
			if (bit == 1) {
				opCS_->Transform_KernelTransform_->set_transformation_sparsity(1,i);
				if (bMatlab_)
					mexPrintf("KernelTransform - DCT sparse - dim: %i \n", i); mexEvalString("drawnow;");
			}
		}
	}

	// configure KernelTransformation - FFT
	if (iKernel_FFT_dim_ != 0){
		for(int i = 0; i < 7; i++){
			int bit = (iKernel_FFT_dim_ & (1 << i)) >> i;
			if (bit == 1){
				opCS_->Transform_KernelTransform_->set_transformation_fft(i);
				if (bMatlab_)
					mexPrintf("KernelTransform - FFT - dim: %i \n", i); mexEvalString("drawnow;");
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
				if (bMatlab_)
					mexPrintf("KernelTransform - fftBA - dim: %i \n", i); mexEvalString("drawnow;");
			}
		}
		opCS_->Transform_fftBA_->set_active();
	}

	/*-------------------------------------------------------------------------
	---------------------------------- fftAA ----------------------------------
	--------------------------------------------------------------------------*/
	if (kSpaceOut_ == true){
		if (bMatlab_)
					mexPrintf("output kSpace data: true\n"); mexEvalString("drawnow;");
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,0);
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,1);
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,2);
		opCS_->Transform_fftAA_->set_transformation_fft(0);
		opCS_->Transform_fftAA_->set_transformation_fft(1);
		opCS_->Transform_fftAA_->set_transformation_fft(2);
		opCS_->Transform_fftAA_->set_active();
	}
	else{
		if (bMatlab_)
					mexPrintf("output kSpace data: false\n"); mexEvalString("drawnow;");
	}
}

GADGET_FACTORY_DECLARE(CS_LAB)
