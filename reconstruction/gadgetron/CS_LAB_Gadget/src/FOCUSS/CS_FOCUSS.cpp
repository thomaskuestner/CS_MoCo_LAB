/*
file name	: 	CS_FOCUSS.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.1

date		: 	13.02.2015

description	: 	implementation of the class "CS_FOCUSS.h"

references	:	-

notes		:	methods, which are included in this file, are same for each inherited sub-class
*/

#include "CS_FOCUSS.h"

namespace Gadgetron{

// read the XML configuration parameters
int CS_FOCUSS::process_config(ACE_Message_Block* mb){
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
			GlobalVar::instance()->iNOuter_ = OuterIterations.value();
			GlobalVar::instance()->iNInner_ = InnerIterations.value();
			GlobalVar::instance()->bESPRActiveCS_ = CSESPReSSo.value();
			GlobalVar::instance()->cfLambda_ = lambda.value();
			GlobalVar::instance()->cfLambdaESPReSSo_ = lambdaESPReSSo.value();
			int iDimFFT = fftSparseDim.value();
			int iDimDCTSparse = dctSparseDim.value();
			int iDimPCASparse = pcaSparseDim.value();
			int iDimKernelFFT = kernelFftDim.value();
			int iScrambleDim = scrambleDim.value();
			int iTransformFFTBA = transformFftBaDim.value();
			int ikSpaceOut = kSpaceOutDim.value();
			int iNorm_ = norm.value();
		#else
			GADGET_DEBUG1("XML Control enabled..\n");
			GlobalVar::instance()->iNOuter_ = this->get_int_value("OuterIterations");
			GlobalVar::instance()->iNInner_ = this->get_int_value("InnerIterations");
			GlobalVar::instance()->bESPRActiveCS_ = this->get_int_value("CSESPReSSo");
			GlobalVar::instance()->cfLambda_ = this->get_double_value("lambda");
			GlobalVar::instance()->cfLambdaESPReSSo_ = this->get_double_value("lambdaESPReSSo");
			int iDimFFT = this->get_int_value("fftSparseDim");
			int iDimDCTSparse = this->get_int_value("dctSparseDim");
			int iDimPCASparse = this->get_int_value("pcaSparseDim");
			int iDimKernelFFT = this->get_int_value("kernelFftDim");
			int iScrambleDim = this->get_int_value("scrambleDim");
			int iTransformFFTBA = this->get_int_value("transformFftBaDim");
			int ikSpaceOut = this->get_int_value("kSpaceOutDim");
			iNorm_ = this->get_int_value("norm");
		#endif

		// update global parameters
		GlobalVar::instance()->iDimFFT_ = iDimFFT;
		GlobalVar::instance()->iDimDCTSparse_ = iDimDCTSparse;
		GlobalVar::instance()->iDimPCASparse_ = iDimPCASparse;
		GlobalVar::instance()->iDimKernelFFT_ = iDimKernelFFT;
		GlobalVar::instance()->iScrambleDim_ = iScrambleDim;
		GlobalVar::instance()->iTransformFFTBA_ = iTransformFFTBA;
		GlobalVar::instance()->ikSpaceOut_ = ikSpaceOut;
	}
	
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		GDEBUG("lambda is %f \n", GlobalVar::instance()->cfLambda_.real());
		GDEBUG("Lambda ESPReSSo is %f \n", GlobalVar::instance()->cfLambdaESPReSSo_.real());
		GDEBUG("Fully Sampled is %f \n", GlobalVar::instance()->fFullySampled_);
		GDEBUG("bESPRActiveCS is %i \n", GlobalVar::instance()->bESPRActiveCS_);
		GDEBUG("kSpaceOutDim is %i \n", GlobalVar::instance()->ikSpaceOut_);
		GDEBUG("transformFftBaDim is %i \n", GlobalVar::instance()->iTransformFFTBA_);
		GDEBUG("kernelFftDim is %i \n", GlobalVar::instance()->iDimKernelFFT_);
		GDEBUG("pcaSparseDim is %i \n", GlobalVar::instance()->iDimPCASparse_);
		GDEBUG("dctSparseDim is %i \n", GlobalVar::instance()->iDimDCTSparse_);
		GDEBUG("fftSparseDim is %i  \n", GlobalVar::instance()->iDimFFT_);
		GDEBUG("scrambleDim is %i \n", GlobalVar::instance()->iScrambleDim_);
		GDEBUG("InnerIterations is %i \n", GlobalVar::instance()->iNInner_);
		GDEBUG("OuterIterations is %i \n", GlobalVar::instance()->iNOuter_);
	#else
		GADGET_DEBUG2("lambda is %f \n", GlobalVar::instance()->cfLambda_);
		GADGET_DEBUG2("Lambda ESPReSSo is %f \n", GlobalVar::instance()->cfLambdaESPReSSo_);
		GADGET_DEBUG2("Fully Sampled is %f \n", GlobalVar::instance()->fFullySampled_);
		GADGET_DEBUG2("bESPRActiveCS is %i \n", GlobalVar::instance()->bESPRActiveCS_);
		GADGET_DEBUG2("kSpaceOutDim is %i \n", GlobalVar::instance()->ikSpaceOut_);
		GADGET_DEBUG2("transformFftBaDim is %i \n", GlobalVar::instance()->iTransformFFTBA_);
		GADGET_DEBUG2("kernelFftDim is %i \n", GlobalVar::instance()->iDimKernelFFT_);
		GADGET_DEBUG2("pcaSparseDim is %i \n", GlobalVar::instance()->iDimPCASparse_);
		GADGET_DEBUG2("dctSparseDim is %i \n", GlobalVar::instance()->iDimDCTSparse_);
		GADGET_DEBUG2("fftSparseDim is %i  \n", GlobalVar::instance()->iDimFFT_);
		GADGET_DEBUG2("scrambleDim is %i \n", GlobalVar::instance()->iScrambleDim_);
		GADGET_DEBUG2("InnerIterations is %i \n", GlobalVar::instance()->iNInner_);
		GADGET_DEBUG2("OuterIterations is %i \n", GlobalVar::instance()->iNOuter_);
	#endif

	if (GlobalVar::instance()->iNInner_ <= 0) GlobalVar::instance()->iNInner_ = 20;
	if (GlobalVar::instance()->iNOuter_ <= 0) GlobalVar::instance()->iNOuter_ = 2;

	// p-value for the lp-norm
	fP_ = .5;

	// convergence boundary
	fEpsilon_ = (float)1e-6;

	// setup of the transformation parameters - sparsity dim, fft dim, ..
	fSetupTransformation();

	return GADGET_OK;
};

// set several variables
void CS_FOCUSS::fInitVal(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1){

	// initialize the global vector variables
	if (GlobalVar::instance()->iDimPCASparse_ != 0){
		for (int iI = 0; iI < 6; iI++){
			GlobalVar::instance()->vbStatPrinc_.push_back(false);
			GlobalVar::instance()->KLTVec_.push_back(new hoNDKLT_CS< std::complex< float > > ());
		}
	}
}

// init transformations
void CS_FOCUSS::fSetupTransformation(){

	// instantiate transformation objects
	Transform_KernelTransform_	= new Transform();
	Transform_fftBA_		= new Transform();
	Transform_fftAA_		= new Transform();

	int dim;
	/*-------------------------------------------------------------------------
	---------------------------- KernelTransform ------------------------------
	--------------------------------------------------------------------------*/
	// KernelTransform - active
	Transform_KernelTransform_->set_active();

	// configure KernelTransformation - sparsifying transform
	// check FFT entry
	dim = GlobalVar::instance()->iDimFFT_;
	if (dim != 0){
		for(int i = 0; i < 7; i++){
			int bit = (dim & (1 << i)) >> i;
			if (bit == 1) {
				Transform_KernelTransform_->set_transformation_sparsity(0,i);

				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("KernelTransform - FFT sparse - dim: %i \n", i);
				#else
					GADGET_DEBUG2("KernelTransform - FFT sparse - dim: %i \n", i);
				#endif
			}
		}
	}

	// check DCT entry
	dim = GlobalVar::instance()->iDimDCTSparse_;
	if (dim != 0){
		for(int i = 0; i < 7; i++){
			int bit = (dim & (1 << i)) >> i;
			if (bit == 1) {
				Transform_KernelTransform_->set_transformation_sparsity(1,i);
				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("KernelTransform - DCT sparse - dim: %i \n", i);
				#else
					GADGET_DEBUG2("KernelTransform - DCT sparse - dim: %i \n", i);
				#endif
			}
		}
	}

	// check PCA entry 
	dim = GlobalVar::instance()->iDimPCASparse_;
	if (dim != 0)
	{
		for(int i = 0; i < 7; i++){
			int bit = (dim & (1 << i)) >> i;
			if (bit == 1){
				Transform_KernelTransform_->set_transformation_sparsity(2,i);
				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("KernelTransform - PCA sparse - dim: %i \n", i);
				#else
					GADGET_DEBUG2("KernelTransform - PCA sparse - dim: %i \n", i);
				#endif				
			}
		}
	}	

	// configure KernelTransformation - FFT
	dim = GlobalVar::instance()->iDimKernelFFT_;
	if (dim != 0){
		for(int i = 0; i < 7; i++){
			int bit = (dim & (1 << i)) >> i;
			if (bit == 1){
				int bit2 = (GlobalVar::instance()->iScrambleDim_ & (1 << i)) >> i;
				if (bit2 == 1){
					Transform_KernelTransform_->set_transformation_fft(i, true);
				}
				else{
					Transform_KernelTransform_->set_transformation_fft(i, false);
				}	
				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("KernelTransform - FFT - dim: %i - scrambling: %i\n", i, bit2);					
				#else
					GADGET_DEBUG2("KernelTransform - FFT - dim: %i - scrambling: %i\n", i, bit2);
				#endif
			}
		}
	}

	/*-------------------------------------------------------------------------
	---------------------------------- fftBA ----------------------------------
	--------------------------------------------------------------------------*/
	// configure fftBA - transform dimension before start FOCUSS
	dim = GlobalVar::instance()->iTransformFFTBA_;
	if (dim != 0){
		for(int i = 0; i < 7; i++){
			int bit = (dim & (1 << i)) >> i;
			if (bit == 1){
				Transform_fftBA_->set_transformation_fft(i, true);
				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("fftBA - dim: %i \n", i);
				#else
					GADGET_DEBUG2("fftBA - dim: %i \n", i);
				#endif
			}
		}
		Transform_fftBA_->set_active();
	}

	// configure fftAA - output image or k-space
	dim = GlobalVar::instance()->ikSpaceOut_;
	if (dim != 0){
		Transform_fftAA_->set_transformation_sparsity(0,0);
		Transform_fftAA_->set_transformation_sparsity(0,1);
		Transform_fftAA_->set_transformation_sparsity(0,2);
		Transform_fftAA_->set_transformation_fft(0, true);
		Transform_fftAA_->set_transformation_fft(1, true);
		Transform_fftAA_->set_transformation_fft(2, true);
		Transform_fftAA_->set_active();
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Transform_fftAA_ - active: %i \n", Transform_fftAA_->get_active());
		#else
			GADGET_DEBUG2("Transform_fftAA_ - active: %i \n", Transform_fftAA_->get_active());
		#endif
	}
}
}

