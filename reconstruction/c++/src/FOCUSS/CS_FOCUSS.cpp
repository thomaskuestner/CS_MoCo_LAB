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
#include "GlobalVar_FOCUSS.h"
#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron{

// read the XML configuration parameters
int CS_FOCUSS::process_config(ACE_Message_Block* mb){

	// how to calculate the beta value
	iCGResidual_ = this->get_int_value("CG Beta");

	// maximum number of FOCUSS iterations
	iNOuter_ = this->get_int_value("OuterIterations");
	if (iNOuter_ <= 0) iNOuter_ = 2;

	// maximum number of CG iterations
	iNInner_ = this->get_int_value("InnerIterations");
	if (iNInner_ <= 0) iNInner_ = 20;

	// p-value for the lp-norm
	fP_ = .5;

	// use ESPReSSo-constraint for pure CS data
	bESPRActiveCS_ = this->get_bool_value("CS - ESPReSSo");

	// convergence boundary
	fEpsilon_ = (float)1e-6;

	// setup of the transformation parameters - sparsity dim, fft dim, ..
	fSetupTransformation();

	return GADGET_OK;
};

// set several variables
void CS_FOCUSS::fInitVal(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1){

	// initialize the global vector variables
	for(int i = 0; i < 20; i++){
		GlobalVar_FOCUSS::instance()->vbStatPrinc_.push_back(false);
		GlobalVar_FOCUSS::instance()->vfPrincipleComponents_.push_back(new hoNDArray<std::complex<float> > ());
	}

	iVDMap_				= m1->getObjectPtr()->user_int[7];
	fFullySampled_		= m1->getObjectPtr()->user_float[5];
	cfLambdaESPReSSo_	= m1->getObjectPtr()->user_float[6];
	cfLambda_			= m1->getObjectPtr()->user_float[7];
	iESPReSSoDirection_ = m1->getObjectPtr()->user_int[7];
	fPartialFourierVal_ = m1->getObjectPtr()->user_float[3];
}

void CS_FOCUSS::fSetupTransformation(){

	// instantiate transformation objects
	Transform_KernelTransform_	    = new Transform();
	Transform_fftBA_				= new Transform();
	Transform_fftAA_				= new Transform();

	int dim;
	/*-------------------------------------------------------------------------
	---------------------------- KernelTransform ------------------------------
	--------------------------------------------------------------------------*/
	// KernelTransform - active
	Transform_KernelTransform_->set_active();

	// configure KernelTransformation - sparsifying transform
	// check FFT entry
	if ((dim = this->get_int_value("FFT_Sparse")) != 0){
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
	if ((dim = this->get_int_value("DCT_Sparse")) != 0){
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
	/*if ((dim = this->get_int_value("PCA_Sparse")) != 0){
		for(int i = 0; i < 7; i++){
			int bit = (dim & (1 << i)) >> i;
			if (bit == 1){
				Transform_KernelTransform_->set_transformation_sparsity(2,i);
				GADGET_DEBUG2("KernelTransform - PCA sparse - dim: %i \n", i);
			}
		}
	}*/

	// configure KernelTransformation - FFT
	if ((dim = this->get_int_value("Kernel_FFT_dim")) != 0){
		for(int i = 0; i < 7; i++){
			int bit = (dim & (1 << i)) >> i;
			if (bit == 1){
				Transform_KernelTransform_->set_transformation_fft(i);
				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("KernelTransform - FFT - dim: %i \n", i);
				#else
					GADGET_DEBUG2("KernelTransform - FFT - dim: %i \n", i);
				#endif
			}
		}
	}

	/*-------------------------------------------------------------------------
	---------------------------------- fftBA ----------------------------------
	--------------------------------------------------------------------------*/
	// configure fftBA - transform dimension before start FOCUSS
	if ((dim = this->get_int_value("Transform_fftBA_dim")) != 0){
		for(int i = 0; i < 7; i++){
			int bit = (dim & (1 << i)) >> i;
			if (bit == 1){
				Transform_fftBA_->set_transformation_fft(i);
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
	if ((dim = this->get_int_value("kSpaceOut")) != 0){
		Transform_fftAA_->set_transformation_sparsity(0,0);
		Transform_fftAA_->set_transformation_sparsity(0,1);
		Transform_fftAA_->set_transformation_sparsity(0,2);
		Transform_fftAA_->set_transformation_fft(0);
		Transform_fftAA_->set_transformation_fft(1);
		Transform_fftAA_->set_transformation_fft(2);
		Transform_fftAA_->set_active();
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Transform_fftAA_ - active: %i \n", Transform_fftAA_->get_active());
		#else
			GADGET_DEBUG2("Transform_fftAA_ - active: %i \n", Transform_fftAA_->get_active());
		#endif
	}
}

}
