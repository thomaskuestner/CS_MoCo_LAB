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
/*
int CS_LAB::process_config(ACE_Message_Block* mb){

	
};*/


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
