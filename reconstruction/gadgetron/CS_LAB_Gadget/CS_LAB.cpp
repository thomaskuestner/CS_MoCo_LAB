#include "CS_LAB.h"

using namespace Gadgetron;

void CS_LAB::fExternalControl()
{
	// algorithm: 0 - FOCUSS
	if (iAlgorithm_ == 0) {
		// evaluate dimension and create suitable class object
		if (iDataset_ == 2 && iTime_ == 0) {
			opCS_ = new CS_FOCUSS_2D();
		} else if (iDataset_ == 2 && iTime_) {
// 			opCS_ = new CS_FOCUSS_2Dt();
		} else if (iDataset_ == 3) {
			opCS_ = new CS_FOCUSS_3D();
		} else if (iDataset_ == 4) {
			opCS_ = new CS_FOCUSS_4D();
		}
	} else {
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
	if (iFFT_Sparse_ != 0) {
		for(int i = 0; i < 7; i++) {
			int bit = (iFFT_Sparse_ & (1 << i)) >> i;
			if (bit == 1) {
				opCS_->Transform_KernelTransform_->set_transformation_sparsity(0,i);
			}
		}
	}

	// check DCT entry
	if (iDCT_Sparse_ != 0) {
		for(int i = 0; i < 7; i++) {
			int bit = (iDCT_Sparse_ & (1 << i)) >> i;
			if (bit == 1) {
				opCS_->Transform_KernelTransform_->set_transformation_sparsity(1,i);
			}
		}
	}

	// configure KernelTransformation - FFT
	if (iKernel_FFT_dim_ != 0) {
		for(int i = 0; i < 7; i++) {
			int bit = (iKernel_FFT_dim_ & (1 << i)) >> i;
			if (bit == 1) {
				int bit2 = (GlobalVar::instance()->iScrambleDim_ & (1 << i)) >> i;
				if (bit2 == 1) {
					opCS_->Transform_KernelTransform_->set_transformation_fft(i, true);
				} else {
					opCS_->Transform_KernelTransform_->set_transformation_fft(i, false);
				}
			}
		}
	}

	/*-------------------------------------------------------------------------
	---------------------------------- fftBA ----------------------------------
	--------------------------------------------------------------------------*/
	// configure fftBA - transform dimension before start FOCUSS
	if (iFFTBA_!= 0) {
		for(int i = 0; i < 7; i++) {
			int bit = (iFFTBA_ & (1 << i)) >> i;
			if (bit == 1) {
				opCS_->Transform_fftBA_->set_transformation_fft(i, true);
			}
		}
		opCS_->Transform_fftBA_->set_active();
	}

	/*-------------------------------------------------------------------------
	---------------------------------- fftAA ----------------------------------
	--------------------------------------------------------------------------*/
	if (kSpaceOut_ == true) {
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,0);
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,1);
		opCS_->Transform_fftAA_->set_transformation_sparsity(0,2);
		opCS_->Transform_fftAA_->set_transformation_fft(0, true);
		opCS_->Transform_fftAA_->set_transformation_fft(1, true);
		opCS_->Transform_fftAA_->set_transformation_fft(2, true);
		opCS_->Transform_fftAA_->set_active();
	}
}

int CS_LAB::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	// get dimension of the incoming data object
	std::vector<size_t> vDims = *m2->getObjectPtr()->get_dimensions();

	// copy GadgetContainer and init with m2 data
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* tmp_m2 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
	tmp_m2->getObjectPtr()->create(vDims);
	memcpy(tmp_m2->getObjectPtr()->get_data_ptr(), m2->getObjectPtr()->get_data_ptr(), m2->getObjectPtr()->get_number_of_elements()*sizeof(std::complex< float >));
	
	// evaluate dimension and create suitable class object
	if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) == 1 && vDims.at(3) == 1) {
		opCS_ = new CS_FOCUSS_2D();
		GINFO("Incoming data is 2D - starting 2D FOCUSS reconstruction\n");
	} else if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) == 1 && vDims.at(3) > 1) {
		//opCS_ = new CS_FOCUSS_2Dt();
		GINFO("Incoming data is 2Dt - starting 2Dt FOCUSS reconstruction\n");
	} else if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) > 1 && vDims.at(3) == 1) {
		// squeeze array due to x,y,z,c dimension of 3D FOCUSS class
		sum_dim(*tmp_m2->getObjectPtr(), 3, *tmp_m2->getObjectPtr());
		
		opCS_ = new CS_FOCUSS_3D();
		GINFO("Incoming data is 3D - starting 3D FOCUSS reconstruction\n");
	} else if (vDims.at(0) > 1 && vDims.at(1) > 1 && vDims.at(2) > 1 && vDims.at(3) > 1) {
		opCS_ = new CS_FOCUSS_4D();

		GINFO("Incoming data is 4D - starting 4D FOCUSS reconstruction\n");
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

	// process data in class member function
	opCS_->process(m1, tmp_m2);

	//Now pass on image
	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}

	return GADGET_OK;
}
GADGET_FACTORY_DECLARE(CS_LAB)
