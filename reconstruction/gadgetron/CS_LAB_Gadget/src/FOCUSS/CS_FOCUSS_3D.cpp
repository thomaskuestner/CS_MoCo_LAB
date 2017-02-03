/*
file name	: 	CS_FOCUSS_3D.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.2

date		: 	18.12.2015

description	: 	implementation of the class "CS_FOCUSS_3D" (file CS_FOCUSS.h)

reference	:
				- Gorodnitsky I. and Rao B.: "Sparse Signal Reconstruction from Limited Data Using FOCUSS: A Re-weigted Minimum Norm Algorithm", IEEE Transaction on Signal Processing 1997:45(3)
				- K�stner T. et al.: "ESPReSSo: A Compressed Sensing Partial k-Space Acquisition and Reconstruction", Proc. ISMRM 2014
*/

#include "CS_FOCUSS.h"

using namespace Gadgetron;
/*
// read config file
int CS_FOCUSS_3D::process_config(ACE_Message_Block* mb){

	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		bXMLControl_ = bXMLControl.value();
	#else
		bXMLControl_ = this->get_bool_value("XMLControl");
	#endif	

	if (bXMLControl_) {

		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1		
			iNOuter_ = OuterIterations.value();
		  	//iCGResidual_ = iCGResidual.value();
			iNInner_ = InnerIterations.value();
			//bESPRActiveCS_ = bESPRActiveCS.value();
			//cfLambda_ = lambda.value();
		#else
		  	// how to calculate the beta value
			iCGResidual_ = this->get_int_value("CG Beta");	

			// maximum number of FOCUSS iterations	
		  	iNOuter_ = this->get_int_value("OuterIterations");
		  	if (iNOuter_ <= 0) iNOuter_ = 2;	

			// maximum number of CG iterations	
		  	iNInner_ = this->get_int_value("InnerIterations");
		  	if (iNInner_ <= 0) iNInner_ = 20;
		
			// use ESPReSSo-constraint for pure CS data
			bESPRActiveCS_ = this->get_bool_value("CS - ESPReSSo");

			// FOCUSS lambda
			cfLambda_ = this->get_double_value("lambda");
		#endif
	}

	// p-value for the lp-norm
	fP_ = .5;

	// convergence boundary
	fEpsilon_ = (float)1e-6;

	// setup of the transformation parameters - sparsity dim, fft dim, ..
	fSetupTransformation();

	return GADGET_OK;
};*/

//--------------------------------------------------------------------------
//------------- process - CG-FOCUSS with additional constraints ------------
//--------------------------------------------------------------------------
int CS_FOCUSS_3D::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
	//------------------------------------------------------------------------
	//------------------------------ initial ---------------------------------
	//------------------------------------------------------------------------
	//--- set variables - store incoming data - permute incoming data --------
	//------------------------------------------------------------------------
	if (bDebug_)
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Starting FOCUSS reconstruction\n");
		#else
			GADGET_DEBUG1("Starting FOCUSS reconstruction\n");
		#endif

	// init member values based on header information
	fInitVal(m1);

	// declare recon output
	hoNDArray<std::complex<float> >  hacfOutput;

	// FOCUSS reconstruction - this function is also called by the Matlab implementation
	fRecon(*m2->getObjectPtr(), hacfOutput);

	// new GadgetContainer
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* cm2 = new GadgetContainerMessage<hoNDArray< std::complex<float> > >();

	// concatenate data with header
	m1->cont(cm2);

	// create output
	try{cm2->getObjectPtr()->create(&vtDim_);}
	catch (std::runtime_error &err){
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GEXCEPTION(err,"Unable to allocate new image array\n");
		#else
			GADGET_DEBUG_EXCEPTION(err,"Unable to allocate new image array\n");
		#endif
		m1->release();
		return -1;
	}

	// copy data
	memcpy(cm2->getObjectPtr()->get_data_ptr(), hacfOutput.begin(), sizeof(std::complex<float>)*hacfOutput.get_number_of_elements());

	// pass data on stream or to control object
	if (!bControl_){
		//Now pass on image
		if (this->next()->putq(m1) < 0) {
			return GADGET_FAIL;
		}
	}
	return GADGET_OK;
}

//--------------------------------------------------------------------------
//------------- 3D FOCUSS reconstruction - accessible from MATLAB ----------
//--------------------------------------------------------------------------
int CS_FOCUSS_3D::fRecon(hoNDArray<std::complex<float> >  &hacfInput, hoNDArray<std::complex<float> >  &hacfRecon){

	// input dimensions
	vtDim_ = *hacfInput.get_dimensions();

	// number of channels
	iNChannels_ = (int)vtDim_[3];

	// if Matlab is used, initialize singleton variables
	/*if (bMatlab_){
		for (int iI = 0; iI < 20; iI++){
			GlobalVar_FOCUSS::instance()->vbStatPrinc_.push_back(false);
			GlobalVar_FOCUSS::instance()->vfPrincipleComponents_.push_back(new hoNDArray<std::complex<float> > ());
		}
	}*/

	// const complex values
	const std::complex<float> cfZero(0.0);
	const std::complex<float> cfOne(1.0);

	// store incoming data in array
	hoNDArray<std::complex<float> >  hacfKSpace = hacfInput;

	// permute kSpace: x-y-z-c -> y-z-x-c
	std::vector<size_t> vtDimOrder; vtDimOrder.push_back(1); vtDimOrder.push_back(2); vtDimOrder.push_back(0); vtDimOrder.push_back(3);
	hacfKSpace = *permute(&hacfKSpace, &vtDimOrder,false);

	// update dim_ vector
	vtDim_.clear(); vtDim_ = *hacfKSpace.get_dimensions();

	//------------------------------------------------------------------------
	//-------------------------- sampling mask -------------------------------
	//------------------------------------------------------------------------
	hoNDArray<std::complex<float> >  hacfFullMask(hacfKSpace.get_dimensions()); hacfFullMask.fill(cfZero);
	hoNDArray<bool> habFullMask(hacfKSpace.get_dimensions()); habFullMask.fill(false);
	pcfPtr_ = hacfKSpace.get_data_ptr();
	pcfPtr2_ = hacfFullMask.get_data_ptr();
	pbPtr_ = habFullMask.get_data_ptr();
	for (int i = 0; i < hacfKSpace.get_number_of_elements(); i++){
		if (pcfPtr_[i] != cfZero){
			pcfPtr2_[i] = cfOne;
			pbPtr_[i] = true;
		}
	}

	//-------------------------------------------------------------------------
	//---------------- iFFT x direction - x ky kz ^= v (n�) -------------------
	//-------------------------------------------------------------------------
	if (Transform_fftBA_->get_active()){
		if (!bMatlab_ && bDebug_)
			#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				GDEBUG("FFT in read direction..\n");
			#else
				GADGET_DEBUG1("FFT in read direction..\n");
			#endif
		else if(bMatlab_ && bDebug_){
			// mexPrintf("FFT in read direction..\n");mexEvalString("drawnow;");
		}
		Transform_fftBA_->FTransform(hacfKSpace);
	}
	hoNDArray<std::complex<float> >  hacfWWindowed = hacfKSpace;

	//------------------------------------------------------------------------
	//---------------------------- windowing ---------------------------------
	//------------------------------------------------------------------------
	fGetCalibrationSize(habFullMask);
	fWindowing(hacfWWindowed);

	//------------------------------------------------------------------------
	//-------------------------- ESPReSSo init -------------------------------
	//------------------------------------------------------------------------
	//------ find symmetrical, conjugate sampled part - build filter ---------
	//------------------------------------------------------------------------
	if (bESPRActiveCS_ || fPartialFourierVal_ < 1.0)
		fInitESPReSSo(habFullMask);

	//-------------------------------------------------------------------------
	//----------------------- initial estimate --------------------------------
	//-------------------------------------------------------------------------
	if (!bMatlab_ && bDebug_)
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("Prepare initial estimate..\n");
		#else
			GADGET_DEBUG1("Prepare initial estimate..\n");
		#endif
	else if(bMatlab_ && bDebug_){
		// mexPrintf("Prepare initial estimate..\n"); mexEvalString("drawnow;");
	}

	// W in x-y-z-cha --> new base
	Transform_KernelTransform_->FTransform(hacfWWindowed);

	// W = abs(FTrafo(W)).^p - p = .5
	fAbsPow(hacfWWindowed, fP_);

	// calculate energy and divide windowed image for initial estimate
	hoNDArray<std::complex<float> >  hacfTotEnergy(hacfWWindowed.get_dimensions());
	pcfPtr_ = hacfTotEnergy.get_data_ptr();
	pcfPtr2_ = hacfWWindowed.get_data_ptr();
	for (int iCha = 0; iCha < iNChannels_; iCha++){
		size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iCha;
		hoNDArray<std::complex<float> >  hacfEnergyPerChannel(vtDim_[0], vtDim_[1], vtDim_[2], hacfWWindowed.get_data_ptr()+ tOffset, false);
		float fTmp = fCalcEnergy(hacfEnergyPerChannel);
		if (!bMatlab_ && bDebug_)
			#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				GDEBUG("energy in channel[%i]: %e..\n",iCha, fTmp);
			#else
				GADGET_DEBUG2("energy in channel[%i]: %e..\n",iCha, fTmp);
			#endif
		else if(bMatlab_ && bDebug_){
			// mexPrintf("energy in channel[%i]: %e..\n",iCha, fTmp); mexEvalString("drawnow;");
		}
		// fill channel
		#pragma  omp parallel for
		for (long i = 0; i < vtDim_[0]*vtDim_[1]*vtDim_[2]; i++){
			pcfPtr_[i+tOffset] = std::complex<float>(fTmp);
			pcfPtr2_[i+tOffset] /= fTmp;
		}
	}

	//-------------------------------------------------------------------------
	//------------------- iterative calculation -------------------------------
	//-------------------------------------------------------------------------
	// initial estimates for CG - all zero (except g_old)
	hoNDArray<std::complex<float> >  hacfQ(hacfWWindowed.get_dimensions());
	hoNDArray<std::complex<float> >  hacfRho = hacfQ, hacfG_old = hacfQ, hacfD = hacfQ, hacfRho_fft = hacfQ, hacfE = hacfQ, hacfG = hacfQ, hacfE_ifft = hacfQ, hacfBeta = hacfQ, hacfZ = hacfQ, hacfAlpha = hacfQ, hacfGradient_ESPReSSo = hacfQ;

	// outer loop - FOCUSS iterations
	for (int iOuter = 0; iOuter < iNOuter_; iOuter++){
		if (!bMatlab_ && bDebug_)
			#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				GDEBUG("FOCUSS loop: %i\n", iOuter);
			#else
				GADGET_DEBUG2("FOCUSS loop: %i\n", iOuter);
			#endif
		else if(bMatlab_ && bDebug_){
			// mexPrintf("FOCUSS loop: %i\n", iOuter);mexEvalString("drawnow;");
		}

		// reset initial values
		hacfRho.fill(cfZero); hacfD.fill(cfZero); hacfQ.fill(cfZero); hacfG_old.fill(std::complex<float>(1.0,1.0));

		// inner loop for CG
		for (int iInner = 0; iInner < iNInner_; iInner++){
			try{
				if (!bMatlab_ && bDebug_)
					#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
						GDEBUG("CG Loop: %i\n", iInner);
					#else
						GADGET_DEBUG2("CG Loop: %i\n", iInner);
					#endif
				else if(bMatlab_ && bDebug_){
					// mexPrintf("CG Loop: %i\n", iInner);	mexEvalString("drawnow;");
				}

				// rho: x-y-z ---> x-ky-kz
				hacfRho_fft = hacfRho;
				Transform_KernelTransform_->BTransform(hacfRho_fft);

				// e = v - Phi*F*rho - e: x-ky-kz
				fAminusBmultC(hacfKSpace,hacfFullMask,hacfRho_fft,hacfE);

				//l2 norm calculation - check epsilon
				std::vector<float> vfVec;
				for (int iCha = 0; iCha < iNChannels_; iCha++){
					size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iCha;
					hoNDArray<std::complex<float> >  eCha(vtDim_[0], vtDim_[1], vtDim_[2], hacfE.get_data_ptr()+ tOffset, false);
					vfVec.push_back(abs(dot(&eCha, &eCha)));
					vfVec[iCha] = std::sqrt(vfVec[iCha]);
					if (!bMatlab_ && bDebug_)
						#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
							GDEBUG("||e|| ch. %i  =  %e\n", iCha, vfVec[iCha]);
						#else
							GADGET_DEBUG2("||e|| ch. %i  =  %e\n", iCha, vfVec[iCha]);
						#endif
					else if(bMatlab_ && bDebug_){
						// mexPrintf("||e|| ch. %i  =  %e\n", iCha, vfVec[iCha]);mexEvalString("drawnow;");
					}
				}

				// how many channels are converged
				int iNom = 0;
				for (int iI = 0; iI < vfVec.size(); iI++){
					if (vfVec[iI] > fEpsilon_) iNom++;
				}
				if (!bMatlab_ && bDebug_)
					#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
						GDEBUG("number of non converged channels - %i\n", iNom);
					#else
						GADGET_DEBUG2("number of non converged channels - %i\n", iNom);
					#endif
				else if(bMatlab_ && bDebug_){
					// mexPrintf("number of non converged channels - %i\n", iNom);
					// mexEvalString("drawnow;");
				}

				// if all channels converged -> stop calculation
				if (iNom == 0) break;

				// e: x-ky-kz --> x-y-z
				hacfE_ifft = hacfE;
				Transform_KernelTransform_->FTransform(hacfE_ifft);
			}
			catch(...){
				if (!bMatlab_ && bDebug_)
					#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
						GDEBUG("Exception in first part..\n");
					#else
						GADGET_DEBUG1("Exception in first part..\n");
					#endif
				else if(bMatlab_ && bDebug_){
					// mexPrintf("Exception in first part..\n");mexEvalString("drawnow;");
				}
			}

			//------------------------------------------------------------------------
			//---------------------------- constraints -------------------------------
			//------------------------------------------------------------------------
			// emphasize conjugate similarity
			try{
				if (cfLambdaESPReSSo_ != cfZero && (bESPRActiveCS_ || fPartialFourierVal_ < 1.0)){
					hacfGradient_ESPReSSo = hacfRho;
					fGradESPReSSo(hacfGradient_ESPReSSo, hacfFullMask, hacfKSpace, hacfWWindowed, hacfQ);
				}
				else
					hacfGradient_ESPReSSo.fill(cfZero);
			}
			catch(...){
				if (!bMatlab_ && bDebug_)
					#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
						GDEBUG("Exception in ESPReSSo constraint\n");
					#else
						GADGET_DEBUG1("Exception in ESPReSSo constraint\n");
					#endif
				else  if(bMatlab_ && bDebug_){
					// mexPrintf("Exception in ESPReSSo constraint\n");
					// mexEvalString("drawnow;");
				}
				hacfGradient_ESPReSSo.fill(cfZero);
			}

			//-------------calculate gradient -------------------------
			// G = -conj(W).*IFFT(e)+Lambda.*Q + LambdaESPReSSo.*gradESPReSSo
			try{
				fCalcGradient(hacfWWindowed, hacfE_ifft, cfLambda_, hacfQ, cfLambdaESPReSSo_, hacfGradient_ESPReSSo, hacfG);
			}
			catch(...){
				if (!bMatlab_ && bDebug_)
					#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
						GDEBUG("Exception in gradient calculation\n");
					#else
						GADGET_DEBUG1("Exception in gradient calculation\n");
					#endif
				else if(bMatlab_ && bDebug_){
					// mexPrintf("Exception in gradient calculation\n");
					// mexEvalString("drawnow;");
				}
				hacfG.fill(cfZero);
			}

			//------------------- cg beta - Polak-Ribiere -------------------------------
			std::complex<float> fBetaCha (0.0);
			pcfPtr_ = hacfBeta.get_data_ptr();

			// loop over channels
			for (int iCha = 0; iCha < iNChannels_; iCha++){

				// fill sub array with data from higher order data array
				size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iCha;
				hoNDArray<std::complex<float> >  hacfSubArrayG_old(vtDim_[0], vtDim_[1], vtDim_[2], hacfG_old.get_data_ptr()+ tOffset, false);
				hoNDArray<std::complex<float> >  hacfSubArrayG(vtDim_[0], vtDim_[1], vtDim_[2], hacfG.get_data_ptr()+ tOffset, false);
				std::complex<float> fNumerator(0.0), fDenominator(0.0), fRightTerm(0.0);

				// calculate nominator
				pcfPtr2_ = hacfSubArrayG.get_data_ptr();
				for (long iI = 0; iI < hacfSubArrayG.get_number_of_elements(); iI++){
					fNumerator +=  pcfPtr2_[iI]*pcfPtr2_[iI];
				}

				// calculate denominator
				pcfPtr2_ = hacfSubArrayG_old.get_data_ptr();
				for (long iI = 0; iI < hacfSubArrayG.get_number_of_elements(); iI++){
					fDenominator +=  pcfPtr2_[iI]*pcfPtr2_[iI];
				}
				if (abs(fDenominator) != 0) fBetaCha = fNumerator / fDenominator;

				// fill part of the 3D array
				#pragma  omp parallel for
				for (long lI = 0; lI < vtDim_[0]*vtDim_[1]*vtDim_[2]; lI++)
					pcfPtr_[lI+tOffset] = fBetaCha;
			}
			//--------------------------------------------------------------------------

			// d = beta.*d - G and g_old = G
			fAmultBminusC(hacfBeta, hacfD, hacfG, hacfD);
			hacfG_old = hacfG;

			// z = Phi.*FFT(W.*d) - x-ky-kz
			multiply(hacfWWindowed, hacfD, hacfZ);
			Transform_KernelTransform_->BTransform(hacfZ);
			fMultiply(hacfZ, hacfFullMask);

			//---------------------------- cg alpha -------------------------------------
			//alpha(:,:,:,c) = (z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:));
			pcfPtr_ = hacfAlpha.get_data_ptr();
			for (int iCha = 0; iCha < iNChannels_; iCha++){
				std::complex<float> fAlphaCha (0.0);
				// fill sub array with data from higher order data array
				size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iCha;
				hoNDArray<std::complex<float> >  hacfSubArrayE(vtDim_[0], vtDim_[1], vtDim_[2], hacfE.get_data_ptr()+ tOffset, false);
				hoNDArray<std::complex<float> >  hacfSubArrayZ(vtDim_[0], vtDim_[1], vtDim_[2], hacfZ.get_data_ptr()+ tOffset, false);
				std::complex<float> fNumerator(0.0), fDenominator(0.0);
				// calculate nominator
				pcfPtr2_ = hacfSubArrayE.get_data_ptr();
				pcfPtr3_ = hacfSubArrayZ.get_data_ptr();
				for (long iI = 0; iI < hacfSubArrayE.get_number_of_elements(); iI++){
					fNumerator +=  std::conj(pcfPtr3_[iI])*pcfPtr2_[iI];
				}

				// calculate denominator
				for (long iI = 0; iI < hacfSubArrayZ.get_number_of_elements(); iI++){
					fDenominator +=  std::conj(pcfPtr3_[iI])*pcfPtr3_[iI];
				}
				if (abs(fDenominator) != 0) fAlphaCha = fNumerator / fDenominator;

				// fill 3D alpha array
				#pragma  omp parallel for
				for (long lI = 0; lI < vtDim_[0]*vtDim_[1]*vtDim_[2]; lI++)
					pcfPtr_[lI+tOffset] = fAlphaCha;
			}
			//--------------------------------------------------------------------------

			// q = q + alpha.*d
			fAplusBmultC(hacfQ, hacfAlpha, hacfD, hacfQ);

			// rho = W.*q
			multiply(hacfWWindowed, hacfQ, hacfRho);
		}

		// W = abs(rho).^p - p = .5 normalization
		hacfWWindowed = hacfRho;
		fAbsPowDivide(hacfWWindowed, fP_, hacfTotEnergy);
	}

	// rho = W.*q
	multiply(hacfWWindowed, hacfQ, hacfRho);

	//rho = kernelBTrafo(rho) -> x-y-z cart
	Transform_KernelTransform_->KernelBTransform(hacfRho);

	// output k-space for further k-space filtering or image
	if (Transform_fftAA_->get_active()){
		Transform_fftAA_->BTransform(hacfRho);
	}

	// permute output - rho: y-z-x-c -> x-y-z-c
	vtDimOrder.clear(); vtDimOrder.push_back(2); vtDimOrder.push_back(0); vtDimOrder.push_back(1); vtDimOrder.push_back(3);
	hacfRho = *permute(&hacfRho, &vtDimOrder,false);
	vtDim_.clear(); vtDim_ = *hacfRho.get_dimensions();

	// set output and return
	hacfRecon = hacfRho;

	if (!bMatlab_ && bDebug_)
			#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				GDEBUG("FOCUSS done..\n");
			#else
				GADGET_DEBUG1("FOCUSS done..\n");
			#endif
		else if(bMatlab_ && bDebug_){
			// mexPrintf("FOCUSS done..\n");
			// mexEvalString("drawnow;");
		}

	return GADGET_OK;
}

//--------------------------------------------------------------------------
//------------------ gradient of the ESPReSSo constraint -------------------
//--------------------------------------------------------------------------
void CS_FOCUSS_3D::fGradESPReSSo(hoNDArray<std::complex<float> > & hacfRho, hoNDArray<std::complex<float> > &hacfFullMask, hoNDArray<std::complex<float> > &hacfKSpace, hoNDArray<std::complex<float> > &hacfW, hoNDArray<std::complex<float> > &hacfQ){

	//--------------------------------------------------------------------------
	//------------------ transform to Cartesian k-space image ------------------
	//--------------------------------------------------------------------------
	// new base --> Cartesian base (ky-kz-x-c)
	hoNDArray<std::complex<float> >  hacfRhoTmp = hacfRho;
	Transform_KernelTransform_->BTransform(hacfRhoTmp);

	if (fAllZero(hacfRhoTmp))
		hacfRhoTmp = hacfKSpace;

	// back to kSpace - (ky-kz-x-c) --> (ky-kz-kx-c)
	if (Transform_fftBA_->get_active()) Transform_fftBA_->BTransform(hacfRhoTmp);

	//--------------------------------------------------------------------------
	//--------------------------- espresso_reconGrad ---------------------------
	//--------------------------------------------------------------------------

	// if error catched - return zero for non-affecting the reconstruction result with false calculations
	try{

		// flip in espresso direction
		hoNDArray<std::complex<float> >  hacfRhoFlipped = hacfRhoTmp;
		pcfPtr_ = hacfRhoTmp.get_data_ptr();
		pcfPtr2_ = hacfRhoFlipped.get_data_ptr();
		if (iESPReSSoDirection_ == 1)
				#pragma omp parallel for
				for (int iCha = 0; iCha < vtDim_[3]; iCha++)
					for (int idX = 0; idX < vtDim_[2]; idX++)
						for (int idZ = 0; idZ < vtDim_[1]; idZ++)
							for (int idY = 0; idY < vtDim_[0]; idY++)
								pcfPtr2_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pcfPtr_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];

		// partition encoding direction
		else if(iESPReSSoDirection_ == 2)
			#pragma omp parallel for
			for (int iCha = 0; iCha < vtDim_[3]; iCha++)
				for (int idX = 0; idX < vtDim_[2]; idX++)
					for (int idY = 0; idY < vtDim_[0]; idY++)
						for (int idZ = 0; idZ < vtDim_[1]; idZ++)
							pcfPtr2_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pcfPtr_[idY + (vtDim_[1]-1 - idZ)*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];


		// filter kSpace for low-resolution and smooth transitions (k-space times filter - s(kx,ky,kz)*H(kx,ky,kz))
		multiply(hacfRhoFlipped, hacfFilter_, hacfRhoFlipped);
		hoNDFFT<float>::instance()->fftshift3D(hacfRhoFlipped);

		// iFFT for initial image (all three spatial dimensions) - (ky,kz,kx)-->(y,z,x)
		hoNDFFT<float>::instance()->ifft3(hacfRhoFlipped);

		// get phase image - phase = exp(2i * angle(rhoTmp))
		hoNDArray<std::complex<float> >  hacfPhase = hacfRhoFlipped;
		pcfPtr_ = hacfPhase.get_data_ptr();
		#pragma omp parallel for
		for (int iI = 0; iI < hacfPhase.get_number_of_elements(); iI++)
			pcfPtr_[iI] = std::polar(1.0, 2.0*std::arg(pcfPtr_[iI]));

		// kSpaceCombi
		hoNDArray<std::complex<float> >  hacfKSpaceCombi = hacfKSpace;

		//--------------------------------------------------------------------------
		//--------------------------- conjugate symmetry ---------------------------
		//--------------------------------------------------------------------------
		// get data pointers
		pcfPtr2_ = hacfKSpaceCombi.get_data_ptr();
		pbPtr2_ = habMaskConj_.get_data_ptr();
		// mapping conjugate points
		// ESPReSSo direction is phase encoding direction
		// ESPReSSo direction is phase encoding direction
		if (iESPReSSoDirection_ == 1){
			for (int iCha = 0; iCha < vtDim_[3]; iCha++)
				for (int idX = 0; idX < vtDim_[2]; idX++)
					for (int idZ = 0; idZ < vtDim_[1]; idZ++)
						for (int idY = 0; idY < vtDim_[0]; idY++){
							if (habMaskConj_(idY, idZ, idX, iCha)==true){
								pcfPtr2_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = std::conj(pcfPtr2_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]]);
							}
						}
		}
		// ESPReSSo direction is partition encoding direction
		else if(iESPReSSoDirection_ == 2){
			for (int iCha = 0; iCha < vtDim_[3]; iCha++)
				for (int idX = 0; idX < vtDim_[2]; idX++)
					for (int idY = 0; idY < vtDim_[0]; idY++)
						for (int idZ = 0; idZ < vtDim_[1]; idZ++){
							if (habMaskConj_(idY, idZ, idX, iCha)==true){
								pcfPtr2_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = std::conj(pcfPtr2_[idY + (vtDim_[1]-1 - idZ)*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]]);
							}
						}
		}

		//--------------------------------------------------------------------------
		//------------------------- get ESPReSSo gradient --------------------------
		//--------------------------------------------------------------------------
		Transform_KernelTransform_->FTransform(hacfKSpaceCombi);
		fESPReSSoOut(hacfW, hacfQ, hacfPhase, hacfKSpaceCombi, hacfRho);
	}
	catch(...){
		if (!bMatlab_ && bDebug_)
			#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				GDEBUG("Error occured in ESPReSSo gradient calculation...return 0\n");
			#else
				GADGET_DEBUG1("Error occured in ESPReSSo gradient calculation...return 0\n");
			#endif
		else if(bMatlab_ && bDebug_){
			// mexPrintf("Error occured in ESPReSSo gradient calculation...return 0\n");mexEvalString("drawnow;");
		}
		clear(hacfRho);
	}
	return;
}

//--------------------------------------------------------------------------
//------------------- initialize the ESPReSSo constraint -------------------
//--------------------------------------------------------------------------
void CS_FOCUSS_3D::fInitESPReSSo(hoNDArray<bool>& habFullMask){

	// temporal array
	hoNDArray<bool> habSamplPtrSym;

	// get center points
	int iCenterX = (int)(vtDim_[2]/2), iCenterY = (int)(vtDim_[0]/2), iCenterZ = (int)(vtDim_[1]/2);

	if (cfLambdaESPReSSo_ != std::complex<float>(0.0)){

		// ESPReSSo acquisition
		if (fPartialFourierVal_ < 1.0 && fPartialFourierVal_ > 0.0){

			//-------------------------------------------------------------------------
			//------------------------- Upper / Lower ?  ------------------------------
			//------- check if upper or lower region in the k-space is sampled --------
			//-------------------------------------------------------------------------
			int iNumFoundUpper = 0, iNumFoundLower = 0;

			// check only for one channel - sampling mask is same for all channels
			pbPtr_ = habFullMask.get_data_ptr();

			// ESPReSSo direction is phase encoding direction
			if (iESPReSSoDirection_ == 1){
				for (int idY = (int)(fPartialFourierVal_*vtDim_[0])-1; idY < vtDim_[0]; idY++){
					if (pbPtr_[idY + vtDim_[0]*iCenterZ + vtDim_[0]*vtDim_[1]*iCenterX] != false){
						iNumFoundUpper++;
					}
				}
				for (int idY = (int)(fPartialFourierVal_*vtDim_[0])-1; idY >= 0; idY--){
					if (pbPtr_[idY + vtDim_[0]*iCenterZ + vtDim_[0]*vtDim_[1]*iCenterX] != false){
						iNumFoundLower++;
					}
				}
				(iNumFoundLower > iNumFoundUpper) ? bESPReSSoIsLower_ = true : bESPReSSoIsLower_ = false;
			}

			// ESPReSSo direction is partition encoding direction
			else{
				for (int idZ = (int)(fPartialFourierVal_*vtDim_[1])-1; idZ < vtDim_[1]; idZ++){
					if (pbPtr_[iCenterY + idZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] != false){
						iNumFoundUpper++;
					}
				}
				for (int idZ = (int)(fPartialFourierVal_*vtDim_[1])-1; idZ >= 0; idZ--){
					if (pbPtr_[iCenterY + idZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] != false){
						iNumFoundLower++;
					}
				}
				(iNumFoundLower > iNumFoundUpper) ? bESPReSSoIsLower_ = true : bESPReSSoIsLower_ = false;
			}

			//--------------------------------------------------------------------------
			//------------------------- symmetrical sampling pattern -------------------
			//- sampling mask is mirrored at the center line depending on ESPReSSo dir.-
			//--------------------------------------------------------------------------
			if (!bMatlab_ && bDebug_)
				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("compute symmetrical sampling pattern\n");
				#else
					GADGET_DEBUG1("compute symmetrical sampling pattern\n");
				#endif
			else if(bMatlab_ && bDebug_){
				// mexPrintf("compute symmetrical sampling pattern\n");mexEvalString("drawnow;");
			}

			habSamplPtrSym = habFullMask;
			pbPtr_ = habSamplPtrSym.get_data_ptr();

			// lower half sampled
			if (bESPReSSoIsLower_){

				// ESPReSSo direction is phase encoding direction
				if (iESPReSSoDirection_ == 1)
					#pragma omp parallel for
					for (int iCha = 0; iCha < vtDim_[3]; iCha++)
						for (int idX = 0; idX < vtDim_[2]; idX++)
							for (int idZ = 0; idZ < vtDim_[1]; idZ++)
								for (int idY = vtDim_[0]-1; idY > vtDim_[0]/2; idY--)
									pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];

				// ESPReSSo direction is partition encoding direction
				else
					#pragma omp parallel for
					for (int iCha = 0; iCha < vtDim_[3]; iCha++)
						for (int idX = 0; idX < vtDim_[2]; idX++)
							for (int idY = 0; idY < vtDim_[0]; idY++)
								for (int idZ = vtDim_[1]-1; idZ > vtDim_[1]/2; idZ--)
									pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr_[idY + (vtDim_[1]-1 - idZ)*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
			}

			// upper half sampled
			else{
				// ESPReSSo direction is phase encoding direction
				if (iESPReSSoDirection_ == 1)
					#pragma omp parallel for
					for (int iCha = 0; iCha < vtDim_[3]; iCha++)
						for (int idX = 0; idX < vtDim_[2]; idX++)
							for (int idZ = 0; idZ < vtDim_[1]; idZ++)
								for (int idY = 0; idY < vtDim_[0]/2; idY++)
									pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];


				// ESPReSSo direction is partition encoding direction
				else
					#pragma omp parallel for
					for (int iCha = 0; iCha < vtDim_[3]; iCha++)
						for (int idX = 0; idX < vtDim_[2]; idX++)
							for (int idY = 0; idY < vtDim_[0]; idY++)
								for (int idZ = 0; idZ < vtDim_[1]/2; idZ++)
									pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr_[idY + (vtDim_[1]-1 - idZ)*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
			}

			//-------------------------------------------------------------------------
			//--------------------------- conjugate sampling pattern ------------------
			//---- conjugation - get only unsampled points for complex conjugation ----
			//-------------------------------------------------------------------------
			if (!bMatlab_ && bDebug_)
				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("compute conjugate sampling pattern\n");
				#else
					GADGET_DEBUG1("compute conjugate sampling pattern\n");
				#endif
			else if(bMatlab_ && bDebug_){
				// mexPrintf("compute conjugate sampling pattern\n");
				// mexEvalString("drawnow;");
			}

			habMaskConj_.create(habFullMask.get_dimensions());
			pbPtr_ = habMaskConj_.get_data_ptr();
			pbPtr2_ = habFullMask.get_data_ptr();
			pbPtr3_ = habSamplPtrSym.get_data_ptr();

			// get unsampled points - XOR gating
			#pragma omp parallel for
			for (int i = 0; i < habMaskConj_.get_number_of_elements(); i++)
				pbPtr_[i] = pbPtr2_[i]^pbPtr3_[i];
		}

		// ESPReSSo constraint for pure CS data without Partial Fourier technique
		else{

			//--------------------------------------------------------------------------
			//----------- ESPReSSo for Compressed Sensing w/o Partial Fourier ----------
			//--------------------------------------------------------------------------
			habMaskRight_ = habFullMask;
			pbPtr_ = habMaskRight_.get_data_ptr();

			// set upper half zero
			#pragma omp parallel for
			for (int iCha = 0; iCha < vtDim_[3]; iCha++)
				for (int idX = 0; idX < vtDim_[2]; idX++)
					for (int idY = 0; idY < vtDim_[0]/2; idY++)
						for (int idZ = 0; idZ < vtDim_[1]; idZ++)
							pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = false;

			// set lower half zero
			habMaskLeft_ = habFullMask;
			pbPtr_ = habMaskLeft_.get_data_ptr();
			#pragma omp parallel for
			for (int iCha = 0; iCha < vtDim_[3]; iCha++)
				for (int idX = 0; idX < vtDim_[2]; idX++)
					for (int idY = vtDim_[0]/2; idY < vtDim_[0]; idY++)
						for (int idZ = 0; idZ < vtDim_[1]; idZ++)
							pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = false;

			// logical OR with kSpaceCenter
			pbPtr_  = habMaskRight_.get_data_ptr();
			pbPtr2_ = habKSpaceCenter_.get_data_ptr();
			pbPtr3_ = habMaskLeft_.get_data_ptr();
			#pragma omp parallel for
			for (long i = 0; i < habMaskRight_.get_number_of_elements(); i++)
				if (pbPtr2_[i] == true){
					pbPtr_[i] = true;
					pbPtr3_[i] = true;
				}

			// flip both masks (lmaskRight and lmaskLeft) and write to habMaskConj_ and habMaskConj2_
			habMaskConj_.create(habMaskRight_.get_dimensions());
			habMaskConj2_.create(habMaskRight_.get_dimensions());
			bool* pcfPtr_habMaskConj_ = habMaskConj_.get_data_ptr();
			bool* pcfPtr_habMaskConj2_= habMaskConj2_.get_data_ptr();
			#pragma omp parallel for
			for (int iCha = 0; iCha < vtDim_[3]; iCha++)
				for (int idX = 0; idX < vtDim_[2]; idX++)
					for (int idZ = 0; idZ < vtDim_[1]; idZ++)
						for (int idY = 0; idY < vtDim_[0]; idY++){
							pcfPtr_habMaskConj_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] =  pbPtr_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
							pcfPtr_habMaskConj2_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr3_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
						}

			pbPtr_ = habFullMask.get_data_ptr();
			pbPtr2_ = habMaskConj_.get_data_ptr();
			pbPtr3_ = habMaskConj2_.get_data_ptr();

			// get unsampled points - XOR gating
			#pragma omp parallel for
			for (int i = 0; i < habMaskConj_.get_number_of_elements(); i++){
				pbPtr2_[i] = pbPtr_[i]^pbPtr2_[i] || pbPtr_[i]^pbPtr3_[i];
			}

			// set values for symmetrical kSpace part / filter
			iESPReSSoDirection_ = 1;
			bESPReSSoIsLower_ = true;
		}

		//-------------------------------------------------------------------------
		//-------------------- symmetrical kSpace part / filter --------------------
		//--------------------------------------------------------------------------
		// get indices of center kSpace lines - find symmetrical kSpace part and calc filter
		if (!bMatlab_ && bDebug_)
			#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				GDEBUG("find symmetrical kSpace part..\n");
			#else
				GADGET_DEBUG1("find symmetrical kSpace part..\n");
			#endif
		else if(bMatlab_ && bDebug_){
			// mexPrintf("find symmetrical kSpace part..\n");mexEvalString("drawnow;");
		}

		// vectors for storing the outer sampled lines
		std::vector<int> viKSpaceLines_1;
		std::vector<int> viKSpaceLines_2;

		// loop over the symmetrical pattern beginning at the center line to the first unsampled point
		// --> difference between center line and first zero is half the symmetrical sampled lines
		// store the upper and lower symmetrical sampled line in the kSpaceLine vector

		// ESPReSSo is active
		if (fPartialFourierVal_ != 1.0){
			pbPtr_ = habSamplPtrSym.get_data_ptr();

			// phase encoding direction
			if (iESPReSSoDirection_ == 1){
				for (int iZ = 0; iZ < vtDim_[1]; iZ++){

					// lower region is sampled
					if (bESPReSSoIsLower_ == true){
						if (fPartialFourierVal_ != .5){
							viKSpaceLines_1.push_back(fPartialFourierVal_*vtDim_[0]-iCenterY);
							viKSpaceLines_2.push_back(fPartialFourierVal_*vtDim_[0]);
						}
						else{
							for(int iY = 0; iY < iCenterY; iY++){
								if (pbPtr_[iY + iZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] == true){
									// number of line in k-space
									viKSpaceLines_1.push_back(iY);
									// number of "mirrored" line
									viKSpaceLines_2.push_back(vtDim_[0]-iY);
									break;
								}
								else if(iY==iCenterY-1){
									// number of line in k-space
									viKSpaceLines_1.push_back(iY-1);
									// number of "mirrored" line
									viKSpaceLines_2.push_back(vtDim_[0]-iY+1);
								}
							}
						}
					}
					// upper region is sampled
					else{
						if (fPartialFourierVal_ != .5){
							viKSpaceLines_1.push_back(fPartialFourierVal_*vtDim_[0]-iCenterY);
							viKSpaceLines_2.push_back(fPartialFourierVal_*vtDim_[0]);
						}
						else{
							for(int iY = vtDim_[0]; iY > iCenterY; iY--){
								if (pbPtr_[iY + iZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] == true){
									// number of line in k-space
									viKSpaceLines_1.push_back(vtDim_[0]-iY);
									// number of "mirrored" line
									viKSpaceLines_2.push_back(iY);
									break;
								}
								else if(iY==iCenterY+1){
									// number of line in k-space
									viKSpaceLines_1.push_back(iY-1);
									// number of "mirrored" line
									viKSpaceLines_2.push_back(vtDim_[0]-iY+1);
								}
							}
						}
					}
				}
			}
			// partition encoding direction
			else{
				for (int iY = 0; iY < vtDim_[0]; iY++){

					// lower region is sampled
					if (bESPReSSoIsLower_ == true){
						if (fPartialFourierVal_ != .5){
							viKSpaceLines_1.push_back(fPartialFourierVal_*vtDim_[1]-iCenterZ);
							viKSpaceLines_2.push_back(fPartialFourierVal_*vtDim_[1]);
						}
						else{
							for(int iZ = 0; iZ < iCenterZ; iZ++){
								if (pbPtr_[iY + iZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] == true){
									// number of line in k-space
									viKSpaceLines_1.push_back(iZ);
									// number of "mirrored" line
									viKSpaceLines_2.push_back(vtDim_[1]-iZ);
									break;
								}
								else if(iZ==iCenterZ-1){
									// number of line in k-space
									viKSpaceLines_1.push_back(iZ-1);
									// number of "mirrored" line
									viKSpaceLines_2.push_back(vtDim_[1]-iZ+1);
								}
							}
						}
					}
					// upper region is sampled
					else{
						if (fPartialFourierVal_ != .5){
							viKSpaceLines_1.push_back(fPartialFourierVal_*vtDim_[1]-iCenterZ);
							viKSpaceLines_2.push_back(fPartialFourierVal_*vtDim_[1]);
						}
						else{
							for(int iZ = vtDim_[1]; iZ > iCenterZ; iZ--){
								if (pbPtr_[iY + iZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] == true){
									// number of line in k-space
									viKSpaceLines_1.push_back(vtDim_[1]-iZ);
									// number of "mirrored" line
									viKSpaceLines_2.push_back(iZ);
									break;
								}
								else if(iZ==iCenterZ+1){
									// number of line in k-space
									viKSpaceLines_1.push_back(iZ-1);
									// number of "mirrored" line
									viKSpaceLines_2.push_back(vtDim_[1]-iZ+1);
								}
							}
						}
					}
				}
			}
		}

		// additional complex conjugate symmetry constraint for normal Compressed Sensing data
		else{
			int iLine1 = .75*vtDim_[0]-iCenterY;
			int iLine2 = .75*vtDim_[0];
			for (int idZ = 0; idZ < vtDim_[1]; idZ++){
				viKSpaceLines_1.push_back(iLine1);
				viKSpaceLines_2.push_back(iLine2);
			}
		}

		// if no line is found --> filter width will be 0 --> take default values from window calculation
		if (viKSpaceLines_1.size() == 0 || viKSpaceLines_2.size()==0){
			viKSpaceLines_1.clear(); viKSpaceLines_2.clear();

			if (!bMatlab_ && bDebug_)
				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("Estimate filter size by calculated window");
				#else
					GADGET_DEBUG1("Estimate filter size by calculated window");
				#endif
			else if(bMatlab_ && bDebug_){
				// mexPrintf("Estimate filter size by calculated window");	mexEvalString("drawnow;");
			}

			// phase encoding direction
			if (iESPReSSoDirection_ == 1)
				for (int iZ = 0; iZ < vtDim_[1]; iZ++){
					int iLine1 = .75*vtDim_[0]-iCenterY;
					int iLine2 = .75*vtDim_[0];
					viKSpaceLines_1.push_back(iLine1);
					viKSpaceLines_2.push_back(iLine2);
				}

			// partition encoding direction
			else
				for (int iY = 0; iY < vtDim_[0]; iY++){
					int iLine1 = .75*vtDim_[1]-iCenterZ;
					int iLine2 = .75*vtDim_[1];
					viKSpaceLines_1.push_back(iLine1);
					viKSpaceLines_2.push_back(iLine2);
				}
		}
		if (!bMatlab_ && bDebug_)
			#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				GDEBUG("calculate filter..\n");
			#else
				GADGET_DEBUG1("calculate filter..\n");
			#endif
		else if(bMatlab_ && bDebug_){
			// mexPrintf("calculate filter..\n");mexEvalString("drawnow;");
		}

		// create filter array (inital all ones)
		hacfFilter_.create(habFullMask.get_dimensions());
		hacfFilter_.fill(std::complex<float>(1.0));
		hoNDArray<std::complex<float> >  hacfFilt1D(vtDim_[0], vtDim_[1]);
		hoNDArray<std::complex<float> >  hacfFilt1DFullArray(habFullMask.get_dimensions());

		// loop over the 3 spatial dimensions
		for (int iDim = 0; iDim < 3; iDim++){
			hacfFilt1D.fill(std::complex<float>(0.0)); hacfFilt1DFullArray.fill(std::complex<float>(0.0));

			// filter for ESPReSSo direction (Hanning)
			if (iDim == iESPReSSoDirection_){
				// phase encoding direction
				if (iESPReSSoDirection_ == 1){
					for (int iZ = 0; iZ < vtDim_[1]; iZ++){

						// get filter coefficients (passed parameter is window width)
						std::vector<float> HanningCoeff = fGetHanningWindow(viKSpaceLines_2.at(iZ)-viKSpaceLines_1.at(iZ));

						// place coefficients in the filt1D array - loop over filter "lines"
						for (int iI = 0, iY = viKSpaceLines_1.at(iZ); iY < viKSpaceLines_2.at(iZ); iY++, iI++)
							hacfFilt1D(iY, iZ) = std::complex<float>(HanningCoeff.at(iI));
					}
				}
				// partition encoding direction
				else if (iESPReSSoDirection_ == 2)
					for (int iY = 0; iY < vtDim_[0]; iY++){
						std::vector<float> vfHanningCoeff = fGetHanningWindow(viKSpaceLines_2.at(iY)-viKSpaceLines_1.at(iY));

						for (int iI = 0, iZ = viKSpaceLines_1.at(iY); iZ < viKSpaceLines_2.at(iY); iZ++, iI++)
							hacfFilt1D(iY, iZ) = std::complex<float>(vfHanningCoeff.at(iI));
					}

				// repetition in x-direction for getting a full 3D array
				for (int iCha = 0; iCha < vtDim_[3]; iCha++)
					for (int iX = 0; iX < vtDim_[2]; iX++)
						for (int iY = 0; iY < vtDim_[0]; iY++)
							for (int iZ = 0; iZ < vtDim_[1]; iZ++)
								hacfFilt1DFullArray(iY, iZ, iX, iCha) = hacfFilt1D(iY, iZ);

			}

			// filter for other direction (Hamming) - same procedure as in ESPReSSo direction
			else{
				// x-direction
				if (iDim == 0){
					std::vector<float> vfHammingCoeff = fGetHammingWindow(vtDim_[2]);
					for (int iCha = 0; iCha < vtDim_[3]; iCha++)
						for (int idX = 0; idX < vtDim_[2]; idX++)
							for (int idY = 0; idY < vtDim_[0]; idY++)
								for (int idZ = 0; idZ < vtDim_[1]; idZ++)
									hacfFilt1DFullArray(idY, idZ, idX, iCha) = std::complex<float>(vfHammingCoeff.at(idX));
				}

				// y-direction
				if (iDim == 1){
					std::vector<float> vfHammingCoeff = fGetHammingWindow(vtDim_[0]);
					for (int iCha = 0; iCha < vtDim_[3]; iCha++)
						for (int idX = 0; idX < vtDim_[2]; idX++)
							for (int idY = 0; idY < vtDim_[0]; idY++)
								for (int idZ = 0; idZ < vtDim_[1]; idZ++)
									hacfFilt1DFullArray(idY, idZ, idX, iCha) = std::complex<float>(vfHammingCoeff.at(idY));
				}

				// z-direction
				if (iDim == 2){
					std::vector<float> vfHammingCoeff = fGetHammingWindow(vtDim_[1]);
					for (int iCha = 0; iCha < vtDim_[3]; iCha++)
						for (int idX = 0; idX < vtDim_[2]; idX++)
							for (int idY = 0; idY < vtDim_[0]; idY++)
								for (int idZ = 0; idZ < vtDim_[1]; idZ++)
									hacfFilt1DFullArray(idY, idZ, idX, iCha) = std::complex<float>(vfHammingCoeff.at(idZ));
				}
			}
			hacfFilter_ *= hacfFilt1DFullArray;
		}

		// check maximum value in filter array
		pcfPtr_ = hacfFilter_.get_data_ptr();
		for (long iI = 0; iI < hacfFilter_.get_number_of_elements(); iI++){
			if (abs(pcfPtr_[iI]) > 1.0){
				pcfPtr_[iI] = std::complex<float>(0.0);
			}
		}
	}
}

//--------------------------------------------------------------------------
//---------------------------- windowing -----------------------------------
//--------------------------------------------------------------------------
void CS_FOCUSS_3D::fWindowing(hoNDArray< std::complex< float > > & hacfWWindowed){

	// array with mask
	hoNDArray< std::complex< float > >  hacfMask3D(hacfWWindowed.get_dimensions()); hacfMask3D.fill(std::complex< float >(0.0));

	// get calibration mask
	std::vector< size_t > vStart, vSize;
	for (int iI = 0; iI < 2; iI++){
		if (viCalibrationSize_.at(iI) % 2){
			vStart.push_back(std::floor((float)vtDim_[iI]/2)+std::ceil(-(float)viCalibrationSize_.at(iI)/2));
		}
		else{
			vStart.push_back(std::floor((float)vtDim_[iI]/2)+std::ceil(-(float)viCalibrationSize_.at(iI)/2));
		}
	}
	vStart.push_back(0); vStart.push_back(0);
	for (int iY = vStart.at(0); iY < vStart.at(0)+viCalibrationSize_.at(0); iY++)
		for (int iZ = vStart.at(1); iZ < vStart.at(1)+viCalibrationSize_.at(1); iZ++)
			for (int iX = vStart.at(2); iX < vStart.at(2)+viCalibrationSize_.at(2); iX++){
				for (int iC = vStart.at(3); iC < vStart.at(3)+viCalibrationSize_.at(3); iC++)
					hacfMask3D(iY, iZ, iX, iC) = std::complex< float >(1.0);
			}

	// windowing W
	hacfWWindowed *= hacfMask3D;

	// get kSpaceCenter
	habKSpaceCenter_.create(hacfWWindowed.get_dimensions());
	habKSpaceCenter_.fill(true);
	pbPtr_ = habKSpaceCenter_.get_data_ptr();
	pcfPtr_ = hacfMask3D.get_data_ptr();
	#pragma omp parallel for
	for (long lI = 0; lI < hacfMask3D.get_number_of_elements(); lI++)
		if(pcfPtr_[lI] == std::complex< float >(0.0))
			pbPtr_[lI] = false;

	if (!bMatlab_ && bDebug_)
		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GDEBUG("data windowed for initial estimate and kSpaceCenter found..\n");
		#else
			GADGET_DEBUG1("data windowed for initial estimate and kSpaceCenter found..\n");
		#endif
	else if(bMatlab_ && bDebug_){
		// mexPrintf("data windowed for initial estimate and kSpaceCenter found..\n");mexEvalString("drawnow;");
	}
}

//--------------------------------------------------------------------------
//---------------------- get calibration size ------------------------------
//--------------------------------------------------------------------------
void CS_FOCUSS_3D::fGetCalibrationSize(const hoNDArray<bool> &habArray){
	int iSY = 2, iSZ = 2;
	bool bYflag = false, bZflag = false;
	std::vector<size_t> vtDim = *habArray.get_dimensions();

	hoNDArray<bool> habMask_helper;

	while(!(bYflag && bZflag)){
		if (!bYflag){
			fCropArrYZ(habArray, iSY+1, iSZ, habMask_helper);
			if(fAllOne(habMask_helper))
				iSY++;
			else
				bYflag = true;
		}

		if (!bZflag){
			fCropArrYZ(habArray, iSY, iSZ+1, habMask_helper);
			if(fAllOne(habMask_helper))
				iSZ++;
			else
				bZflag = true;
		}

		if (iSY == vtDim[0])
			bYflag = true;
		if (iSZ == vtDim[1])
			bZflag = true;
	}

	// push values on calibration size vector
	viCalibrationSize_.push_back(iSY);
	viCalibrationSize_.push_back(iSZ);
	viCalibrationSize_.push_back(vtDim[2]);
	viCalibrationSize_.push_back(vtDim[3]);

	for (std::vector<int>::size_type i = 0; i != viCalibrationSize_.size(); i++){
		if (!bMatlab_ && bDebug_)
			#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				GDEBUG("calibration size: %i..\n", viCalibrationSize_.at(i));
			#else
				GADGET_DEBUG2("calibration size: %i..\n", viCalibrationSize_.at(i));
			#endif
		else if(bMatlab_ && bDebug_){
			// mexPrintf("calibration size: %i..\n", viCalibrationSize_.at(i));mexEvalString("drawnow;");
		}
	}
}

GADGET_FACTORY_DECLARE(CS_FOCUSS_3D)
