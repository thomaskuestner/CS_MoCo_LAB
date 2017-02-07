/*
file name	: 	CS_FOCUSS_2D.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.1

date		: 	13.02.2015

description	: 	implementation of the class "CS_FOCUSS_2D" (file CS_FOCUSS.h)

notes		:	no Total Variation and ESPReSSo constraint
*/

#include "CS_FOCUSS.h"

using namespace Gadgetron;

int CS_FOCUSS_2D::process_config(ACE_Message_Block* mb){

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
			cfLambda_ = lambda.value();			
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



//--------------------------------------------------------------------------
//------------- process - CG-FOCUSS with additional constraints ------------
//--------------------------------------------------------------------------
int CS_FOCUSS_2D::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
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


int CS_FOCUSS_2D::fRecon(hoNDArray<std::complex<float> >  &hacfInput, hoNDArray<std::complex<float> >  &hacfRecon){

	// input dimensions
	vtDim_ = *hacfInput.get_dimensions();

	// number of channels
	iNChannels_ = (int)vtDim_[2];

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

	// permute kSpace: x-y-c -> y-x-c
	std::vector<size_t> vtDimOrder; vtDimOrder.push_back(1); vtDimOrder.push_back(0); vtDimOrder.push_back(2);
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
	for (int i = 0; i < hacfKSpace.get_number_of_elements(); i++)
		if (pcfPtr_[i] != cfZero){
			pcfPtr2_[i] = cfOne;
			pbPtr_[i] = true;
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
	if (!bMatlab_ && bDebug_)
			#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
				GDEBUG("get calib size..\n");
			#else
				GADGET_DEBUG1("get calib size..\n");
			#endif
		else if(bMatlab_ && bDebug_){
			// mexPrintf("get calib size..\n");mexEvalString("drawnow;");
		}
	fGetCalibrationSize(habFullMask);
	fWindowing(hacfWWindowed);

	/*-------------------------------------------------------------------------
	------------------------- initial estimate --------------------------------
	--------------------------------------------------------------------------*/
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
		size_t tOffset = vtDim_[0]*vtDim_[1]*iCha;
		hoNDArray<std::complex<float> >  hacfEnergyPerChannel(vtDim_[0], vtDim_[1], hacfWWindowed.get_data_ptr()+ tOffset, false);
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
		for (long i = 0; i < vtDim_[0]*vtDim_[1]; i++){
			pcfPtr_[i+tOffset] = std::complex<float>(fTmp);
			pcfPtr2_[i+tOffset] /= fTmp;
		}
	}
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
		GDEBUG("Prepare initial estimate..\n");
	#else
		GADGET_DEBUG1("Prepare initial estimate..\n");
	#endif

	/*-------------------------------------------------------------------------
	--------------------- iterative calculation -------------------------------
	--------------------------------------------------------------------------*/
	// initial estimates for CG - all zero (except g_old)
	hoNDArray<std::complex<float> >  hacfQ(hacfWWindowed.get_dimensions());
	hoNDArray<std::complex<float> >  hacfRho = hacfQ, hacfG_old = hacfQ, hacfD = hacfQ, hacfRho_fft = hacfQ, hacfE = hacfQ, hacfG = hacfQ, hacfE_ifft = hacfQ, hacfBeta = hacfQ, hacfZ = hacfQ, hacfAlpha = hacfQ, hacfGradient_ESPReSSo = hacfQ;

	// outer loop for FOCUSS
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

				// rho: x-y ---> x-ky
				hacfRho_fft = hacfRho;
				Transform_KernelTransform_->BTransform(hacfRho_fft);

				// e = v - Phi*F*rho - e: x-ky
				fAminusBmultC(hacfKSpace,hacfFullMask,hacfRho_fft,hacfE);

				//l2 norm calculation - check epsilon
				std::vector<float> vfVec;
				for (int iCha = 0; iCha < iNChannels_; iCha++){
					size_t tOffset = vtDim_[0]*vtDim_[1]*iCha;
					hoNDArray<std::complex<float> >  eCha(vtDim_[0], vtDim_[1], hacfE.get_data_ptr()+ tOffset, false);
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

				// e: x-ky --> x-y
				hacfE_ifft = hacfE;
				Transform_KernelTransform_->FTransform(hacfE_ifft);

				//------------------------------------------------------------------------
				//---------------------------- constraints -------------------------------
				//------------------------------------------------------------------------

				// emphasize conjugate similarity - not used in 2D
				hoNDArray<std::complex<float> >  hacfGradient_ESPReSSo = hacfRho; hacfGradient_ESPReSSo.fill(0.0);

				//----------------- gradient -------------------------
				// G = -conj(W).*IFFT(e)+Lambda.*Q
				fCalcGradient(hacfWWindowed, hacfE_ifft, cfLambda_, hacfQ, cfLambdaESPReSSo_, hacfGradient_ESPReSSo, hacfG);

				//------------------- cg beta - Polak-Ribiere -------------------------------
				std::complex<float> fBetaCha (0.0);
				pcfPtr_ = hacfBeta.get_data_ptr();

				// loop over channels
				for (int iCha = 0; iCha < iNChannels_; iCha++){

					// fill sub array with data from higher order data array
					size_t tOffset = vtDim_[0]*vtDim_[1]*iCha;
					hoNDArray<std::complex<float> >  hacfSubArrayG_old(vtDim_[0], vtDim_[1], hacfG_old.get_data_ptr()+ tOffset, false);
					hoNDArray<std::complex<float> >  hacfSubArrayG(vtDim_[0], vtDim_[1], hacfG.get_data_ptr()+ tOffset, false);
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
					for (long lI = 0; lI < vtDim_[0]*vtDim_[1]; lI++)
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
					size_t tOffset = vtDim_[0]*vtDim_[1]*iCha;
					hoNDArray<std::complex<float> >  hacfSubArrayE(vtDim_[0], vtDim_[1], hacfE.get_data_ptr()+ tOffset, false);
					hoNDArray<std::complex<float> >  hacfSubArrayZ(vtDim_[0], vtDim_[1], hacfZ.get_data_ptr()+ tOffset, false);
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
					for (long lI = 0; lI < vtDim_[0]*vtDim_[1]; lI++)
						pcfPtr_[lI+tOffset] = fAlphaCha;
				}
				//--------------------------------------------------------------------------

				// q = q + alpha.*d
				fAplusBmultC(hacfQ, hacfAlpha, hacfD, hacfQ);

				// rho = W.*q
				multiply(hacfWWindowed, hacfQ, hacfRho);
			}
			catch(...){
				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
					GDEBUG("exception in inner loop..\n");
				#else
					GADGET_DEBUG1("exception in inner loop..\n");
				#endif
			}
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

	// permute output - rho: y-x-c -> x-y-c
	vtDimOrder.clear(); vtDimOrder.push_back(1); vtDimOrder.push_back(0); vtDimOrder.push_back(2);
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
//---------------------------- windowing -----------------------------------
//--------------------------------------------------------------------------
void CS_FOCUSS_2D::fWindowing(hoNDArray<std::complex<float> > & hacfWWindowed){
	// array with mask
	hoNDArray<std::complex<float> >  hacfMask2D(hacfWWindowed.get_dimensions()); hacfMask2D.fill(std::complex<float>(0.0));

	// get calibration mask
	std::vector<size_t> vStart, vSize;
	for (int iI = 0; iI < 1; iI++){
		if (viCalibrationSize_.at(iI) % 2){
			vStart.push_back(std::floor((float)vtDim_[iI]/2)+std::ceil(-(float)viCalibrationSize_.at(iI)/2));
		}
		else{
			vStart.push_back(std::floor((float)vtDim_[iI]/2)+std::ceil(-(float)viCalibrationSize_.at(iI)/2));
		}
	}
	vStart.push_back(0); vStart.push_back(0);
	for (int iY = vStart.at(0); iY < vStart.at(0)+viCalibrationSize_.at(0); iY++)
		for (int iX = vStart.at(1); iX < vStart.at(1)+viCalibrationSize_.at(1); iX++){
			for (int iC = vStart.at(2); iC < vStart.at(2)+viCalibrationSize_.at(2); iC++)
				hacfMask2D(iY, iX, iC) = std::complex<float>(1.0);
		}

	// windowing W
	hacfWWindowed *= hacfMask2D;

	// get kSpaceCenter
	habKSpaceCenter_.create(hacfWWindowed.get_dimensions());
	habKSpaceCenter_.fill(true);
	pbPtr_ = habKSpaceCenter_.get_data_ptr();
	pcfPtr_ = hacfMask2D.get_data_ptr();
	#pragma omp parallel for
	for (long lI = 0; lI < hacfMask2D.get_number_of_elements(); lI++)
		if(pcfPtr_[lI] == std::complex<float>(0.0))
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
void CS_FOCUSS_2D::fGetCalibrationSize(const hoNDArray<bool> &habArray){
	int iSY = 2;
	bool bYflag = false;
	std::vector<size_t> vtDim = *habArray.get_dimensions();
	bool * pbArray = habArray.get_data_ptr();

	while(!(bYflag)){
		if (!bYflag){
			for (int iY = std::ceil((float)vtDim[0]/2)-iSY+1; iY < std::ceil((float)vtDim[0]/2)+iSY+1; iY++){
				if (!habArray(iY, std::ceil((float)vtDim[1]/2)))//!pbArray)
					bYflag = true;
				else
					iSY++;
			}
		}
		if (iSY == vtDim[0])
			bYflag = true;
	}

	// push values on calibration size vector
	viCalibrationSize_.push_back(iSY);
	viCalibrationSize_.push_back(vtDim[1]);
	viCalibrationSize_.push_back(vtDim[2]);

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

GADGET_FACTORY_DECLARE(CS_FOCUSS_2D)
