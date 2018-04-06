///*
//file name	: 	CS_FOCUSS_2Dt.cpp
//
//author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
//
//version		: 	1.0
//
//date		: 	13.02.2015
//
//description	: 	implementation of the class "CS_FOCUSS_2Dt" (file CS_FOCUSS.h)
//
//reference	:	MATLAB implementation of K�stner, T.
//*/
//
//#include "CS_FOCUSS.h"
//
//using namespace Gadgetron;
////--------------------------------------------------------------------------
////------------- process - CG-FOCUSS with additional constraints ------------
////--------------------------------------------------------------------------
//int CS_FOCUSS_2Dt::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
//{
//	//------------------------------------------------------------------------
//	//------------------------------ initial ---------------------------------
//	//------------------------------------------------------------------------
//	//--- set variables - store incoming data - permute incoming data --------
//	//------------------------------------------------------------------------
//	GDEBUG("Starting FOCUSS reconstruction\n");
//
//	// init member values based on header information
//	fInitVal(m1);
//
//	// dimensions of incoming data - only handle 3D
//	vtDim_ = *m2->getObjectPtr()->get_dimensions();
//	iNChannels_ = (int)vtDim_[3];
//
//	hoNDArray<std::complex<float> >  hacfOutput;
//	fRecon(*m2->getObjectPtr(), hacfOutput);
//
//	// put calculated image on stream
//	*m2->getObjectPtr() = hacfOutput;
//
//	// pass data on stream or to control object
//	if (!bControl_){
//		//Now pass on image
//		if (this->next()->putq(m1) < 0) {
//			return GADGET_FAIL;
//		}
//	}
//	return GADGET_OK;
//}
//
//int CS_FOCUSS_2Dt::fRecon(hoNDArray<std::complex<float> >  &hacfInput, hoNDArray<std::complex<float> >  &hacfRecon){
//	// input dimensions
//	vtDim_ = *hacfInput.get_dimensions();
//
//	// number of channels
//	iNChannels_ = (int)vtDim_[3];
//
//	// if Matlab is used, initialize singleton variables
//	if (bMatlab_){
//		for (int iI = 0; iI < 20; iI++){
//			GlobalVar_FOCUSS::instance()->vbStatPrinc_.push_back(false);
//			GlobalVar_FOCUSS::instance()->vfPrincipleComponents_.push_back(new hoNDArray<std::complex<float> > ());
//		}
//	}
//
//	// const complex values
//	const std::complex<float> cfZero(0.0);
//	const std::complex<float> cfOne(1.0);
//
//	// store incoming data in array
//	hoNDArray<std::complex<float> >  hacfKSpace = hacfInput;
//
//	//------------------------------------------------------------------------
//	//-------------------------- sampling mask -------------------------------
//	//------------------------------------------------------------------------
//	hoNDArray<std::complex<float> >  hacfFullMask(hacfKSpace.get_dimensions()); hacfFullMask.fill(cfZero);
//	hoNDArray<bool> habFullMask(hacfKSpace.get_dimensions()); habFullMask.fill(false);
//	pcfPtr_ = hacfKSpace.get_data_ptr();
//	pcfPtr2_ = hacfFullMask.get_data_ptr();
//	pbPtr_ = habFullMask.get_data_ptr();
//	for (int i = 0; i < hacfKSpace.get_number_of_elements(); i++)
//		if (pcfPtr_[i] != cfZero){
//			pcfPtr2_[i] = cfOne;
//			pbPtr_[i] = true;
//		}
//
//	//-------------------------------------------------------------------------
//	//---------------- iFFT x direction - x ky kz ^= v (n�) -------------------
//	//-------------------------------------------------------------------------
//	if (Transform_fftBA_->get_active()){
// 		GDEBUG("FFT in read direction..\n");
//		Transform_fftBA_->FTransform(hacfKSpace);
//	}
//	hoNDArray<std::complex<float> >  hacfWWindowed = hacfKSpace;
//
//	//------------------------------------------------------------------------
//	//---------------------------- windowing ---------------------------------
//	//------------------------------------------------------------------------
//	fGetCalibrationSize(habFullMask);
//	fWindowing(hacfWWindowed);
//
//	/*-------------------------------------------------------------------------
//	------------------------- initial estimate --------------------------------
//	--------------------------------------------------------------------------*/
// 	GDEBUG("Prepare initial estimate..\n");
//
//	// W in x-y-z-cha --> new base
//	Transform_KernelTransform_->FTransform(hacfWWindowed);
//
//	// W = abs(FTrafo(W)).^p - p = .5
//	fAbsPow(hacfWWindowed, fP_);
//
//	// calculate energy and divide windowed image for initial estimate
//	hoNDArray<std::complex<float> >  hacfTotEnergy(hacfWWindowed.get_dimensions());
//	pcfPtr_ = hacfTotEnergy.get_data_ptr();
//	pcfPtr2_ = hacfWWindowed.get_data_ptr();
//	for (int iCha = 0; iCha < iNChannels_; iCha++){
//		size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iCha;
//		hoNDArray<std::complex<float> >  hacfEnergyPerChannel(vtDim_[0], vtDim_[1], vtDim_[2], hacfWWindowed.get_data_ptr()+ tOffset, false);
//		float fTmp = fCalcEnergy(hacfEnergyPerChannel);
// 		GDEBUG("energy in channel[%i]: %e..\n",iCha, fTmp);
//		// fill channel
//		#pragma  omp parallel for
//		for (long i = 0; i < vtDim_[0]*vtDim_[1]*vtDim_[2]; i++){
//			pcfPtr_[i+tOffset] = std::complex<float>(fTmp);
//			pcfPtr2_[i+tOffset] /= fTmp;
//		}
//	}
//	GDEBUG("Prepare initial estimate..\n");
//
//	/*-------------------------------------------------------------------------
//	--------------------- iterative calculation -------------------------------
//	--------------------------------------------------------------------------*/
//	// initial estimates for CG - all zero (except g_old)
//	hoNDArray<std::complex<float> >  hacfQ(hacfWWindowed.get_dimensions());
//	hoNDArray<std::complex<float> >  hacfRho = hacfQ, hacfG_old = hacfQ, hacfD = hacfQ, hacfRho_fft = hacfQ, hacfE = hacfQ, hacfG = hacfQ, hacfE_ifft = hacfQ, hacfBeta = hacfQ, hacfZ = hacfQ, hacfAlpha = hacfQ, hacfGradient_ESPReSSo = hacfQ;
//
//	// outer loop for FOCUSS
//	for (int iOuter = 0; iOuter < iNOuter_; iOuter++){
// 		GDEBUG("FOCUSS loop: %i\n", iOuter);
//
//		// reset initial values
//		hacfRho.fill(cfZero); hacfD.fill(cfZero); hacfQ.fill(cfZero); hacfG_old.fill(std::complex<float>(1.0,1.0));
//
//		// inner loop for CG
//		for (int iInner = 0; iInner < iNInner_; iInner++){
//			try{
// 				GDEBUG("CG Loop: %i\n", iInner);
//
//				// rho: x-y-z ---> x-ky-kz
//				hacfRho_fft = hacfRho;
//				Transform_KernelTransform_->BTransform(hacfRho_fft);
//
//				// e = v - Phi*F*rho - e: x-ky-kz
//				fAminusBmultC(hacfKSpace,hacfFullMask,hacfRho_fft,hacfE);
//
//				//l2 norm calculation - check epsilon
//				std::vector<float> vfVec;
//				for (int iCha = 0; iCha < iNChannels_; iCha++){
//					size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iCha;
//					hoNDArray<std::complex<float> >  eCha(vtDim_[0], vtDim_[1], vtDim_[2], hacfE.get_data_ptr()+ tOffset, false);
//					vfVec.push_back(abs(dot(&eCha, &eCha)));
//					vfVec[iCha] = std::sqrt(vfVec[iCha]);
// 					GDEBUG("||e|| ch. %i  =  %e\n", iCha, vfVec[iCha]);
//				}
//
//				// how many channels are converged
//				int iNom = 0;
//				for (int iI = 0; iI < vfVec.size(); iI++){
//					if (vfVec[iI] > fEpsilon_) iNom++;
//				}
// 				GDEBUG("number of non converged channels - %i\n", iNom);
//
//				// if all channels converged -> stop calculation
//				if (iNom == 0) break;
//
//				// e: x-ky-kz --> x-y-z
//				hacfE_ifft = hacfE;
//				Transform_KernelTransform_->FTransform(hacfE_ifft);
//
//				//------------------------------------------------------------------------
//				//---------------------------- constraints -------------------------------
//				//------------------------------------------------------------------------
//
//				// emphasize conjugate similarity - not used in 2Dt
//				hoNDArray<std::complex<float> >  hacfGradient_ESPReSSo = hacfRho; hacfGradient_ESPReSSo.fill(0.0);
//
//				//----------------- gradient -------------------------
//				// G = -conj(W).*IFFT(e)+Lambda.*Q
//				fCalcGradient(hacfWWindowed, hacfE_ifft, cfLambda_, hacfQ, cfLambdaESPReSSo_, hacfGradient_ESPReSSo, hacfG);
//
//				//------------------- cg beta - Polak-Ribiere -------------------------------
//				std::complex<float> fBetaCha (0.0);
//				pcfPtr_ = hacfBeta.get_data_ptr();
//
//				// loop over channels
//				for (int iCha = 0; iCha < iNChannels_; iCha++){
//
//					// fill sub array with data from higher order data array
//					size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iCha;
//					hoNDArray<std::complex<float> >  hacfSubArrayG_old(vtDim_[0], vtDim_[1], vtDim_[2], hacfG_old.get_data_ptr()+ tOffset, false);
//					hoNDArray<std::complex<float> >  hacfSubArrayG(vtDim_[0], vtDim_[1], vtDim_[2], hacfG.get_data_ptr()+ tOffset, false);
//					std::complex<float> fNumerator(0.0), fDenominator(0.0), fRightTerm(0.0);
//
//					// calculate nominator
//					pcfPtr2_ = hacfSubArrayG.get_data_ptr();
//					for (long iI = 0; iI < hacfSubArrayG.get_number_of_elements(); iI++){
//						fNumerator +=  pcfPtr2_[iI]*pcfPtr2_[iI];
//					}
//
//					// calculate denominator
//					pcfPtr2_ = hacfSubArrayG_old.get_data_ptr();
//					for (long iI = 0; iI < hacfSubArrayG.get_number_of_elements(); iI++){
//						fDenominator +=  pcfPtr2_[iI]*pcfPtr2_[iI];
//					}
//					if (abs(fDenominator) != 0) fBetaCha = fNumerator / fDenominator;
//
//					// fill part of the 3D array
//					#pragma  omp parallel for
//					for (long lI = 0; lI < vtDim_[0]*vtDim_[1]*vtDim_[2]; lI++)
//						pcfPtr_[lI+tOffset] = fBetaCha;
//				}
//				//--------------------------------------------------------------------------
//
//				// d = beta.*d - G and g_old = G
//				fAmultBminusC(hacfBeta, hacfD, hacfG, hacfD);
//				hacfG_old = hacfG;
//
//				// z = Phi.*FFT(W.*d) - x-ky-kz
//				multiply(hacfWWindowed, hacfD, hacfZ);
//				Transform_KernelTransform_->BTransform(hacfZ);
//				fMultiply(hacfZ, hacfFullMask);
//
//				//---------------------------- cg alpha -------------------------------------
//				//alpha(:,:,:,c) = (z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:));
//				pcfPtr_ = hacfAlpha.get_data_ptr();
//				for (int iCha = 0; iCha < iNChannels_; iCha++){
//					std::complex<float> fAlphaCha (0.0);
//					// fill sub array with data from higher order data array
//					size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iCha;
//					hoNDArray<std::complex<float> >  hacfSubArrayE(vtDim_[0], vtDim_[1], vtDim_[2], hacfE.get_data_ptr()+ tOffset, false);
//					hoNDArray<std::complex<float> >  hacfSubArrayZ(vtDim_[0], vtDim_[1], vtDim_[2], hacfZ.get_data_ptr()+ tOffset, false);
//					std::complex<float> fNumerator(0.0), fDenominator(0.0);
//					// calculate nominator
//					pcfPtr2_ = hacfSubArrayE.get_data_ptr();
//					pcfPtr3_ = hacfSubArrayZ.get_data_ptr();
//					for (long iI = 0; iI < hacfSubArrayE.get_number_of_elements(); iI++){
//						fNumerator +=  std::conj(pcfPtr3_[iI])*pcfPtr2_[iI];
//					}
//
//					// calculate denominator
//					for (long iI = 0; iI < hacfSubArrayZ.get_number_of_elements(); iI++){
//						fDenominator +=  std::conj(pcfPtr3_[iI])*pcfPtr3_[iI];
//					}
//					if (abs(fDenominator) != 0) fAlphaCha = fNumerator / fDenominator;
//
//					// fill 3D alpha array
//					#pragma  omp parallel for
//					for (long lI = 0; lI < vtDim_[0]*vtDim_[1]*vtDim_[2]; lI++)
//						pcfPtr_[lI+tOffset] = fAlphaCha;
//				}
//				//--------------------------------------------------------------------------
//
//				// q = q + alpha.*d
//				fAplusBmultC(hacfQ, hacfAlpha, hacfD, hacfQ);
//
//				// rho = W.*q
//				multiply(hacfWWindowed, hacfQ, hacfRho);
//			}
//			catch(...){
//				GDEBUG("exception in inner loop..\n");
//			}
//		}
//
//		// W = abs(rho).^p - p = .5 normalization
//		hacfWWindowed = hacfRho;
//		fAbsPowDivide(hacfWWindowed, fP_, hacfTotEnergy);
//	}
//
//	// rho = W.*q
//	multiply(hacfWWindowed, hacfQ, hacfRho);
//
//	//rho = kernelBTrafo(rho) -> x-y-z cart
//	Transform_KernelTransform_->KernelBTransform(hacfRho);
//
//	// output k-space for further k-space filtering or image
//	if (Transform_fftAA_->get_active()){
//		Transform_fftAA_->BTransform(hacfRho);
//	}
//
//	// set output and return
//	hacfRecon = hacfRho;
//
// 	GDEBUG("FOCUSS done..\n");
//
//	return GADGET_OK;
//}
//
////--------------------------------------------------------------------------
////---------------------------- windowing -----------------------------------
////--------------------------------------------------------------------------
//void CS_FOCUSS_2Dt::fWindowing(hoNDArray<std::complex<float> > & hacfWWindowed){
//	// array with mask
//	hoNDArray<std::complex<float> >  hacfMask3D(hacfWWindowed.get_dimensions()); hacfMask3D.fill(std::complex<float>(0.0));
//
//	// get calibration mask
//	std::vector<size_t> vStart, vSize;
//	for (int iI = 0; iI < 2; iI++){
//		if (viCalibrationSize_.at(iI) % 2){
//			vStart.push_back(std::floor((float)vtDim_[iI]/2)+std::ceil(-(float)viCalibrationSize_.at(iI)/2));
//		}
//		else{
//			vStart.push_back(std::floor((float)vtDim_[iI]/2)+std::ceil(-(float)viCalibrationSize_.at(iI)/2));
//		}
//	}
//	vStart.push_back(0); vStart.push_back(0);
//	for (int iY = vStart.at(0); iY < vStart.at(0)+viCalibrationSize_.at(0); iY++)
//		for (int iZ = vStart.at(1); iZ < vStart.at(1)+viCalibrationSize_.at(1); iZ++)
//			for (int iX = vStart.at(2); iX < vStart.at(2)+viCalibrationSize_.at(2); iX++){
//				for (int iC = vStart.at(3); iC < vStart.at(3)+viCalibrationSize_.at(3); iC++)
//					hacfMask3D(iY, iZ, iX, iC) = std::complex<float>(1.0);
//			}
//
//	// windowing W
//	hacfWWindowed *= hacfMask3D;
//
//	// get kSpaceCenter
//	habKSpaceCenter_.create(hacfWWindowed.get_dimensions());
//	habKSpaceCenter_.fill(true);
//	pbPtr_ = habKSpaceCenter_.get_data_ptr();
//	pcfPtr_ = hacfMask3D.get_data_ptr();
//	#pragma omp parallel for
//	for (long lI = 0; lI < hacfMask3D.get_number_of_elements(); lI++)
//		if(pcfPtr_[lI] == std::complex<float>(0.0))
//			pbPtr_[lI] = false;
//
// 	GDEBUG("data windowed for initial estimate and kSpaceCenter found..\n");
//}
//
//// get calibration size
//void CS_FOCUSS_2Dt::fGetCalibrationSize(const hoNDArray<bool> &habArray){
//	int iSY = 2;
//	bool bYflag = false;
//	std::vector<size_t> vtDim = *habArray.get_dimensions();
//	bool * pbArray = habArray.get_data_ptr();
//
//	while(!(bYflag)){
//		if (!bYflag){
//			for (int iY = std::ceil((float)vtDim[0]/2)-iSY+1; iY < std::ceil((float)vtDim[0]/2)+iSY+1; iY++){
//				if (!pbArray)
//					bYflag = true;
//				else
//					iSY++;
//			}
//		}
//		if (iSY == vtDim[0])
//			bYflag = true;
//	}
//
//	// push values on calibration size vector
//	viCalibrationSize_.push_back(iSY);
//	viCalibrationSize_.push_back(vtDim[1]);
//	viCalibrationSize_.push_back(vtDim[2]);
//	viCalibrationSize_.push_back(vtDim[3]);
//
//	for (std::vector<int>::size_type i = 0; i != viCalibrationSize_.size(); i++){
// 		GDEBUG("calibration size: %i..\n", viCalibrationSize_.at(i));
//	}
//}
//
//GADGET_FACTORY_DECLARE(CS_FOCUSS_2Dt)
