#include "CS_FOCUSS.h"

using namespace Gadgetron;

//--------------------------------------------------------------------------
//------------- process - CG-FOCUSS with additional constraints ------------
//--------------------------------------------------------------------------
int CS_FOCUSS_4D::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	//------------------------------------------------------------------------
	//------------------------------ initial ---------------------------------
	//------------------------------------------------------------------------
	//--- set variables - store incoming data - permute incoming data --------
	//------------------------------------------------------------------------
	GINFO("Starting FOCUSS reconstruction\n");

	// init member values based on header information
	fInitVal(m1);

	// declare recon output
	hoNDArray<std::complex<float> > hacfOutput;

	// FOCUSS reconstruction - this function is also called by the Matlab implementation
	fRecon(*m2->getObjectPtr(), hacfOutput);

	// new GadgetContainer
	GadgetContainerMessage<hoNDArray<std::complex<float> > > *cm2 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();

	// concatenate data with header
	m1->cont(cm2);

	// create output
	try {
		cm2->getObjectPtr()->create(&vtDim_);
	} catch (std::runtime_error &err) {
		GEXCEPTION(err,"Unable to allocate new image array\n");
		m1->release();

		return GADGET_FAIL;
	}

	// copy data
	memcpy(cm2->getObjectPtr()->get_data_ptr(), hacfOutput.begin(), sizeof(std::complex<float>)*hacfOutput.get_number_of_elements());

	// pass data on stream or to control object
	if (!bControl_) {
		//Now pass on image
		if (this->next()->putq(m1) < 0) {
			return GADGET_FAIL;
		}
	}

	return GADGET_OK;
}

//--------------------------------------------------------------------------
//------------- 4D FOCUSS reconstruction - accessible from MATLAB ----------
//--------------------------------------------------------------------------
int CS_FOCUSS_4D::fRecon(hoNDArray<std::complex<float> > &hacfInput, hoNDArray<std::complex<float> > &hacfRecon)
{
	iNorm_ = 1;

	// input dimensions (x-y-z-t-c)
	vtDim_ = *hacfInput.get_dimensions();

	// number of channels
	iNChannels_ = (int)vtDim_[4];

	// const complex values
	const std::complex<float> cfZero(0.0);
	const std::complex<float> cfOne(1.0);

	// store incoming data in array
	hoNDArray<std::complex<float> > hacfKSpace = hacfInput;

	// permute kSpace: kx-ky-kz-t-c -> t-ky-kz-kx-c
	std::vector<size_t> vtDimOrder;
	vtDimOrder.push_back(3);
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(0);
	vtDimOrder.push_back(4);
	hacfKSpace = *permute(&hacfKSpace, &vtDimOrder, false);

	// update dim_ vector
	vtDim_.clear();
	vtDim_ = *hacfKSpace.get_dimensions();

	//------------------------------------------------------------------------
	//-------------------------- sampling mask -------------------------------
	//------------------------------------------------------------------------
	hoNDArray<std::complex<float> > hacfFullMask(hacfKSpace.get_dimensions());
	hoNDArray<bool> habFullMask(hacfKSpace.get_dimensions());

	for (size_t i = 0; i < hacfKSpace.get_number_of_elements(); i++) {
		if (hacfKSpace.at(i) != cfZero) {
			hacfFullMask.at(i) = cfOne;
			habFullMask.at(i) = true;
		} else {
			hacfFullMask.at(i) = cfZero;
			habFullMask.at(i) = false;
		}
	}

	//-------------------------------------------------------------------------
	//---------------- iFFT x direction - x ky kz ^= v (nü) -------------------
	//-------------------------------------------------------------------------
	if (Transform_fftBA_->get_active()) {
		GINFO("FFT in read direction..\n");
		Transform_fftBA_->FTransform(hacfKSpace);
	}

	hoNDArray<std::complex<float> > hacfWWindowed = hacfKSpace;

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
	if (GlobalVar::instance()->bESPRActiveCS_ || fPartialFourierVal_ < 1.0) {
		fInitESPReSSo(habFullMask);
	}

	//-------------------------------------------------------------------------
	//----------------------- initial estimate --------------------------------
	//-------------------------------------------------------------------------
	GINFO("Prepare initial estimate..\n");

	// W in x-y-z-cha --> new base
	Transform_KernelTransform_->FTransform(hacfWWindowed);

	// W = abs(FTrafo(W)).^p - p = .5
	fAbsPow(hacfWWindowed, fP_);

	// calculate energy and divide windowed image for initial estimate
	hoNDArray<std::complex<float> > hacfTotEnergy(hacfWWindowed.get_dimensions());
	for (int iCha = 0; iCha < iNChannels_; iCha++) {
		size_t tOffset = hacfWWindowed.get_size(0)*hacfWWindowed.get_size(1)*hacfWWindowed.get_size(2)*hacfWWindowed.get_size(3)*iCha;
		hoNDArray<std::complex<float> > hacfEnergyPerChannel(hacfWWindowed.get_size(0), hacfWWindowed.get_size(1), hacfWWindowed.get_size(2), hacfWWindowed.get_size(3), hacfWWindowed.get_data_ptr() + tOffset, false);
		float channel_max_energy;

		if (iNorm_ == 0) {
			channel_max_energy = fCalcEnergy(hacfEnergyPerChannel);
		} else if (iNorm_ == 1) {
			// hacfEnergyPerChannel only contains float (imag=0), so taking real part only is enough
			channel_max_energy = abs(hacfEnergyPerChannel.at(amax(hacfEnergyPerChannel)).real());
		} else {
			channel_max_energy = 1.0;
		}

		GDEBUG("energy in channel[%i]: %e..\n", iCha, channel_max_energy);

		// fill channel
		#pragma omp parallel for
		for (size_t i = 0; i < hacfEnergyPerChannel.get_number_of_elements(); i++) {
			hacfTotEnergy.at(i+tOffset) = std::complex<float>(channel_max_energy);
			hacfWWindowed.at(i+tOffset) /= std::complex<float>(channel_max_energy);
		}
	}

	//-------------------------------------------------------------------------
	//------------------- iterative calculation -------------------------------
	//-------------------------------------------------------------------------
	// initial estimates for CG - all zero (except g_old)
	hoNDArray<std::complex<float> > hacfQ(hacfWWindowed.get_dimensions());
	hacfQ.fill(cfZero);

	hoNDArray<std::complex<float> > hacfRho = hacfQ, hacfG_old = hacfQ, hacfD = hacfQ, hacfRho_fft = hacfQ, hacfE = hacfQ, hacfG = hacfQ, hacfE_ifft = hacfQ, hacfBeta = hacfQ, hacfZ = hacfQ, hacfAlpha = hacfQ, hacfGradient_ESPReSSo = hacfQ;

	// outer loop - FOCUSS iterations
	for (int iOuter = 0; iOuter < GlobalVar::instance()->iNOuter_; iOuter++) {
		GDEBUG("FOCUSS loop: %i\n", iOuter);

		// reset initial values
		hacfRho.fill(cfZero);
		hacfD.fill(cfZero);
		hacfQ.fill(cfZero);
		hacfG_old.fill(std::complex<float>(1.0,1.0));

		// inner loop for CG
		for (int iInner = 0; iInner < GlobalVar::instance()->iNInner_; iInner++) {
			try {
				GDEBUG("CG Loop: %i\n", iInner);

				// rho: x-y-z ---> x-ky-kz
				hacfRho_fft = hacfRho;
				Transform_KernelTransform_->BTransform(hacfRho_fft);

				// e = v - Phi*F*rho - e: x-ky-kz
				multiply(hacfFullMask, hacfRho_fft, hacfE);
				subtract(hacfKSpace, hacfE, hacfE);

				//l2 norm calculation - check epsilon
				std::vector<float> vfVec;
				for (int iCha = 0; iCha < iNChannels_; iCha++) {
					size_t tOffset = hacfE.get_size(0)*hacfE.get_size(1)*hacfE.get_size(2)*hacfE.get_size(3)*iCha;
					hoNDArray<std::complex<float> > eCha(hacfE.get_size(0), hacfE.get_size(1), hacfE.get_size(2), hacfE.get_size(3), hacfE.get_data_ptr() + tOffset, false);
					vfVec.push_back(abs(dot(&eCha, &eCha)));
					vfVec[iCha] = std::sqrt(vfVec[iCha]);

					GDEBUG("||e|| ch. %i  =  %e\n", iCha, vfVec[iCha]);
				}

				// how many channels are converged
				unsigned int iNom = 0;
				for (size_t iI = 0; iI < vfVec.size(); iI++) {
					if (vfVec[iI] > fEpsilon_) {
						iNom++;
					}
				}

				GDEBUG("number of non converged channels - %i\n", iNom);

				// if all channels converged -> stop calculation
				if (iNom == 0) {
					break;
				}
			} catch (...) {
				GERROR("Exception in first part..\n");
			}

			//------------------------------------------------------------------------
			//---------------------------- constraints -------------------------------
			//------------------------------------------------------------------------
			// emphasize conjugate similarity
			try {
				if (GlobalVar::instance()->cfLambdaESPReSSo_ != cfZero
					&& (GlobalVar::instance()->bESPRActiveCS_ || fPartialFourierVal_ < 1.0)
				) {
					hacfGradient_ESPReSSo = hacfRho;
					fGradESPReSSo(hacfGradient_ESPReSSo, hacfFullMask, hacfKSpace, hacfWWindowed, hacfQ);
				} else {
					hacfGradient_ESPReSSo.fill(cfZero);
				}
			} catch (...) {
				GERROR("Exception in ESPReSSo constraint\n");
				hacfGradient_ESPReSSo.fill(cfZero);
			}

			//-------------calculate gradient -------------------------
			// G = -conj(W).*IFFT(e)+Lambda.*Q + LambdaESPReSSo.*gradESPReSSo
			try {
				// e: x-ky-kz --> x-y-z
				hacfE_ifft = hacfE;
				Transform_KernelTransform_->FTransform(hacfE_ifft);

				fCalcGradient(hacfWWindowed, hacfE_ifft, GlobalVar::instance()->cfLambda_, hacfQ, GlobalVar::instance()->cfLambdaESPReSSo_, hacfGradient_ESPReSSo, hacfG);
			} catch (...) {
				GERROR("Exception in gradient calculation\n");
				hacfG.fill(cfZero);
			}

			//------------------- cg beta - Polak-Ribiere -------------------------------
			#pragma omp parallel for
			for (int iCha = 0; iCha < iNChannels_; iCha++) {
				// fill sub array with data from higher order data array
				size_t tOffset = hacfG.get_size(0)*hacfG.get_size(1)*hacfG.get_size(2)*hacfG.get_size(3)*iCha;
				hoNDArray<std::complex<float> > hacfSubArrayG_old(hacfG.get_size(0), hacfG.get_size(1), hacfG.get_size(2), hacfG.get_size(3), hacfG_old.get_data_ptr() + tOffset, false);
				hoNDArray<std::complex<float> > hacfSubArrayG(hacfG.get_size(0), hacfG.get_size(1), hacfG.get_size(2), hacfG.get_size(3), hacfG.get_data_ptr() + tOffset, false);
				float fBetaCha = 0.0;
				float fNumerator = 0.0;
				float fDenominator = 0.0;

				// calculate nominator
				for (size_t iI = 0; iI < hacfSubArrayG.get_number_of_elements(); iI++) {
					fNumerator += std::pow(hacfSubArrayG.at(iI).real(), 2) + std::pow(hacfSubArrayG.at(iI).imag(), 2);		// = std::conj(hacfSubArrayG.at(iI))*hacfSubArrayG.at(iI)
				}

				// calculate denominator
				for (size_t iI = 0; iI < hacfSubArrayG.get_number_of_elements(); iI++) {
					fDenominator += std::pow(hacfSubArrayG_old.at(iI).real(), 2) + std::pow(hacfSubArrayG_old.at(iI).imag(), 2);		// = std::conj(hacfSubArrayG_old.at(iI))*hacfSubArrayG_old.at(iI)
				}

				if (abs(fDenominator) != 0) {
					fBetaCha = fNumerator / fDenominator;
				}

				// fill sub-array of the 4D array
				#pragma omp parallel for
				for (size_t lI = 0; lI < hacfSubArrayG.get_number_of_elements(); lI++) {
					hacfBeta.at(lI+tOffset) = std::complex<float>(fBetaCha);
				}
			}
			//--------------------------------------------------------------------------

			// d = beta.*d - G and g_old = G
			multiply(hacfBeta, hacfD, hacfD);
			subtract(hacfD, hacfG, hacfD);
			hacfG_old = hacfG;

			// z = Phi.*FFT(W.*d) - x-ky-kz
			multiply(hacfWWindowed, hacfD, hacfZ);
			Transform_KernelTransform_->BTransform(hacfZ);
			multiply(hacfZ, hacfFullMask, hacfZ);

			//---------------------------- cg alpha -------------------------------------
			//alpha(:,:,:,c) = (z_helper(:)'*e_helper(:))/(z_helper(:)'*z_helper(:));
			#pragma omp parallel for
			for (int iCha = 0; iCha < iNChannels_; iCha++) {
				std::complex<float> fAlphaCha = 0.0;

				// fill sub array with data from higher order data array
				size_t tOffset = hacfE.get_size(0)*hacfE.get_size(1)*hacfE.get_size(2)*hacfE.get_size(3)*iCha;
				hoNDArray<std::complex<float> > hacfSubArrayE(hacfE.get_size(0), hacfE.get_size(1), hacfE.get_size(2), hacfE.get_size(3), hacfE.get_data_ptr() + tOffset, false);
				hoNDArray<std::complex<float> > hacfSubArrayZ(hacfE.get_size(0), hacfE.get_size(1), hacfE.get_size(2), hacfE.get_size(3), hacfZ.get_data_ptr() + tOffset, false);
				std::complex<float> fNumerator = 0.0;
				std::complex<float> fDenominator = 0.0;

				// calculate nominator
				for (size_t iI = 0; iI < hacfSubArrayE.get_number_of_elements(); iI++) {
					fNumerator += std::conj(hacfSubArrayZ.at(iI))*hacfSubArrayE.at(iI);
				}

				// calculate denominator
				for (size_t iI = 0; iI < hacfSubArrayZ.get_number_of_elements(); iI++) {
					fDenominator += std::conj(hacfSubArrayZ.at(iI))*hacfSubArrayZ.at(iI);
				}

				if (abs(fDenominator) != 0) {
					fAlphaCha = fNumerator / fDenominator;
				}

				// fill 3D alpha array
				#pragma omp parallel for
				for (size_t lI = 0; lI < hacfSubArrayE.get_number_of_elements(); lI++) {
					hacfAlpha.at(lI+tOffset) = fAlphaCha;
				}
			}
			//--------------------------------------------------------------------------

			// q = q + alpha.*d
			multiply(hacfAlpha, hacfD, hacfD);
			add(hacfQ, hacfD, hacfQ);

			// rho = W.*q
			multiply(hacfWWindowed, hacfQ, hacfRho);
		}

		// W = abs(rho).^p - p = .5 normalization
		hacfWWindowed = hacfRho;
		fAbsPow(hacfWWindowed, fP_);

		if (iNorm_ == 1) {
			#pragma omp parallel for
			for (int iCha = 0; iCha < iNChannels_; iCha++) {
				// re-calculate channel energy if iNorm == self scaling
				size_t tOffset = hacfWWindowed.get_size(0)*hacfWWindowed.get_size(1)*hacfWWindowed.get_size(2)*hacfWWindowed.get_size(3)*iCha;
				hoNDArray<std::complex<float> > hacfEnergyPerChannel(hacfWWindowed.get_size(0), hacfWWindowed.get_size(1), hacfWWindowed.get_size(2), hacfWWindowed.get_size(3), hacfWWindowed.get_data_ptr()+ tOffset, false);
				// hacfEnergyPerChannel only contains float (imag=0), so taking real part only is enough
				float channel_max_energy = abs(hacfEnergyPerChannel.at(amax(hacfEnergyPerChannel)).real());

				// don't divide by zero
				if (channel_max_energy == 0.0) {
					GWARN("channel_max_energy is 0!\n");
					continue;
				}

				// perform division at all elements in channel
				hacfEnergyPerChannel /= std::complex<float>(channel_max_energy);
			}
		} else if (iNorm_ == 0) {
			// perform division with total energies
			divide(hacfWWindowed, hacfTotEnergy, hacfWWindowed);
		} else {
			GERROR("Illegal state of iNorm_ (=%d)!", iNorm_);
			return -1;
		}
	}

	// rho = W.*q
	multiply(hacfWWindowed, hacfQ, hacfRho);

	//rho = kernelBTrafo(rho) -> x-y-z cart
	Transform_KernelTransform_->KernelBTransform(hacfRho);

	// output k-space for further k-space filtering or image
	if (Transform_fftAA_->get_active()) {
		Transform_fftAA_->BTransform(hacfRho);
	}

	// permute output - rho: t-y-z-x-c -> x-y-z-t-c
	vtDimOrder.clear();
	vtDimOrder.push_back(3);
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(0);
	vtDimOrder.push_back(4);
	hacfRho = *permute(&hacfRho, &vtDimOrder,false);

	vtDim_.clear();
	vtDim_ = *hacfRho.get_dimensions();

	// set output and return
	hacfRecon = hacfRho;

	GINFO("FOCUSS done..\n");

	return GADGET_OK;
}

//--------------------------------------------------------------------------
//------------------ gradient of the ESPReSSo constraint -------------------
//--------------------------------------------------------------------------
void CS_FOCUSS_4D::fGradESPReSSo(hoNDArray<std::complex<float> > &hacfRho, hoNDArray<std::complex<float> > &hacfFullMask, hoNDArray<std::complex<float> > &hacfKSpace, hoNDArray<std::complex<float> > &hacfW, hoNDArray<std::complex<float> > &hacfQ)
{
	//--------------------------------------------------------------------------
	//------------------ transform to Cartesian k-space image ------------------
	//--------------------------------------------------------------------------
	// new base --> Cartesian base
	hoNDArray<std::complex<float> > hacfRhoTmp = hacfRho;
	Transform_KernelTransform_->BTransform(hacfRhoTmp);

	if (fAllZero(hacfRhoTmp)) {
		hacfRhoTmp = hacfKSpace;
	}

	// back to kSpace
	if (Transform_fftBA_->get_active()) {
		Transform_fftBA_->BTransform(hacfRhoTmp);
	}

	//--------------------------------------------------------------------------
	//--------------------------- espresso_reconGrad ---------------------------
	//--------------------------------------------------------------------------
	// permute array: t-y-z-x-c -> y-z-x-c-t
	std::vector<size_t> vtDimOrder;
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(3);
	vtDimOrder.push_back(4);
	vtDimOrder.push_back(0);
	hacfRhoTmp = *permute(&hacfRhoTmp, &vtDimOrder,false);
	vtDim_ = *hacfRhoTmp.get_dimensions();
	hacfKSpace = *permute(&hacfKSpace, &vtDimOrder, false);

	hoNDArray< std::complex<float> > hacfKSpaceCombi(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], vtDim_[4]);
	hacfKSpaceCombi.fill( std::complex< float >(0.0, 0.0) );

	hoNDArray< std::complex<float> > hacfPhase(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], vtDim_[4]);
	hacfPhase.fill( std::complex< float >(0.0, 0.0) );

	// if error catched - return zero for non-affecting the reconstruction result with false calculations
	try {
		// loop over time states/gates
		for (size_t idT = 0; idT < vtDim_[0]; idT++) {
			// get 3D+c sub-array
			size_t lOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*vtDim_[3]*idT;
			hoNDArray<std::complex<float> > hacfSub3D(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], hacfRhoTmp.get_data_ptr()+lOffset, false);
			hoNDArray<std::complex<float> > hacfSubKSpace(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], hacfKSpace.get_data_ptr()+lOffset, false);
			hoNDArray<std::complex<float> > hacfFilter3D(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], hacfFilter_.get_data_ptr()+lOffset, false);
			hoNDArray<bool> habMaskConj3D(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], habMaskConj_.get_data_ptr()+lOffset, false);

			// flip in espresso direction
			hoNDArray<std::complex<float> > hacfRhoFlipped3D = hacfSub3D;
			pcfPtr_ = hacfSub3D.get_data_ptr();
			pcfPtr2_ = hacfRhoFlipped3D.get_data_ptr();
			if (iESPReSSoDirection_ == 1) {//ESPReSSo direction - y
					#pragma omp parallel for
					for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
						for (size_t idX = 0; idX < vtDim_[2]; idX++) {
							for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
								for (size_t idY = 0; idY < vtDim_[0]; idY++) {
									pcfPtr2_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pcfPtr_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
								}
							}
						}
					}
			} else if (iESPReSSoDirection_ == 2) {//ESPReSSo directon - z
				// partition encoding direction
				#pragma omp parallel for
				for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
					for (size_t idX = 0; idX < vtDim_[2]; idX++) {
						for (size_t idY = 0; idY < vtDim_[0]; idY++) {
							for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
								pcfPtr2_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pcfPtr_[idY + (vtDim_[1]-1 - idZ)*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
							}
						}
					}
				}
			}

			// filter kSpace for low-resolution and smooth transitions (k-space times filter - s(kx,ky,kz)*H(kx,ky,kz))
			multiply(hacfRhoFlipped3D, hacfFilter3D, hacfRhoFlipped3D);
			hoNDFFT_CS<float>::instance()->fftshift3D(hacfRhoFlipped3D);

			// iFFT for initial image (all three spatial dimensions) - (ky,kz,kx)-->(y,z,x)
			hoNDFFT_CS<float>::instance()->ifft3(hacfRhoFlipped3D);

			// get phase image - phase = exp(2i * angle(rhoTmp))
			hoNDArray<std::complex<float> > hacfPhase3D = hacfRhoFlipped3D;
			pcfPtr_ = hacfPhase3D.get_data_ptr();

			#pragma omp parallel for
			for (size_t iI = 0; iI < hacfPhase3D.get_number_of_elements(); iI++) {
				pcfPtr_[iI] = std::polar(1.0, 2.0*std::arg(pcfPtr_[iI]));
			}

			// kSpaceCombi
			hoNDArray<std::complex<float> > hacfKSpaceCombi3D = hacfSubKSpace;

			//--------------------------------------------------------------------------
			//--------------------------- conjugate symmetry ---------------------------
			//--------------------------------------------------------------------------
			// get data pointers
			pcfPtr2_ = hacfKSpaceCombi3D.get_data_ptr();
			pbPtr2_ = habMaskConj3D.get_data_ptr();

			// mapping conjugate points
			// ESPReSSo direction is phase encoding direction
			if (iESPReSSoDirection_ == 1) {
				for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
					for (size_t idX = 0; idX < vtDim_[2]; idX++) {
						for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
							for (size_t idY = 0; idY < vtDim_[0]; idY++) {
								if (habMaskConj_(idY, idZ, idX, iCha)==true) {
									pcfPtr2_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = std::conj(pcfPtr2_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]]);
								}
							}
						}
					}
				}
			} else if (iESPReSSoDirection_ == 2) {
				// ESPReSSo direction is partition encoding direction
				for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
					for (size_t idX = 0; idX < vtDim_[2]; idX++) {
						for (size_t idY = 0; idY < vtDim_[0]; idY++) {
							for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
								if (habMaskConj_(idY, idZ, idX, iCha)==true) {
									pcfPtr2_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = std::conj(pcfPtr2_[idY + (vtDim_[1]-1 - idZ)*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]]);
								}
							}
						}
					}
				}
			}

			memcpy(hacfKSpaceCombi.get_data_ptr() + lOffset, hacfKSpaceCombi3D.get_data_ptr(), sizeof(std::complex<float>)*hacfKSpaceCombi3D.get_number_of_elements());
			memcpy(hacfPhase.get_data_ptr() + lOffset, hacfPhase3D.get_data_ptr(), sizeof(std::complex<float>)*hacfPhase3D.get_number_of_elements());
		}

		// permute array: y-z-x-c-t -> t-y-z-x-c
		vtDimOrder.clear();
		vtDimOrder.push_back(4);
		vtDimOrder.push_back(0);
		vtDimOrder.push_back(1);
		vtDimOrder.push_back(2);
		vtDimOrder.push_back(3);
		hacfKSpaceCombi = *permute(&hacfKSpaceCombi, &vtDimOrder,false);
		vtDim_ = *hacfKSpaceCombi.get_dimensions();
		hacfPhase= *permute(&hacfPhase, &vtDimOrder, false);

		//--------------------------------------------------------------------------
		//------------------------- get ESPReSSo gradient --------------------------
		//--------------------------------------------------------------------------
		Transform_KernelTransform_->FTransform(hacfKSpaceCombi);
		fESPReSSoOut(hacfW, hacfQ, hacfPhase, hacfKSpaceCombi, hacfRho);
	} catch (...) {
		GERROR("Error occured in ESPReSSo gradient calculation.\n");
		hacfRho.clear();//clear(hacfRho);

		return;
	}

	return;
}

//--------------------------------------------------------------------------
//------------------- initialize the ESPReSSo constraint -------------------
//--------------------------------------------------------------------------
void CS_FOCUSS_4D::fInitESPReSSo(hoNDArray<bool>& habFullMask)
{
	// permute array: t-y-z-x-c -> y-z-x-c-t
	std::vector<size_t> vtDimOrder;
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(3);
	vtDimOrder.push_back(4);
	vtDimOrder.push_back(0);
	habFullMask = *permute(&habFullMask, &vtDimOrder,false);
	vtDim_ = *habFullMask.get_dimensions();

	// loop over states/gates
	for (size_t idT = 0; idT < vtDim_[4]; idT++) {
		// get 3D+c sub-array
		size_t lOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*vtDim_[3]*idT;
		hoNDArray<bool> habFullMask3D(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], habFullMask.get_data_ptr()+lOffset, false);
		hoNDArray<bool> habKSpaceCenter3D(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], habKSpaceCenter_.get_data_ptr()+lOffset, false);

		// temporal array
		hoNDArray<bool> habSamplPtrSym;

		// get center points
		int iCenterX = (int)(vtDim_[2]/2);
		int iCenterY = (int)(vtDim_[0]/2);
		int iCenterZ = (int)(vtDim_[1]/2);

		if (GlobalVar::instance()->cfLambdaESPReSSo_ != std::complex<float>(0.0)) {
			// ESPReSSo acquisition
			if (fPartialFourierVal_ < 1.0 && fPartialFourierVal_ > 0.0) {
				//-------------------------------------------------------------------------
				//------------------------- Upper / Lower ?  ------------------------------
				//------- check if upper or lower region in the k-space is sampled --------
				//-------------------------------------------------------------------------
				int iNumFoundUpper = 0;
				int iNumFoundLower = 0;

				// check only for one channel - sampling mask is same for all channels
				pbPtr_ = habFullMask3D.get_data_ptr();

				// ESPReSSo direction is phase encoding direction
				if (iESPReSSoDirection_ == 1) {
					for (int idY = (int)(fPartialFourierVal_*vtDim_[0])-1; static_cast<unsigned int>(idY) < vtDim_[0]; idY++) {
						if (pbPtr_[idY + vtDim_[0]*iCenterZ + vtDim_[0]*vtDim_[1]*iCenterX] != false) {
							iNumFoundUpper++;
						}
					}

					for (int idY = (int)(fPartialFourierVal_*vtDim_[0])-1; idY >= 0; idY--) {
						if (pbPtr_[idY + vtDim_[0]*iCenterZ + vtDim_[0]*vtDim_[1]*iCenterX] != false) {
							iNumFoundLower++;
						}
					}
				} else {
					// ESPReSSo direction is partition encoding direction
					for (int idZ = (int)(fPartialFourierVal_*vtDim_[1])-1; static_cast<unsigned int>(idZ) < vtDim_[1]; idZ++) {
						if (pbPtr_[iCenterY + idZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] != false) {
							iNumFoundUpper++;
						}
					}

					for (int idZ = (int)(fPartialFourierVal_*vtDim_[1])-1; idZ >= 0; idZ--) {
						if (pbPtr_[iCenterY + idZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] != false) {
							iNumFoundLower++;
						}
					}
				}

				bESPReSSoIsLower_ = (iNumFoundLower > iNumFoundUpper);

				//--------------------------------------------------------------------------
				//------------------------- symmetrical sampling pattern -------------------
				//- sampling mask is mirrored at the center line depending on ESPReSSo dir.-
				//--------------------------------------------------------------------------
				GINFO("compute symmetrical sampling pattern\n");

				habSamplPtrSym = habFullMask3D;
				pbPtr_ = habSamplPtrSym.get_data_ptr();

				// lower half sampled
				if (bESPReSSoIsLower_) {
					// ESPReSSo direction is phase encoding direction
					if (iESPReSSoDirection_ == 1) {
						#pragma omp parallel for
						for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
							for (size_t idX = 0; idX < vtDim_[2]; idX++) {
								for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
									for (size_t idY = vtDim_[0]-1; idY > vtDim_[0]/2; idY--) {
										pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
									}
								}
							}
						}
					} else {
						// ESPReSSo direction is partition encoding direction
						#pragma omp parallel for
						for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
							for (size_t idX = 0; idX < vtDim_[2]; idX++) {
								for (size_t idY = 0; idY < vtDim_[0]; idY++) {
									for (size_t idZ = vtDim_[1]-1; idZ > vtDim_[1]/2; idZ--) {
										pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr_[idY + (vtDim_[1]-1 - idZ)*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
									}
								}
							}
						}
					}
				} else {
					// upper half sampled
					if (iESPReSSoDirection_ == 1) {
						// ESPReSSo direction is phase encoding direction
						#pragma omp parallel for
						for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
							for (size_t idX = 0; idX < vtDim_[2]; idX++) {
								for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
									for (size_t idY = 0; idY < vtDim_[0]/2; idY++) {
										pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
									}
								}
							}
						}
					} else {
						// ESPReSSo direction is partition encoding direction
						#pragma omp parallel for
						for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
							for (size_t idX = 0; idX < vtDim_[2]; idX++) {
								for (size_t idY = 0; idY < vtDim_[0]; idY++) {
									for (size_t idZ = 0; idZ < vtDim_[1]/2; idZ++) {
										pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr_[idY + (vtDim_[1]-1 - idZ)*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
									}
								}
							}
						}
					}
				}

				//-------------------------------------------------------------------------
				//--------------------------- conjugate sampling pattern ------------------
				//---- conjugation - get only unsampled points for complex conjugation ----
				//-------------------------------------------------------------------------
				GINFO("compute conjugate sampling pattern\n");

				hoNDArray<bool> habMaskConj3D(habFullMask3D.get_dimensions());
				pbPtr_ = habMaskConj3D.get_data_ptr();
				pbPtr2_ = habFullMask3D.get_data_ptr();
				pbPtr3_ = habSamplPtrSym.get_data_ptr();

				// get unsampled points - XOR gating
				#pragma omp parallel for
				for (size_t i = 0; i < habMaskConj3D.get_number_of_elements(); i++) {
					pbPtr_[i] = pbPtr2_[i]^pbPtr3_[i];
				}
			} else {
				// ESPReSSo constraint for pure CS data without Partial Fourier technique

				//--------------------------------------------------------------------------
				//----------- ESPReSSo for Compressed Sensing w/o Partial Fourier ----------
				//--------------------------------------------------------------------------
				hoNDArray<bool> habMaskRight3D = habFullMask3D;
				pbPtr_ = habMaskRight3D.get_data_ptr();

				// set upper half zero
				#pragma omp parallel for
				for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
					for (size_t idX = 0; idX < vtDim_[2]; idX++) {
						for (size_t idY = 0; idY < vtDim_[0]/2; idY++) {
							for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
								pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = false;
							}
						}
					}
				}

				// set lower half zero
				hoNDArray<bool> habMaskLeft3D = habFullMask3D;
				pbPtr_ = habMaskLeft3D.get_data_ptr();

				#pragma omp parallel for
				for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
					for (size_t idX = 0; idX < vtDim_[2]; idX++) {
						for (size_t idY = vtDim_[0]/2; idY < vtDim_[0]; idY++) {
							for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
								pbPtr_[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = false;
							}
						}
					}
				}

				// logical OR with kSpaceCenter
				pbPtr_ = habMaskRight3D.get_data_ptr();
				pbPtr2_ = habKSpaceCenter3D.get_data_ptr();
				pbPtr3_ = habMaskLeft3D.get_data_ptr();

				#pragma omp parallel for
				for (size_t i = 0; i < habMaskRight3D.get_number_of_elements(); i++) {
					if (pbPtr2_[i] == true) {
						pbPtr_[i] = true;
						pbPtr3_[i] = true;
					}
				}

				// flip both masks (lmaskRight and lmaskLeft) and write to habMaskConj_ and habMaskConj2_
				hoNDArray<bool> habMaskConj3D(habMaskRight3D.get_dimensions());
				hoNDArray<bool> habMaskConj2_3D(habMaskRight3D.get_dimensions());
				bool *pcfPtr_habMaskConj3D = habMaskConj3D.get_data_ptr();
				bool *pcfPtr_habMaskConj2_3D= habMaskConj2_3D.get_data_ptr();

				#pragma omp parallel for
				for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
					for (size_t idX = 0; idX < vtDim_[2]; idX++) {
						for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
							for (size_t idY = 0; idY < vtDim_[0]; idY++) {
								pcfPtr_habMaskConj3D[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
								pcfPtr_habMaskConj2_3D[idY + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]] = pbPtr3_[(vtDim_[0]-1 - idY) + idZ*vtDim_[0] + idX*vtDim_[0]*vtDim_[1] + iCha*vtDim_[0]*vtDim_[1]*vtDim_[2]];
							}
						}
					}
				}

				pbPtr_ = habFullMask3D.get_data_ptr();
				pbPtr2_ = habMaskConj3D.get_data_ptr();
				pbPtr3_ = habMaskConj2_3D.get_data_ptr();

				// get unsampled points - XOR gating
				#pragma omp parallel for
				for (size_t i = 0; i < habMaskConj3D.get_number_of_elements(); i++) {
					pbPtr2_[i] = pbPtr_[i]^pbPtr2_[i] || pbPtr_[i]^pbPtr3_[i];
				}

				// fill 4D
				memcpy(habMaskConj_.get_data_ptr() + lOffset, habMaskConj3D.get_data_ptr(), sizeof(std::complex<float>)*habMaskConj3D.get_number_of_elements());

				// set values for symmetrical kSpace part / filter
				iESPReSSoDirection_ = 1;
				bESPReSSoIsLower_ = true;
			}

			//-------------------------------------------------------------------------
			//-------------------- symmetrical kSpace part / filter --------------------
			//--------------------------------------------------------------------------
			// get indices of center kSpace lines - find symmetrical kSpace part and calc filter
			GINFO("find symmetrical kSpace part..\n");

			// vectors for storing the outer sampled lines
			std::vector<int> viKSpaceLines_1;
			std::vector<int> viKSpaceLines_2;

			// loop over the symmetrical pattern beginning at the center line to the first unsampled point
			// --> difference between center line and first zero is half the symmetrical sampled lines
			// store the upper and lower symmetrical sampled line in the kSpaceLine vector

			// ESPReSSo is active
			if (fPartialFourierVal_ != 1.0) {
				pbPtr_ = habSamplPtrSym.get_data_ptr();

				// phase encoding direction
				if (iESPReSSoDirection_ == 1) {
					for (size_t iZ = 0; iZ < vtDim_[1]; iZ++) {
						// lower region is sampled
						if (bESPReSSoIsLower_ == true) {
							if (fPartialFourierVal_ != .5) {
								viKSpaceLines_1.push_back(fPartialFourierVal_*vtDim_[0]-iCenterY);
								viKSpaceLines_2.push_back(fPartialFourierVal_*vtDim_[0]);
							} else {
								for (int iY = 0; iY < iCenterY; iY++) {
									if (pbPtr_[iY + iZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] == true) {
										// number of line in k-space
										viKSpaceLines_1.push_back(iY);

										// number of "mirrored" line
										viKSpaceLines_2.push_back(vtDim_[0]-iY);

										break;
									} else if (iY == iCenterY-1) {
										// number of line in k-space
										viKSpaceLines_1.push_back(iY-1);

										// number of "mirrored" line
										viKSpaceLines_2.push_back(vtDim_[0]-iY+1);
									}
								}
							}
						} else {
							// upper region is sampled
							if (fPartialFourierVal_ != .5) {
								viKSpaceLines_1.push_back(fPartialFourierVal_*vtDim_[0]-iCenterY);
								viKSpaceLines_2.push_back(fPartialFourierVal_*vtDim_[0]);
							} else {
								for (int iY = vtDim_[0]; iY > iCenterY; iY--) {
									if (pbPtr_[iY + iZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] == true) {
										// number of line in k-space
										viKSpaceLines_1.push_back(vtDim_[0]-iY);

										// number of "mirrored" line
										viKSpaceLines_2.push_back(iY);

										break;
									} else if (iY == iCenterY+1) {
										// number of line in k-space
										viKSpaceLines_1.push_back(iY-1);

										// number of "mirrored" line
										viKSpaceLines_2.push_back(vtDim_[0]-iY+1);
									}
								}
							}
						}
					}
				} else {
					// partition encoding direction
					for (size_t iY = 0; iY < vtDim_[0]; iY++) {
						// lower region is sampled
						if (bESPReSSoIsLower_ == true) {
							if (fPartialFourierVal_ != .5) {
								viKSpaceLines_1.push_back(fPartialFourierVal_*vtDim_[1]-iCenterZ);
								viKSpaceLines_2.push_back(fPartialFourierVal_*vtDim_[1]);
							} else {
								for (int iZ = 0; iZ < iCenterZ; iZ++) {
									if (pbPtr_[iY + iZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] == true) {
										// number of line in k-space
										viKSpaceLines_1.push_back(iZ);

										// number of "mirrored" line
										viKSpaceLines_2.push_back(vtDim_[1]-iZ);

										break;
									} else if (iZ==iCenterZ-1) {
										// number of line in k-space
										viKSpaceLines_1.push_back(iZ-1);

										// number of "mirrored" line
										viKSpaceLines_2.push_back(vtDim_[1]-iZ+1);
									}
								}
							}
						} else {
							// upper region is sampled
							if (fPartialFourierVal_ != .5) {
								viKSpaceLines_1.push_back(fPartialFourierVal_*vtDim_[1]-iCenterZ);
								viKSpaceLines_2.push_back(fPartialFourierVal_*vtDim_[1]);
							} else {
								for (int iZ = vtDim_[1]; iZ > iCenterZ; iZ--) {
									if (pbPtr_[iY + iZ*vtDim_[0] + iCenterX*vtDim_[0]*vtDim_[1]] == true) {
										// number of line in k-space
										viKSpaceLines_1.push_back(vtDim_[1]-iZ);

										// number of "mirrored" line
										viKSpaceLines_2.push_back(iZ);

										break;
									} else if (iZ==iCenterZ+1) {
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
			} else {
				// additional complex conjugate symmetry constraint for normal Compressed Sensing data
				int iLine1 = .75*vtDim_[0]-iCenterY;
				int iLine2 = .75*vtDim_[0];

				for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
					viKSpaceLines_1.push_back(iLine1);
					viKSpaceLines_2.push_back(iLine2);
				}
			}

			// if no line is found --> filter width will be 0 --> take default values from window calculation
			if (viKSpaceLines_1.size() == 0 || viKSpaceLines_2.size()==0) {
				viKSpaceLines_1.clear(); viKSpaceLines_2.clear();

				GINFO("Estimate filter size by calculated window\n");

				// phase encoding direction
				if (iESPReSSoDirection_ == 1) {
					for (size_t iZ = 0; iZ < vtDim_[1]; iZ++) {
						int iLine1 = .75*vtDim_[0]-iCenterY;
						int iLine2 = .75*vtDim_[0];
						viKSpaceLines_1.push_back(iLine1);
						viKSpaceLines_2.push_back(iLine2);
					}
				} else {
					// partition encoding direction
					for (size_t iY = 0; iY < vtDim_[0]; iY++) {
						int iLine1 = .75*vtDim_[1]-iCenterZ;
						int iLine2 = .75*vtDim_[1];
						viKSpaceLines_1.push_back(iLine1);
						viKSpaceLines_2.push_back(iLine2);
					}
				}
			}

			GINFO("calculate filter..\n");

			// create filter array (inital all ones)
			hoNDArray< std::complex< float > > hacfFilter3D(habFullMask3D.get_dimensions());
			hacfFilter3D.fill(std::complex<float>(1.0));
			hoNDArray<std::complex<float> > hacfFilt1D(vtDim_[0], vtDim_[1]);
			hoNDArray<std::complex<float> > hacfFilt1DFullArray(habFullMask3D.get_dimensions());

			// loop over the 3 spatial dimensions
			for (int iDim = 0; iDim < 3; iDim++) {
				hacfFilt1D.fill(std::complex<float>(0.0)); hacfFilt1DFullArray.fill(std::complex<float>(0.0));

				// filter for ESPReSSo direction (Hanning)
				if (iDim == iESPReSSoDirection_) {
					// phase encoding direction
					if (iESPReSSoDirection_ == 1) {
						for (size_t iZ = 0; iZ < vtDim_[1]; iZ++) {
							// get filter coefficients (passed parameter is window width)
							std::vector<float> *HanningCoeff = fGetHanningWindow(viKSpaceLines_2.at(iZ)-viKSpaceLines_1.at(iZ));

							// place coefficients in the filt1D array - loop over filter "lines"
							size_t iI;
							int iY;
							for (iI = 0, iY = viKSpaceLines_1.at(iZ); iY < viKSpaceLines_2.at(iZ); iY++, iI++) {
								hacfFilt1D(iY, iZ) = std::complex<float>(HanningCoeff->at(iI));
							}

							delete HanningCoeff;
							HanningCoeff = NULL;
						}
					} else if (iESPReSSoDirection_ == 2)
						// partition encoding direction
						for (size_t iY = 0; iY < vtDim_[0]; iY++) {
							std::vector<float> *vfHanningCoeff = fGetHanningWindow(viKSpaceLines_2.at(iY)-viKSpaceLines_1.at(iY));

							size_t iI;
							int iZ;
							for (iI = 0, iZ = viKSpaceLines_1.at(iY); iZ < viKSpaceLines_2.at(iY); iZ++, iI++) {
								hacfFilt1D(iY, iZ) = std::complex<float>(vfHanningCoeff->at(iI));
							}

							delete vfHanningCoeff;
							vfHanningCoeff = NULL;
						}

					// repetition in x-direction for getting a full 3D array
					for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
						for (size_t iX = 0; iX < vtDim_[2]; iX++) {
							for (size_t iY = 0; iY < vtDim_[0]; iY++) {
								for (size_t iZ = 0; iZ < vtDim_[1]; iZ++) {
									hacfFilt1DFullArray(iY, iZ, iX, iCha) = hacfFilt1D(iY, iZ);
								}
							}
						}
					}
				} else {
					// filter for other direction (Hamming) - same procedure as in ESPReSSo direction

					// x-direction
					if (iDim == 0) {
						std::vector<float> *vfHammingCoeff = fGetHammingWindow(vtDim_[2]);

						for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
							for (size_t idX = 0; idX < vtDim_[2]; idX++) {
								for (size_t idY = 0; idY < vtDim_[0]; idY++) {
									for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
										hacfFilt1DFullArray(idY, idZ, idX, iCha) = std::complex<float>(vfHammingCoeff->at(idX));
									}
								}
							}
						}

						delete vfHammingCoeff;
						vfHammingCoeff = NULL;
					} else if (iDim == 1) {
						// y-direction
						std::vector<float> *vfHammingCoeff = fGetHammingWindow(vtDim_[0]);

						for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
							for (size_t idX = 0; idX < vtDim_[2]; idX++) {
								for (size_t idY = 0; idY < vtDim_[0]; idY++) {
									for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
										hacfFilt1DFullArray(idY, idZ, idX, iCha) = std::complex<float>(vfHammingCoeff->at(idY));
									}
								}
							}
						}

						delete vfHammingCoeff;
						vfHammingCoeff = NULL;
					} else if (iDim == 2) {
						// z-direction
						std::vector<float> *vfHammingCoeff = fGetHammingWindow(vtDim_[1]);

						for (size_t iCha = 0; iCha < vtDim_[3]; iCha++) {
							for (size_t idX = 0; idX < vtDim_[2]; idX++) {
								for (size_t idY = 0; idY < vtDim_[0]; idY++) {
									for (size_t idZ = 0; idZ < vtDim_[1]; idZ++) {
										hacfFilt1DFullArray(idY, idZ, idX, iCha) = std::complex<float>(vfHammingCoeff->at(idZ));
									}
								}
							}
						}

						delete vfHammingCoeff;
						vfHammingCoeff = NULL;
					}
				}

				hacfFilter3D *= hacfFilt1DFullArray;
			}

			// check maximum value in filter array
			pcfPtr_ = hacfFilter3D.get_data_ptr();
			for (size_t iI = 0; iI < hacfFilter3D.get_number_of_elements(); iI++) {
				if (abs(pcfPtr_[iI]) > 1.0) {
					pcfPtr_[iI] = std::complex<float>(0.0);
				}
			}

			memcpy(hacfFilter_.get_data_ptr() + lOffset, hacfFilter3D.get_data_ptr(), sizeof(std::complex<float>)*hacfFilter3D.get_number_of_elements());
		}

		// permute array: y-z-x-c-t -> t-y-z-x-c
		vtDimOrder.clear();
		vtDimOrder.push_back(4);
		vtDimOrder.push_back(0);
		vtDimOrder.push_back(1);
		vtDimOrder.push_back(2);
		vtDimOrder.push_back(3);
		vtDim_ = *habKSpaceCenter_.get_dimensions();
	}
}

//--------------------------------------------------------------------------
//---------------------------- windowing -----------------------------------
//--------------------------------------------------------------------------
void CS_FOCUSS_4D::fWindowing(hoNDArray<std::complex<float> > &hacfWWindowed)
{
	// permute array: t-y-z-x-c -> y-z-x-c-t
	std::vector<size_t> vtDimOrder;
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(3);
	vtDimOrder.push_back(4);
	vtDimOrder.push_back(0);
	hacfWWindowed = *permute(&hacfWWindowed, &vtDimOrder,false);

	vtDim_ = *hacfWWindowed.get_dimensions();

	habKSpaceCenter_.create(hacfWWindowed.get_dimensions());
	habKSpaceCenter_.fill(true);

	// loop over t
	for (size_t idT = 0; idT < vtDim_[4]; idT++) {
		// get 3D+c sub-array
		long lOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*vtDim_[3]*idT;
		hoNDArray<std::complex<float> > hacfWWindowed3D(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], hacfWWindowed.get_data_ptr()+lOffset, false);

		// array with mask
		hoNDArray<std::complex<float> > hacfMask3D(hacfWWindowed3D.get_dimensions());
		hacfMask3D.fill(std::complex<float>(0.0));

		// get calibration mask
		std::vector<size_t> vStart, vSize;
		for (size_t iI = 0; iI < 2; iI++) {
			if (viCalibrationSize_.at(iI + idT*4) % 2) {
				vStart.push_back(std::floor((float)vtDim_[iI]/2)+std::ceil(-(float)viCalibrationSize_.at(iI + idT*4)/2));
			} else {
				vStart.push_back(std::floor((float)vtDim_[iI]/2)+std::ceil(-(float)viCalibrationSize_.at(iI + idT*4)/2));
			}
		}

		vStart.push_back(0);
		vStart.push_back(0);

		for (size_t iY = vStart.at(0); iY < vStart.at(0)+viCalibrationSize_.at(0 + idT*4); iY++) {
			for (size_t iZ = vStart.at(1); iZ < vStart.at(1)+viCalibrationSize_.at(1 + idT*4); iZ++) {
				for (size_t iX = vStart.at(2); iX < vStart.at(2)+viCalibrationSize_.at(2 + idT*4); iX++) {
					for (size_t iC = vStart.at(3); iC < vStart.at(3)+viCalibrationSize_.at(3 + idT*4); iC++) {
						hacfMask3D(iY, iZ, iX, iC) = std::complex<float>(1.0);
					}
				}
			}
		}

		// windowing W
		hacfWWindowed3D *= hacfMask3D;

		// get kSpaceCenter
		hoNDArray<bool> habKSpaceCenter3D(hacfWWindowed3D.get_dimensions());
		habKSpaceCenter3D.fill(true);
		pbPtr_ = habKSpaceCenter3D.get_data_ptr();
		pcfPtr_ = hacfMask3D.get_data_ptr();

		#pragma omp parallel for
		for (size_t lI = 0; lI < hacfMask3D.get_number_of_elements(); lI++) {
			if (pcfPtr_[lI] == std::complex<float>(0.0)) {
				pbPtr_[lI] = false;
			}
		}

		GINFO("data windowed for initial estimate and kSpaceCenter found..\n");

		memcpy(habKSpaceCenter_.get_data_ptr() + lOffset, habKSpaceCenter3D.get_data_ptr(), sizeof(bool)*habKSpaceCenter3D.get_number_of_elements());
		memcpy(hacfWWindowed.get_data_ptr() + lOffset, hacfWWindowed3D.get_data_ptr(), sizeof(std::complex<float>)*hacfWWindowed3D.get_number_of_elements());
	}

	// permute array: y-z-x-c-t -> t-y-z-x-c
	vtDimOrder.clear();
	vtDimOrder.push_back(4);
	vtDimOrder.push_back(0);
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(3);
	hacfWWindowed = *permute(&hacfWWindowed, &vtDimOrder, false);

	vtDim_ = *hacfWWindowed.get_dimensions();
}

//--------------------------------------------------------------------------
//---------------------- get calibration size ------------------------------
//--------------------------------------------------------------------------
void CS_FOCUSS_4D::fGetCalibrationSize(hoNDArray<bool> &habArray)
{
	// permute array: t-y-z-x-c -> y-z-x-c-t
	std::vector<size_t> vtDimOrder;
	vtDimOrder.push_back(1);
	vtDimOrder.push_back(2);
	vtDimOrder.push_back(3);
	vtDimOrder.push_back(4);
	vtDimOrder.push_back(0);
	habArray = *permute(&habArray, &vtDimOrder,false);

	vtDim_ = *habArray.get_dimensions();

	// loop over t
	for (size_t idT = 0; idT < vtDim_[4]; idT++) {
		size_t lOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*vtDim_[3]*idT;
		hoNDArray< bool > hoNDArray3D(vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3], habArray.get_data_ptr() + lOffset, false);

		size_t iSY = 2;
		size_t iSZ = 2;

		bool bYflag = false;
		bool bZflag = false;

		std::vector<size_t> vtDim = *hoNDArray3D.get_dimensions();

		hoNDArray<bool> habMask_helper;

		while (!(bYflag && bZflag)) {
			if (!bYflag) {
				fCropArrYZ(habArray, iSY+1, iSZ, habMask_helper);

				if (fAllOne(habMask_helper)) {
					iSY++;
				} else {
					bYflag = true;
				}
			}

			if (!bZflag) {
				fCropArrYZ(habArray, iSY, iSZ+1, habMask_helper);

				if (fAllOne(habMask_helper)) {
					iSZ++;
				} else {
					bZflag = true;
				}
			}

			if (iSY == vtDim[0]) {
				bYflag = true;
			}

			if (iSZ == vtDim[1]) {
				bZflag = true;
			}
		}

		// push values on calibration size vector
		viCalibrationSize_.push_back(iSY);
		viCalibrationSize_.push_back(iSZ);
		viCalibrationSize_.push_back(vtDim[2]);
		viCalibrationSize_.push_back(vtDim[3]);

		for (std::vector<int>::size_type i = 0; i != viCalibrationSize_.size(); i++) {
			GDEBUG("calibration size: %i..\n", viCalibrationSize_.at(i));
		}
	}
}

GADGET_FACTORY_DECLARE(CS_FOCUSS_4D)
