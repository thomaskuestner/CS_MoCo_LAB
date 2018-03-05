#include "CS_Retro_NavigatorGadget.h"

#include <cmath>

using namespace Gadgetron;

// class constructor
CS_Retro_NavigatorGadget::CS_Retro_NavigatorGadget()
{
}

// class destructor - delete temporal buffer/memory
CS_Retro_NavigatorGadget::~CS_Retro_NavigatorGadget()
{
}

// read flexible data header
int CS_Retro_NavigatorGadget::process_config(ACE_Message_Block *mb)
{
	return GADGET_OK;
}

int CS_Retro_NavigatorGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1,GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m3)
{
	// get gadget property
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	iNavMethod_ = NavigationMethod.value();
#else
	iNavMethod_ = *(get_int_value("NavigationMethod").get());
#endif

	// fetch attribute values from header
	iNoChannels_ = m1->getObjectPtr()->channels;
	iNoNav_		 = m1->getObjectPtr()->user_int[5];
	lNoScans_	 = m3->getObjectPtr()->get_size(1);

	field_of_view_[0] = m1->getObjectPtr()->field_of_view[0];
	field_of_view_[1] = m1->getObjectPtr()->field_of_view[1];
	field_of_view_[2] = m1->getObjectPtr()->field_of_view[2];

	// get navigator signal according to selected method
	// 0: classical
	// 1: PCA
	switch (iNavMethod_) {
	case 0:
		getNav2D(*m2->getObjectPtr());
		break;

	case 1:
		getNav2DPCA(*m2->getObjectPtr());
		break;

	default:
		GERROR("Navigation method %d unknown! Please specify one via gadget property.\n", iNavMethod_);
		return GADGET_FAIL;
	}

	GadgetContainerMessage<hoNDArray<float> > *tmp_m2 = new GadgetContainerMessage<hoNDArray<float> >();
	tmp_m2->getObjectPtr()->create(vNavInt_.size());

	// convert vector to array
	float *fPtr = tmp_m2->getObjectPtr()->get_data_ptr();
	for (size_t iI = 0; iI < vNavInt_.size(); iI++) {
		fPtr[iI] = vNavInt_.at(iI);
	}

	m1->cont(tmp_m2);
	tmp_m2->cont(m3);

	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

// get interpolated navigator signal
void CS_Retro_NavigatorGadget::getNav2D(hoNDArray<std::complex<float> > &aNav)
{
	GDEBUG("\n\n**************************************\n********** get navigator 2D **********\n**************************************\n\n");

	// reconstruct the 1-D projections for all measurements and all channels
	GINFO("domain transformation - k-space to image\n");

	/* MATLAB
	% Reconstruct the 1-D projections for all measurements and all channels
	dImg = fftshift(ifft(ifftshift(dKSpace)));
	dImg = fftshift(ifft(ifftshift(dImg, 3), [], 3), 3);
	dImg = flipdim(dImg, 1); % Invert the RO direction: 1-N -> H-F
	dImg = dImg(iNSamples/4:iNSamples.*3/4 - 1, :, :, :); % RO x t x PE x CH
	*/
	hoNDArray<std::complex<float> > aImg = aNav;
	hoNDFFT_CS<float>::instance()->ifftshift3D(aImg);
	hoNDFFT_CS<float>::instance()->ifft1(aImg);
	hoNDFFT_CS<float>::instance()->fftshift3D(aImg);
	hoNDFFT_CS<float>::instance()->ifft(&aImg, 2, true);
	flip_array(aImg, 0);

	// crop center part due to twofold oversampling
	GINFO("crop center part of two-fold oversampled data..\n");

	std::vector<size_t> vStart, vSize;

	vStart.push_back(static_cast<size_t>(std::floor(0.25*aNav.get_size(0)))-1);
	vStart.push_back(0);
	vStart.push_back(0);
	vStart.push_back(0);

	vSize.push_back(std::floor(0.5*aNav.get_size(0)));
	vSize.push_back(aNav.get_size(1));
	vSize.push_back(aNav.get_size(2));
	vSize.push_back(aNav.get_size(3));

	get_subarray(aImg, vStart, vSize, aImg);

	// get channel power and normalize channel images
	GINFO("calculate channel power and normalize channel images..\n");

	std::vector<float> fPower(iNoChannels_);
	hoNDArray<std::complex<float> > aPower = aImg;
	aPower.fill(std::complex<float>(0,0));
	for (int c = 0; c < iNoChannels_; c++) {
		size_t offset = aImg.get_size(0)*aImg.get_size(1)*aImg.get_size(2)*c;
		hoNDArray<std::complex<float> > SubArray(aImg.get_size(0), aImg.get_size(1), aImg.get_size(2), aImg.get_data_ptr()+offset, false);
		fPower.at(c) = asum(&SubArray);

		// fill part of the 3D array
		#pragma omp parallel for
		for (size_t i = 0; i < aImg.get_size(0)*aImg.get_size(1)*aImg.get_size(2); i++) {
			aPower.at(i+offset) = std::complex<float>(fPower.at(c), fPower.at(c));
		}
	}

	divide(aImg, aPower, aImg);

	GINFO("channel images normalized..\n");

	aPower.clear();

	// get x range of respiratory motion & FFT without scrambling
	hoNDArray<std::complex<float> > aFreq = aImg;
	hoNDFFT_CS<float>::instance()->fft(&aFreq, 1, false);
	multiplyConj(aFreq,aFreq,aPower);

	// conversion from complex float to float
	hoNDArray<float> afPower(aPower.get_dimensions());
	for (size_t i = 0; i < afPower.get_number_of_elements(); i++) {
		afPower[i] = aPower[i].real();
	}

	/*
	dIMGres = 1./(double(dNavPeriod)./1000.*double(iNMeasurements)); % The frequency resolution of dIMG in Hz
	dPower = squeeze(sum(dPower(:, round(1./(5.*dIMGres)):round(1./(3.*dIMGres)), :, :), 2)); % RO x PE x CH
	*/
	float fIMGRes = 1.0/((static_cast<float>(GlobalVar::instance()->iNavPeriod_)/1000.0)*static_cast<float>(iNoNav_)); // frequency resolution of aImg in Hz
	hoNDArray<float> aPowerInChan, aPowerAcrossChan;

	vStart.clear();
	vStart.push_back(0);
	vStart.push_back(std::floor(1.0/(5*fIMGRes)-.5));
	vStart.push_back(0);
	vStart.push_back(0);

	vSize.clear();
	vSize.push_back(afPower.get_size(0));
	vSize.push_back(std::ceil(1.0/(3*fIMGRes)-.5)-std::ceil(1.0/(5*fIMGRes)-.5)+1);
	vSize.push_back(afPower.get_size(2));
	vSize.push_back(afPower.get_size(3));

	get_subarray(afPower, vStart, vSize, afPower);

	sum_dim(afPower, 1, afPower); // RO x PE x CH
	sum_dim(afPower, 1, aPowerInChan); // RO x CH
	sum_dim(aPowerInChan, 1, aPowerAcrossChan); // RO x 1

	// Prevent detection of regions in the abdomen
	for (long i = aPowerAcrossChan.get_number_of_elements(); i > std::floor(aPowerAcrossChan.get_number_of_elements()*.75); i--) {
		aPowerAcrossChan.at(i) = 0;
	}

	// get Gaussian filter kernel
	std::vector<float> vGaussian;
	filter1DGaussian(vGaussian, 20);
	arrayConv(aPowerAcrossChan, vGaussian, 0);

	// find index of maximum
	int iMaxIndex = amax(&aPowerAcrossChan);

	GINFO("data filtered and maximum determined.. iMaxIndex: %i\n", iMaxIndex);

	if ((iMaxIndex < 20) || (iMaxIndex > static_cast<int>(aPowerInChan.get_size(0))-20)) {
		GERROR("Error: iMaxIndex out of bounds..\n");

		return;
	}

	//-------------------------------------------------------------------------
	// sort out channels with no relevant information in target area
	GINFO("get channels which contain most information..\n");

	vStart.clear();
	vStart.push_back(iMaxIndex-20);
	vStart.push_back(0);

	vSize.clear();
	vSize.push_back(41);
	vSize.push_back(iNoChannels_);

	get_subarray(aPowerInChan, vStart, vSize,aPowerInChan);

	GINFO("41 elements around maximum cropped..\n");

	std::vector<float> vGoodChannels;

	// loop over channels
	for (int c = 0; c < iNoChannels_; c++) {
		hoNDArray<float> aTmp;
		size_t offset = aPowerInChan.get_size(0)*c;
		aTmp.create(aPowerInChan.get_size(0), aPowerInChan.get_data_ptr()+offset, false);
		int iTmpInd = amax(&aTmp);
		vGoodChannels.push_back(aTmp.at(iTmpInd));
	}

	int iIndex = std::max_element(vGoodChannels.begin(), vGoodChannels.end())- vGoodChannels.begin();
	float fMax = vGoodChannels.at(iIndex);
	int iNumGood = 0;

	for (size_t i = 0; i < vGoodChannels.size(); i++) {
		if (vGoodChannels.at(i) > .2*fMax) {
			vGoodChannels.at(i) = 1;
			iNumGood++;
		} else {
			vGoodChannels.at(i) = 0;
		}
	}

	for (size_t i = 0; i < vGoodChannels.size(); i++) {
		GDEBUG("vGoodChannels[%i]: %f\n", i,vGoodChannels.at(i));
	}

	//-------------------------------------------------------------------------
	// get the best PE line
	// get sub array of good channels
	GDEBUG("get best PE line..\n");

	hoNDArray<float> aPowerInPE(41,afPower.get_size(1), iNumGood);

	vStart.clear();
	vStart.push_back(iMaxIndex-20);
	vStart.push_back(0);
	vStart.push_back(0);

	vSize.clear();
	vSize.push_back(41);
	vSize.push_back(afPower.get_size(1));
	vSize.push_back(1);

	GDEBUG("vStart: %i, %i, %i, vSize: %i, %i, %i, afPower size: %i, %i, %i\n", vStart.at(0), vStart.at(1), vStart.at(2), vSize.at(0), vSize.at(1), vSize.at(2), afPower.get_size(0), afPower.get_size(1), afPower.get_size(2));

	size_t o = 0; //helper - only non-zero entries of vGoodChannels are of interest
	for (int c = 0; c < iNoChannels_; c++) {
		if (vGoodChannels.at(c) == 1) {
			hoNDArray<float> aTmp;
			aTmp.create(&vSize);
			vStart.at(2) = c;
			get_subarray(afPower, vStart, vSize, aTmp);

			// fill array
			size_t offset = aPowerInPE.get_size(0)*aPowerInPE.get_size(1)*o;
			for (size_t i = 0; i < aTmp.get_number_of_elements(); i++) {
				aPowerInPE.at(i + offset) = aTmp.at(i);
			}

			o++;
		}
	}

	sum_dim(aPowerInPE, 2, aPowerInPE);
	sum_dim(aPowerInPE, 0, aPowerInPE);

	GINFO("\n aPowerInPE\n");
	aPowerInPE.print(std::cout);

	iMaxIndex = amax(&aPowerInPE);
	int iMaxChan = iMaxIndex;

	GINFO("found at %i\n", iMaxIndex);

	//-------------------------------------------------------------------------
	// find best corresponding channels according to best phase encoding position
	GINFO("searching for channels according to best PE line..\n");

	// get sub array
	vStart.clear();
	vStart.push_back(0);
	vStart.push_back(0);
	vStart.push_back(iMaxIndex);
	vStart.push_back(0);

	vSize.clear();
	vSize.push_back(aImg.get_size(0));
	vSize.push_back(aImg.get_size(1));
	vSize.push_back(1);
	vSize.push_back(aImg.get_size(3));

	aFreq.clear();

	get_subarray(aImg, vStart, vSize, aFreq);

	hoNDFFT_CS<float>::instance()->fft(&aFreq, 1, false);
	multiplyConj(aFreq,aFreq,aPower);

	// conversion from complex float to float
	afPower.clear();
	afPower.create(*aPower.get_dimensions());
	for (size_t i = 0; i < afPower.get_number_of_elements(); i++) {
		afPower[i] = aPower[i].real();
	}

	fIMGRes = 1.0/((static_cast<float>(GlobalVar::instance()->iNavPeriod_)/1000.0)*static_cast<float>(iNoNav_)); // frequency resolution of aImg in Hz

	aPowerInChan.clear();
	aPowerAcrossChan.clear();

	vStart.clear();
	vStart.push_back(0);
	vStart.push_back(std::ceil(1.0/(5*fIMGRes)-.5));
	vStart.push_back(0);
	vStart.push_back(0);

	vSize.clear();
	vSize.push_back(afPower.get_size(0));
	vSize.push_back(std::ceil(1.0/(3*fIMGRes)-.5)-std::ceil(1.0/(5*fIMGRes)-.5)+1);
	vSize.push_back(afPower.get_size(2));
	vSize.push_back(afPower.get_size(3));

	GDEBUG("vSize: %i, %i, %i, %i - fIMGRes: %f\n", vSize.at(0), vSize.at(1), vSize.at(2), vSize.at(3), fIMGRes);

	get_subarray(afPower, vStart, vSize, afPower);
	sum_dim(afPower, 1, afPower); // RO x PE x CH
	sum_dim(afPower, 1, aPowerInChan); // RO x CH
	sum_dim(aPowerInChan, 1, aPowerAcrossChan); // RO x 1

	// Prevent detection of regions in the abdomen
	for (long i = aPowerAcrossChan.get_number_of_elements(); i > std::floor(aPowerAcrossChan.get_number_of_elements()*.75); i--) {
		aPowerAcrossChan.at(i) = 0;
	}

	GINFO("filter data with Gaussian kernel..\n");

	// get Gaussian filter kernel
	vGaussian.clear();
	filter1DGaussian(vGaussian, 20);
	arrayConv(aPowerAcrossChan, vGaussian);

	// find index of maximum
	iMaxIndex = amax(&aPowerAcrossChan);
	int dX = iMaxIndex;

	GINFO("found at %i\n", iMaxIndex);

	// get good channels
	vStart.clear();
	vStart.push_back(iMaxIndex-20);
	vStart.push_back(0);

	vSize.clear();
	vSize.push_back(41);
	vSize.push_back(aPowerInChan.get_size(1));

	get_subarray(aPowerInChan, vStart, vSize,aPowerInChan);

	vGoodChannels.clear();

	// loop over channels
	for (int c = 0; c < iNoChannels_; c++) {
		size_t offset = aPowerInChan.get_size(0)*c;

		hoNDArray<float> aTmp;
		aTmp.create(aPowerInChan.get_size(0), aPowerInChan.get_data_ptr()+offset, false);

		int iTmpInd = amax(&aTmp);
		vGoodChannels.push_back(aTmp.at(iTmpInd));
	}

	iIndex = std::max_element(vGoodChannels.begin(), vGoodChannels.end())-vGoodChannels.begin();
	fMax =vGoodChannels.at(iIndex);
	iNumGood = 0;

	for (size_t i = 0; i < vGoodChannels.size(); i++) {
		if (vGoodChannels.at(i) > .2*fMax) {
			vGoodChannels.at(i) = 1;
			iNumGood++;
		} else {
			vGoodChannels.at(i) = 0;
		}
	}

	for (size_t i = 0; i < vGoodChannels.size(); i++) {
		GDEBUG("vGoodChannels[%i]: %f\n", i,vGoodChannels.at(i));
	}

	// get relevant image - dRelevantImg = squeeze(dImg(:, :, dPos, lGoodChannels));
	vStart.clear();
	vStart.push_back(0);
	vStart.push_back(0);
	vStart.push_back(iMaxChan);
	vStart.push_back(0);

	vSize.clear();
	vSize.push_back(aImg.get_size(0));
	vSize.push_back(aImg.get_size(1));
	vSize.push_back(1);
	vSize.push_back(1);

	hoNDArray<std::complex<float> > acRelevantImg(aImg.get_size(0), aImg.get_size(1), iNumGood);

	o = 0; //helper - only non-zero entries of vGoodChannels are of interest
	hoNDArray<std::complex<float> > afTmp;
	for (int c = 0; c < iNoChannels_; c++) {
		if (vGoodChannels.at(c) == 1) {
			vStart.at(3) = c;
			get_subarray(aImg, vStart, vSize, afTmp);

			// fill array
			size_t offset = acRelevantImg.get_size(0)*acRelevantImg.get_size(1)*o;
			for (size_t i = 0; i < afTmp.get_number_of_elements(); i++) {
				acRelevantImg.at(i + offset) = afTmp.at(i);
			}

			o++;
		}
	}

	//-------------------------------------------------------------------------
	GINFO("get SOSImg\n");

	hoNDArray<float> aRelevantImg(*acRelevantImg.get_dimensions());
	multiplyConj(acRelevantImg,acRelevantImg,acRelevantImg);
	// complex float to float datatype
	for (size_t i = 0; i < aRelevantImg.get_number_of_elements(); i++) {
		aRelevantImg[i] = acRelevantImg[i].real();
	}

	// loop over channels
	size_t tOffset = aRelevantImg.get_size(0)*aRelevantImg.get_size(1);
	hoNDArray<float> hafMax(aRelevantImg.get_dimensions());
	float *fPtr = hafMax.get_data_ptr();
	for (size_t iI = 0; iI < aRelevantImg.get_size(2); iI++) {
		hoNDArray<float> tmp(aRelevantImg.get_size(0), aRelevantImg.get_size(1), aRelevantImg.get_data_ptr()+tOffset*iI);
		iMaxIndex = amax(&tmp);
		fMax = tmp.at(iMaxIndex);

		for (size_t iL = 0; iL < tOffset; iL++) {
			fPtr[iL+tOffset*iI] = fMax;
		}
	}

	divide(aRelevantImg, hafMax, aRelevantImg);

	hoNDArray<float> aSOSImg;
	sum_dim(aRelevantImg, 2, aSOSImg);
	sqrt(aSOSImg, aSOSImg);

	// convert array
	hoNDArray<std::complex<float> > cfaSOSImgTest(aSOSImg.get_dimensions());
	cfaSOSImgTest.fill(std::complex<float>(0.0, 0.0));
	std::complex<float> *cfPointer = cfaSOSImgTest.get_data_ptr();
	for (size_t iI = 0; iI < cfaSOSImgTest.get_number_of_elements(); iI++) {
		cfPointer[iI] = std::complex<float>(aSOSImg.at(iI), 0.0);
	}

	//-------------------------------------------------------------------------
	// get navigator
	int iDisplacementMax = 80; // [mm] diaphragm displacement (max. +/- 80mm)
	int iDisplacement = static_cast<int>(std::ceil(static_cast<float>(iDisplacementMax)/(static_cast<float>(field_of_view_[0])/static_cast<float>(aImg.get_size(0))) - 0.5));

	// fill last line in ref image
	hoNDArray<float> aRefImg = aSOSImg;
	aRefImg.fill(0.0);
	for (size_t iI = 0; iI < aSOSImg.get_size(0); iI++) {
		aRefImg.at(iI + aSOSImg.get_size(0)*(aSOSImg.get_size(1)-1)) = aSOSImg.at(iI + aSOSImg.get_size(0)*(aSOSImg.get_size(1)-1));
	}

	// create index vector and navigator vector (filled with zeros)
	std::vector<int> vIdx;
	for (size_t i = 0; i < aRefImg.get_size(1); i++) {
		vIdx.push_back(i);
	}

	vNav_.clear();
	for (size_t i = 0; i < aRefImg.get_size(1); i++) {
		vNav_.push_back(0);
	}

	hoNDArray<float> aRMSImg;
	for (int i = aRefImg.get_size(1)-2; i > 1; i--) {// -1 ; i--){
		aRMSImg.create(aRefImg.get_size(0), 2*iDisplacement+1);
		aRMSImg.fill(0.0);

		//MATLAB: tmp = dRefImg(:,idx(i+1:end))
		hoNDArray<float> tmp;

		vStart.clear();
		vStart.push_back(0);
		vStart.push_back(i+1);

		vSize.clear();
		vSize.push_back(aRefImg.get_size(0));
		vSize.push_back(aRefImg.get_size(1)-i-1);

		get_subarray(aRefImg, vStart, vSize, tmp);

		hoNDArray<float> repTmp(tmp.get_dimensions());

		//MATLAB: dSOSImg(:,idx(i)
		hoNDArray<float> aTmp;

		vStart.clear();
		vStart.push_back(0);
		vStart.push_back(vIdx.at(i));

		vSize.clear();
		vSize.push_back(aSOSImg.get_size(0));
		vSize.push_back(1);

		get_subarray(aSOSImg, vStart, vSize, aTmp);

		hoNDArray<float> aTmp2 = aTmp;
		circshift(aTmp2, -iDisplacement-1, 0);

		for (int l = -iDisplacement; l <= iDisplacement; l++) {
			//MATLAB: circshift(dSOSImg(:,idx(i)), iD)
			circshift(aTmp2, 1, 0);

			//MATLAB: (tmp - repmat(circshift(dSOSImg(:,idx(i)), iD),[1 size(tmp,2)]))
			hoNDArray<float> tmp2(tmp.get_dimensions());
			tmp2.fill(0.0); // result of subtraction

			int N = tmp.get_size(0), LE = tmp.get_size(1);
			float *pA = tmp.begin(), *pB = aTmp2.begin(), *pR = tmp2.begin();

			#pragma omp parallel for default(none) schedule(static) shared(N, pA, pB, pR, LE)
			for (int iL = 0; iL < LE; iL++) {
				for (int iE = 0; iE < N; iE++) {
					pR[iE + N*iL] = pA[iE + N*iL] - pB[iE];
				}
			}

			//MATLAB: (tmp - repmat(circshift(dSOSImg(:,idx(i)), iD),[1 size(tmp,2)])).^2
			multiply(tmp2, tmp2, tmp2);

			//MATLAB: sum((tmp - repmat(circshift(dSOSImg(:,idx(i)), iD),[1 size(tmp,2)])).^2,2)
			hoNDArray<float> hTmp2;
			if (tmp2.get_number_of_dimensions() > 1) {
				std::vector<size_t> vD = *tmp2.get_dimensions();
				vD.pop_back();
				hTmp2.create(&vD);
				hTmp2.fill(0.0);

				float *pNewArray = hTmp2.get_data_ptr();
				float *pOldArray = tmp2.get_data_ptr();
				int N = tmp2.get_size(1);
				int L = hTmp2.get_number_of_elements();

				#pragma omp parallel for default(none) schedule(static) shared(N, L, pNewArray, pOldArray)
				for (int sum_dim = 0; sum_dim < N; sum_dim++) {
					for (int i = 0; i < L; i++) {
						pNewArray[i] += pOldArray[i + L*sum_dim];
					}
				}
			}

			//MATLAB: dRMSImg(:,dDisplacement-iD+1) = sum((tmp - repmat(circshift(dSOSImg(:,idx(i)), iD),[1 size(tmp,2)])).^2,2)
			int iOffset = (iDisplacement-l)*hTmp2.get_size(0);
			memcpy(aRMSImg.begin() + iOffset, hTmp2.begin(), sizeof(float)*hTmp2.get_size(0));
		}

		//MATLAB: sum(dRMSImg(dX-round(dDisplacement/2):dX+round(dDisplacement/2),:))
		vStart.clear();
		vStart.push_back(dX-static_cast<int>((static_cast<float>(iDisplacement)/2)+.5));
		vStart.push_back(0);

		vSize.clear();
		vSize.push_back(iDisplacement+1);
		vSize.push_back(aRMSImg.get_size(1));

		get_subarray(aRMSImg, vStart, vSize, aRMSImg);

		hoNDArray<std::complex<float>> cfaRMSImgTest(aRMSImg.get_dimensions());
		cfaRMSImgTest.fill(std::complex<float>(0.0, 0.0));
		std::complex<float> *cfPointer = cfaRMSImgTest.get_data_ptr();
		for (size_t iI = 0; iI < cfaRMSImgTest.get_number_of_elements(); iI++) {
			cfPointer[iI] = std::complex<float>(aRMSImg.at(iI), 0.0);
		}

		sum_dim(aRMSImg, 0, aRMSImg);

		cfaRMSImgTest.clear();
		cfaRMSImgTest.create(aRMSImg.get_dimensions());
		cfaRMSImgTest.fill(std::complex<float>(0.0, 0.0));
		cfPointer = cfaRMSImgTest.get_data_ptr();
		for (size_t iI = 0; iI < cfaRMSImgTest.get_number_of_elements(); iI++) {
			cfPointer[iI] = std::complex<float>(aRMSImg.at(iI), 0.0);
		}

		//MATLAB: min(sum(dRMSImg(dX-round(dDisplacement/2):dX+round(dDisplacement/2),:)))
		int iMinVal = amin(&aRMSImg);

		//MATLAB: dDisplacement + 1 - dNav(i)
		vNav_.at(i) = iDisplacement - iMinVal;

		//MATLAB: circshift(dSOSImg(:,idx(i)), dNav(i))
		hoNDArray<float> aTmp3 = aTmp;
		circshift(aTmp3, vNav_.at(i), 0);

		if (i%20 == 0) {
			GINFO("Getting Navigator - %.1f %%\n", static_cast<float>(aRefImg.get_size(1)-2-i)/static_cast<float>(aRefImg.get_size(1)-2)*100);
		}

		//MATLAB: dRefImg(:,idx(i)) = circshift(dSOSImg(:,idx(i)), dNav(i))
		memcpy(aRefImg.get_data_ptr()+vIdx.at(i)*aRefImg.get_size(0), aTmp3.get_data_ptr(), sizeof(float)*aTmp3.get_size(0));
	}

	for (size_t i = 0; i < vNav_.size(); i++) {
		vNav_.at(i) *= -1;
	}

	// get Gaussian filter kernel and calculate convolution with navigator data
	vGaussian.clear();
	filter1DGaussian(vGaussian, 5);
	vectorConv(vNav_, vGaussian, 0);

	//-------------------------------------------------------------------------
	// interpolate navigator data signal to TR intervals
	GINFO("interpolation of navigator data to TR intervals..\n");

	for (size_t i = 0; i < vNav_.size(); i++) {
		vNav_.at(i) = -vNav_.at(i);
	}

	int iMin = std::min_element(vNav_.begin(), vNav_.end())-vNav_.begin();
	float fMin = vNav_.at(iMin);
	for (size_t i = 0; i < vNav_.size(); i++) {
		vNav_.at(i) -= fMin;
	}

	// build vector with elements 0..lNoScans_ to interpolate vNavInt_ below
	std::vector<float> vNavIndNew;
	for (long i = 0; i < lNoScans_; i++) {
		vNavIndNew.push_back(i);
	}

	GDEBUG("vNavInd size: %i, vNav_ size: %i, vNavIndNew size: %i\n", GlobalVar::instance()->vNavInd_.size(), vNav_.size(), vNavIndNew.size());

	std::vector<float> vNavInd = GlobalVar::instance()->vNavInd_;
	vNavInt_ = interp1<float>(vNavInd, vNav_, vNavIndNew);

	return;
}

// get interpolated navigator signal by Principal Component Analysis
void CS_Retro_NavigatorGadget::getNav2DPCA(hoNDArray<std::complex<float> > &aNav)
{
	GDEBUG("\n\n**************************************\n********** get navigator 2D **********\n**************************************\n\n");

	// reconstruct the 1-D projections for all measurements and all channels
	GINFO("domain transformation - k-space to image\n");

	/*assumptions:
	* iMeasurementTime_ is the total scan time in seconds
	* lNoScans_ is the same as length(iLC(:,15)) in Matlab
	* */

	size_t iNSamples		= aNav.get_size(0);
	size_t iNMeasurement	= aNav.get_size(1);
	size_t iNavRes	 		= aNav.get_size(2);
	size_t iNChannels		= aNav.get_size(3);

	/* MATLAB
	% Reconstruct the 1-D projections for all measurements and all channels
	dImg = fftshift(ifft(ifftshift(dKSpace)));
	*/
	hoNDArray<std::complex<float> > aImg = aNav;
	hoNDFFT_CS<float>::instance()->ifftshift3D(aImg);
	hoNDFFT_CS<float>::instance()->ifft1(aImg);
	hoNDFFT_CS<float>::instance()->fftshift3D(aImg);

	//Concatenate arrays
	std::vector<size_t> A_dims;
	A_dims.push_back(iNSamples*iNavRes*iNChannels);
	A_dims.push_back(iNMeasurement);

	hoNDArray<std::complex<float> > A;
	A.create(&A_dims);

	// check element bounds -> prevent error in for loop
	if (A.get_number_of_elements() != aImg.get_number_of_elements()) {
		GWARN("Oops! Size of A (=%d) != size of aImg (=%d)\n", A.get_number_of_elements(), aImg.get_number_of_elements());

		if (A.get_number_of_elements() > aImg.get_number_of_elements()) {
			GERROR("Segmentation fault prevented! Check sizes of arrays!\n");
			throw std::runtime_error("A.get_number_of_elements() > aImg.get_number_of_elements()");
		}
	}

	for (size_t i = 0; i < A.get_number_of_elements(); i++) {
		A.at(i) = aImg.at(i);
	}

	hoNDKLT<std::complex<float> > *VT = new hoNDKLT <std::complex<float> >;
	std::vector<size_t> coeff_dims;
	coeff_dims.push_back(iNMeasurement);
	coeff_dims.push_back(iNMeasurement);

	hoNDArray<std::complex<float> > coeff;
	coeff.create(&coeff_dims);

	//Compute PCA based on KLT principal components are saved in coeff in descending order
	VT->prepare(A, static_cast<size_t>(1), static_cast<size_t>(0), true);
	VT->eigen_vector(coeff);

	// delete VT
	// TODO: Check why VT is needed - maybe it should be appear somewhere below?
	delete VT;

	double Fs = iNMeasurement/(lNoScans_*(GlobalVar::instance()->fTR_/1000.0));	// Get the sampling frequency (/1000 because fTR_ is in ms, not in s)

	//fft(result, ...,1);
	// fft only 1 dimensional (first dimension)
	hoNDFFT<float>::instance()->fft1(coeff);

	//nfft2 = 2.^nextpow2(nfft);
	int nfft2 = std::pow(2, std::ceil(log(coeff.get_size(0))/log(2)));
	hoNDArray<std::complex<float> > absresult;
	absresult.create(&coeff_dims);

	//fy = fft(y, nfft2);
	//fy = abs(fy);
	for (size_t i = 0; i < coeff.get_number_of_elements(); i++) {
		absresult.at(i) = abs(coeff.at(i));
	}

	int inspectednr = 15; // only search in the first 15 principal components
	std::vector<size_t> fy_dims;
	fy_dims.push_back(inspectednr);
	fy_dims.push_back(nfft2/2-1);
	hoNDArray<std::complex<float> > fy;
	fy.create(&fy_dims);

	//fy = fy(1:nfft2/2);
	for (size_t x = 0; x < fy.get_size(0); x++) {
		for (size_t i = 0; i < fy.get_size(1); i++) {
			fy.at(i+(x*fy.get_size(1))) = absresult.at(i+(x*fy.get_size(1)));
		}
	}

	//define the search area here
	//Fl = nfft2/Fs * 0.66 Hz;
	//Fu = nfft2/Fs * 1.5 Hz;
	float Fl = nfft2/Fs * 0.66;
	float Fu = nfft2/Fs * 1.5;

	//[value, frequency] = max(fy(floor(Fl):floor(Fu),1));
	//coeff of pca are already in a descending order. Searching only the first 15 columns is basically enough and does not introduce errors.
	std::complex<float> maxvalue = 0;
	int frequency = 0;
	int searcharea = std::floor(Fu) - std::floor(Fl);
	size_t column_number = 0;
	for (size_t x = 0; x < fy.get_size(0); x++) {
		for (int i = 0; i < searcharea; i++) {
			size_t pos = i+std::floor(Fl)+(x*fy.get_size(1));

			// break loop if index position is exceeding the limits of fy
			if (pos > (fy.get_number_of_elements() - 1)) {
				GWARN("Maximum index position exceeded!\n");
				break;
			}

			if (compare_complex_values<float>(maxvalue, fy.at(pos)) < 0) {
				maxvalue = fy.at(pos);
				frequency = i+1;
				column_number = x;
			}
		}
	}

	//frequency = ((Fl)+frequency-2)*Fs/nfft2;
	float realfrequency = (Fl + frequency - 2) * Fs / nfft2;

	//dECG = real(coeff(:,coeffnumber)) - imag(coeff(:,coeffnumber));
	// Fs= 1/dTR changing sampling rate because signal is going to be interpolated
	Fs = 1000.0/GlobalVar::instance()->fTR_;
	std::vector<size_t> dECG_dims;
	dECG_dims.push_back(iNMeasurement);
	dECG_dims.push_back(1);

	hoNDArray<std::complex<float> > dECGtemp;
	dECGtemp.create(&dECG_dims);

	for (size_t i = 0; i < iNMeasurement; i++) {
		dECGtemp.at(i) = coeff.at(i+(column_number * iNMeasurement));
	}

	//get the real and the imag part of the signal and subtract them.
	hoNDArray<std::complex<float> > dECGhoNDArray;
	dECGhoNDArray.create(&dECG_dims);
	hoNDArray<std::complex<float> > realpart;
	realpart.create(&dECG_dims);

	// extract real part, for newer compilers also consider:
	//realpart = real(&dECGtemp);
	for (size_t i = 0; i < iNMeasurement; i++) {
		realpart.at(i) = dECGtemp.at(i).real();
	}

	hoNDArray<std::complex<float> > imaginarypart;
	imaginarypart.create(&dECG_dims);

	// extract imaginary part, for newer compilers also consider:
	//imaginarypart = imag(&dECGtemp);
	for (size_t i = 0; i < iNMeasurement; i++) {
		imaginarypart.at(i) = dECGtemp.at(i).imag();
	}

	subtract(&realpart, &imaginarypart, &dECGhoNDArray);

	// type conversion from complex to float and to vector
	std::vector<std::complex<float> > dECG;
	for (size_t i = 0; i < iNMeasurement; i++) {
		dECG.push_back(real(dECGhoNDArray.at(i)));
	}

	//factor = length(iLC)/size(coeff,1);
	//dECG = fScale(dECG , factor);
	std::vector<std::complex<float> > dECGIndtemp;
	for (long i = 0; i <= lNoScans_; i++) {
		dECGIndtemp.push_back(i);
	}

	std::vector<std::complex<float> > dECGInd;

// 	dECGInd = GlobalVar::instance()->vNavInd_;
	for (size_t i = 0; i < GlobalVar::instance()->vNavInd_.size(); i++) {
		dECGInd.push_back(GlobalVar::instance()->vNavInd_.at(i));
	}

	std::vector<std::complex<float> > dECGInt;
	dECGInt = interp1<std::complex<float> >(dECGInd, dECG, dECGIndtemp);

	// Filter the Signal with a first order butterworth filter

	//===============================================================
	//calculate the numerator and denominator of a first order butterworth filter. (End of calulation is indicated by =======)
	//===============================================================

	//ul = 4*tan(pi*fl/2);
	//uh = 4*tan(pi*fh/2);
	//den = [1 0 0];
	float ul = 4*tan(M_PI*(realfrequency-0.1)/(Fs/2)/2);
	float uh = 4*tan(M_PI*(realfrequency+0.1)/(Fs/2)/2);
	std::vector<std::complex<float> > den;
	den.push_back(1);
	den.push_back(0);
	den.push_back(0);

	//Bandwidth and center frequency
	float Bw = uh - ul;
	float Wn = std::sqrt(ul*uh);

	//t1 = [1+(Wn*(-Bw/Wn)/4) Wn/4; -Wn/4 1];
	//t2 = [1-(Wn*(-Bw/Wn)/4) -Wn/4; Wn/4 1];
	//ad = inv(t2)*t1;

	std::vector<size_t> t1_dims;
	t1_dims.push_back(2);
	t1_dims.push_back(2);
	hoNDArray<std::complex<float> > t1;
	t1.create(&t1_dims);
	hoNDArray<std::complex<float> > t2;
	t2.create(&t1_dims);
	hoNDArray<std::complex<float> > ad;
	ad.create(&t1_dims);

	t1.at(0) = 1+(Wn*(-Bw/Wn)/4);
	t1.at(1) = Wn/4;
	t1.at(2) = -Wn/4;
	t1.at(3) = 1;

	// also transpose
	t2.at(0)= 1-(Wn*(-Bw/Wn)/4);
	t2.at(1)= -Wn/4;
	t2.at(2)= Wn/4;
	t2.at(3)= 1;

	// matlab: ad = inv(t2)*t1;
	std::complex<float> det_t2 = t2.at(0)*t2.at(3)-t2.at(1)*t2.at(2);
	ad.at(0) = (+t2.at(3)*t1.at(0)-t2.at(1)*t1.at(2))/det_t2;
	ad.at(1) = (+t2.at(3)*t1.at(1)-t2.at(1)*t1.at(3))/det_t2;
	ad.at(2) = (-t2.at(2)*t1.at(0)+t2.at(0)*t1.at(2))/det_t2;
	ad.at(3) = (-t2.at(2)*t1.at(1)+t2.at(0)*t1.at(3))/det_t2;

	//%den = poly(ad);
	//e = eig(ad);
	std::vector<size_t> e_dims;
	e_dims.push_back(2);
	e_dims.push_back(1);
	hoNDArray<std::complex<float> > e;
	e.create(&e_dims);

	std::vector<size_t> kern_dims;
	kern_dims.push_back(3);
	kern_dims.push_back(1);
	hoNDArray<std::complex<float> > kern;
	kern.create(&kern_dims);
	std::vector<std::complex<float> > num;
	num.push_back(0);
	num.push_back(0);
	num.push_back(0);

	//eig
	hoNDKLT<std::complex<float> > *Vtemp = new hoNDKLT<std::complex<float> >;

	//get the eigen value
	Vtemp->prepare(ad, static_cast<size_t>(1), static_cast<size_t>(0), false);
	Vtemp->eigen_value(e);

	// delete Vtemp
	// TODO: Check why Vtemp is needed - maybe it should be appear somewhere below?
	delete Vtemp;

	//% Expand recursion formula
	//den(2) = den(2) - e(1)*den(1);
	//den(3) = den(3) - e(2)*den(2);
	//den(2) = den(2) - e(2)*den(1);
	den.at(1) = den.at(1) - e.at(0) * den.at(0);
	den.at(2) = den.at(2) - e.at(1) * den.at(1);
	den.at(1) = den.at(1) - e.at(1) * den.at(0);

	Wn = 2*atan(Wn/4);

	//%  normalize so |H(w)| == 1:
	//%kern = exp(-1i*Wn*(0:2));
	for (int k = 0; k < 3; k++) {
		kern.at(k) = std::exp(std::complex<float>(0.0, -1.0)*static_cast<std::complex<float> >(Wn)*static_cast<std::complex<float> >(k));
	}

	//f = (kern(1)*den(1)+kern(2)*den(2)+kern(3)*den(3))/(kern(1)-kern(3));
	std::complex<float> f = (kern.at(0)*den.at(0)+kern.at(1)*den.at(1)+kern.at(2)*den.at(2))/(kern.at(0)-kern.at(2));
	num.at(0) = f.real();
	num.at(1) = 0;
	num.at(2) = - f.real();

	//===========================================================
	//end of calculating the numerator and denominator of the first order butterworth filter
	//===========================================================

	//filtfilt() equivalent function. b = num and a = den


	//============================================================
	//start of zero phase digital filter function
	//============================================================

	//n    = length(den); always 3 in first order case
	//z(n) = 0;
	//num = num / den(1);
	//den = den / den(1);
	std::vector<std::complex<float> > z;
	z.push_back(0);
	z.push_back(0);
	z.push_back(0);
	num.at(0) = num.at(0)/den.at(0); //num.at(1) is always 0
	num.at(2) = num.at(2)/den.at(0);

	den.at(0) = den.at(0)/den.at(0);
	den.at(1) = den.at(1)/den.at(0);
	den.at(2) = den.at(2)/den.at(0);

	//Y    = zeros(size(X));
	//for m = 1:length(Y)
	//  Y(m) = num(1) * X(m) + z(1);
	//   for i = 2:n
	//      z(i - 1) = num(i) * X(m) + z(i) - den(i) * Y(m);
	//   end
	//end
	std::vector<std::complex<float> > Y;
	for (size_t m = 0; m < dECG.size(); m++) {
		Y.push_back(num.at(0) * dECG.at(m) + z.at(0));

		for (size_t i = 1; i < den.size(); i++) {
			z.at(i - 1) = num.at(i) * dECGInt.at(m) + z.at(i) - den.at(i) * Y.at(m);
		}
	}

	//clear z
	//z(n) = 0;
	z.at(0) = 0;
	z.at(1) = 0;
	z.at(2) = 0;

	//flip vector
	std::reverse(Y.begin(),Y.end());
	dECG.clear();
	dECG = Y;

	//Y    = zeros(size(X));
	// second round filtering (backward)
	//for m = 1:length(Y)
	//   Y(m) = b(1) * X(m) + z(1);
	//   for i = 2:n
	//      z(i - 1) = b(i) * X(m) + z(i) - a(i) * Y(m);
	//   end
	//end
	Y.clear();
	for (size_t m = 0; m < dECG.size(); m++) {
		Y.push_back(num.at(0) * dECG.at(m) + z.at(0));

		for (size_t i = 1; i < den.size(); i++) {
			z.at(i - 1) = num.at(i) * dECGInt.at(m) + z.at(i) - den.at(i) * Y.at(m);
		}
	}

	//flip again
	std::reverse(Y.begin(),Y.end());
	dECG = Y;
	Y.clear();

	//============================================================
	//end of zero phase digital filter function
	//============================================================

	//dECG = diff(dECG);
	std::vector<std::complex<float> > dECGdiff;
	for (size_t i = 0; i < dECG.size()-1; i++) {
		dECGdiff.push_back(dECG.at(i+1)-dECG.at(i));
	}

	return;
}

GADGET_FACTORY_DECLARE(CS_Retro_NavigatorGadget)
