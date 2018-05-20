#include "CS_Retro_NavigatorGadget.h"

#include <armadillo>
#include <algorithm>
#include <cmath>

#include "hoNDArray_math_util.cpp"

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
	// get gadget property
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	min_card_freq_	= MinCardFreq.value();
	max_card_freq_	= MaxCardFreq.value();
	min_resp_freq_	= MinRespFreq.value();
	max_resp_freq_	= MaxRespFreq.value();
	iNavMethod_		= NavigationMethod.value();
#else
	min_card_freq_	= *(get_int_value("MinCardFreq").get());
	max_card_freq_	= *(get_int_value("MaxCardFreq").get());
	min_resp_freq_	= *(get_int_value("MinRespFreq").get());
	max_resp_freq_	= *(get_int_value("MaxRespFreq").get());
	iNavMethod_		= *(get_int_value("NavigationMethod").get());
#endif

	// some basic error checking for parameters
	if (min_card_freq_ < 0 || max_card_freq_ < 0 || min_resp_freq_ < 0 || max_resp_freq_ < 0) {
		GERROR("Given parameters min_card_freq_, max_card_freq_, min_resp_freq_, max_resp_freq_ must not be negative!\n");
		return GADGET_FAIL;
	}

	if (min_card_freq_ >= max_card_freq_) {
		GERROR("max_card_freq_ must be greater than min_card_freq_!\n");
		return GADGET_FAIL;
	}

	if (min_resp_freq_ >= max_resp_freq_) {
		GERROR("max_resp_freq_ must be greater than min_resp_freq_!\n");
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

int CS_Retro_NavigatorGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1,GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m3)
{
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

	// free memory
	m2->cont(NULL);
	m2->release();

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

	// ATTENTION: Work with [s] values
	// dNavPeriod = dNavPeriod/1000; dTR = dTR/1000;
	min_card_freq_ /= 60;
	max_card_freq_ /= 60;
	min_resp_freq_ /= 60;
	max_resp_freq_ /= 60;

	// 1. Step: Restack
	// is already done in C++ code

	// 2. Step: Compute PCA based on KLT principal components are saved in coeff in descending order
	std::vector<size_t> coeff_dims;
	coeff_dims.push_back(iNMeasurement);
	coeff_dims.push_back(iNMeasurement);

	hoNDArray<std::complex<float> > coeff;
	coeff.create(&coeff_dims);
	coeff.delete_data_on_destruct(true);

	// prepare aImg for KLT
	// first: permute
	std::vector<size_t> aImg_new_order;
	aImg_new_order.push_back(0);
	aImg_new_order.push_back(2);
	aImg_new_order.push_back(3);
	aImg_new_order.push_back(1);
	aImg = *permute(&aImg, &aImg_new_order,false);

	// then reshape
	std::vector<size_t> new_aImg_dims;
	new_aImg_dims.push_back(iNSamples*iNavRes*iNChannels);
	new_aImg_dims.push_back(iNMeasurement);
	aImg.reshape(&new_aImg_dims);

	GINFO("Performing KLT (may take a while)\n");

	hoNDKLT<std::complex<float> > VT;
	VT.prepare(aImg, static_cast<size_t>(1), static_cast<size_t>(0), true);
	VT.eigen_vector(coeff);

	GINFO("Continuing with search for motion\n");

	// 3. Step: search for respiratory motion
	// calculate frequency
	//%dFs = 1/((length(iLC)*dTR/1000)/(length(coeff)));
	double f_s = iNMeasurement/(lNoScans_*(GlobalVar::instance()->fTR_/1000.0)); // Get the sampling frequency (/1000 because fTR_ is in ms, not in s)

	// get next base of 2
	//%dFactorfft = 2.^nextpow2(size(dCoeff,1));
	unsigned int factor_fft = std::pow(2, std::ceil(log(coeff.get_size(0))/log(2)));

	// first zero padd coeff
	hoNDArray<std::complex<float> > coeff_padded;
	std::vector<size_t> coeff_padded_dims;
	coeff_padded_dims.push_back(factor_fft);
	coeff_padded_dims.push_back(coeff.get_size(1));
	coeff_padded.create(&coeff_padded_dims);
	for (size_t col = 0; col < coeff.get_size(1); col++) {
		size_t offset_old = col * coeff.get_size(0);
		size_t offset_new = col * coeff_padded.get_size(0);
		memcpy(coeff_padded.get_data_ptr()+offset_new, coeff.get_data_ptr()+offset_old, coeff.get_size(0)*sizeof(coeff.at(0)));
	}

	// and continue with fft
	// note: ifftshift2D is done before zero padding, thus none here
	hoNDFFT_CS<float>::instance()->fft1(coeff_padded);

	//%dCoeffF = abs(dCoeffF);
	hoNDArray<float> coeff_abs;
	coeff_abs.create(&coeff_padded_dims);
	for (size_t i = 0; i < coeff_padded.get_number_of_elements(); i++) {
		coeff_abs.at(i) = abs(coeff_padded.at(i));
	}

	// filtering
	// calculate lower and upper boundary frequency
	// Note: we only need the floor() value, so we build it directly instead of later on.
	//%Fl = dFactorfft/dFs * dCutOffResp(1);
	const unsigned int f_l = std::floor(factor_fft/f_s * min_resp_freq_);
	//%Fu = dFactorfft/dFs * dCutOffResp(2);
	const unsigned int f_u = std::floor(factor_fft/f_s * max_resp_freq_);

	if (f_u <= f_l) {
		GERROR("f_u (=%d) must be greater than f_l (=%d). Please set parameters correct!\n", f_u, f_l);
		throw runtime_error("Illegal parameters min_resp_freq_ or max_resp_freq_!");
	}

	// now get iPeak (dVal can be omitted). iPeak is the column number where the maximum value occures.
	// we could either implement it the Matlab way (crop hoNDArray, search max value per line and max col position)
	// or we do some intelligent data handling to get the position directly:
	float max_val = std::numeric_limits<float>::min();		// initialize as minimum (it could also be 0 because coeff_abs only contains abs values)
	size_t peak_position = 0;		// note: -1 is insufficient because size_t is unsigned, so pos is max. But we will rewrite max_val either.
	for (size_t f = static_cast<size_t>(f_l); f <= static_cast<size_t>(f_u); f++) {
		for (size_t component_number = search_range_min_-1; component_number < search_range_max_; component_number++) {		// be aware: search_range counts MATLAB like (from 1 to length)
			size_t pos = component_number * coeff_abs.get_size(0) + f;

			if (max_val < coeff_abs.at(pos)) {
				max_val = coeff_abs.at(pos);
				peak_position = component_number;
			}
		}
	}

	// now extract navi signal
	//%dRespNavi = real(dCoeff(:,iPeak)) - imag(dCoeff(:,iPeak));
	std::vector<float> resp_navi;

	for (size_t f = 0; f < coeff.get_size(0); f++) {
		size_t pos = peak_position * coeff.get_size(0) + f;
		resp_navi.push_back(coeff.at(pos).real() - coeff.at(pos).imag());
	}

	// and convolute
	//%dRespNavi = conv(dRespNavi, fGaussianLP(5), 'same');
	// first build gaussian low pass
	std::vector<float> gaussian_lowpass = {
		0.0269,
		0.2334,
		0.4794,
		0.2334,
		0.0269,
	};
	// convolute
	resp_navi = arma::conv_to<std::vector<float> >::from(arma::conv(arma::Col<float>(resp_navi), arma::Col<float>(gaussian_lowpass), "same"));

	// Note: Matlab implementation: %dRespNavi = -dRespNavi;
	// is implicitly done during KLT and can be omitted here

	// build vector with elements 0..lNoScans_ to interpolate vNavInt_ below
	std::vector<float> nav_ind_new;
	for (long i = 0; i < lNoScans_; i++) {
		nav_ind_new.push_back(i);
	}

	// get navigator indices
	std::vector<float> nav_ind = GlobalVar::instance()->vNavInd_;

	GDEBUG("nav_ind size: %i, resp_navi size: %i, nav_ind_new size: %i\n", nav_ind.size(), resp_navi.size(), nav_ind_new.size());

	// interpolate to output vector
	vNavInt_ = interp1<float>(nav_ind, resp_navi, nav_ind_new);

	return;
}

/*
 * WARNING: This code has not been tested yet! You will have some fun trying to get it to work.
 * For more implementation details, you can look in the GIT history: CS_Retro_PCANavigatorGadget is you friend.
 * But beware: The code there is erroneous and has already been corrected (watch GIT history if you want to avoid doubled fun)
 * You may also refer to the Matlab implementation.
 */
void CS_Retro_NavigatorGadget::butterworth_filtering(const double fl, const double fh, std::vector<float> &signal)
{
	// Filter the Signal with a first order butterworth filter

	//===============================================================
	//calculate the numerator and denominator of a first order butterworth filter. (End of calulation is indicated by =======)
	//===============================================================

	//ul = 4*tan(pi*fl/2);
	//uh = 4*tan(pi*fh/2);
	float ul = 4*tan(M_PI*fl/2);
	float uh = 4*tan(M_PI*fh/2);

	//Bandwidth and center frequency
	float Bw = uh - ul;
	float Wn = std::sqrt(ul*uh);

	// Matlab:
	//t1 = [1+(Wn*(-Bw/Wn)/4) Wn/4; -Wn/4 1];
	//t2 = [1-(Wn*(-Bw/Wn)/4) -Wn/4; Wn/4 1];
	//ad = inv(t2)*t1;
	// Note: indexing follows mathematical convention: t1[row number][column number]
	float t1[2][2], t2[2][2];
	std::complex<float> ad[2][2];	// note: only real values in ad, but calculation of e later on must be complex

	t1[0][0] = 1+(Wn*(-Bw/Wn)/4);
	t1[0][1] = Wn/4;
	t1[1][0] = -Wn/4;
	t1[1][1] = 1;

	// also transpose
	t2[0][0]= 1-(Wn*(-Bw/Wn)/4);
	t2[0][1]= -Wn/4;
	t2[1][0]= Wn/4;
	t2[1][1]= 1;

	// matlab: ad = inv(t2)*t1;
	float det_t2 = t2[0][0]*t2[1][1]-t2[0][1]*t2[1][0];
	ad[0][0] = (+t2[1][1]*t1[0][0]-t2[0][1]*t1[1][0])/det_t2;
	ad[0][1] = (+t2[1][1]*t1[0][1]-t2[0][1]*t1[1][1])/det_t2;
	ad[1][0] = (-t2[1][0]*t1[0][0]+t2[0][0]*t1[1][0])/det_t2;
	ad[1][1] = (-t2[1][0]*t1[0][1]+t2[0][0]*t1[1][1])/det_t2;

	//%den = poly(ad);
	//e = eig(ad);
	std::complex<float> e[2];

	// computate eigenvalues
	e[0] = (ad[0][0]+ad[1][1]+std::sqrt(std::pow(ad[0][0]+ad[1][1], 2) - std::complex<float>(4)*(ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0])))/std::complex<float>(2);
	e[1] = (ad[0][0]+ad[1][1]-std::sqrt(std::pow(ad[0][0]+ad[1][1], 2) - std::complex<float>(4)*(ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0])))/std::complex<float>(2);

	std::complex<float> kern[3];

	//den = [1 0 0];
	float den[] = {1, 0, 0};

	// In Matlab:
	//% Expand recursion formula
	//den(2) = den(2) - e(1)*den(1);
	//den(3) = den(3) - e(2)*den(2);
	//den(2) = den(2) - e(2)*den(1);
	//
	// Some thoughts:
	//	- e is either real or conjugate complex (see above)
	//	- in complex case: calculation of den becomes real (because e conjugate complex)
	//	- then: imaginary part can always be ignored and the matlab calculation reduces to:
	den[1] = -e[0].real() - e[1].real();
	den[2] = e[0].real() * e[1].real() + std::pow(e[0].imag(), 2);

	Wn = 2*atan(Wn/4);

	//%  normalize so |H(w)| == 1:
	//%kern = exp(-1i*Wn*(0:2));
	for (size_t k = 0; k < ARRAYSIZE(kern); k++) {
		kern[k] = std::exp(std::complex<float>(0.0, -1.0)*std::complex<float>(Wn)*std::complex<float>(k));
	}

	//f = (kern(1)*den(1)+kern(2)*den(2)+kern(3)*den(3))/(kern(1)-kern(3));
	std::complex<float> f = (kern[0]*den[0]+kern[1]*den[1]+kern[2]*den[2])/(kern[0]-kern[2]);

	float num[3] = {f.real(), 0, -f.real()};

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
	float z[3] = {0};

	//Y    = zeros(size(X));
	//for m = 1:length(Y)
	//  Y(m) = num(1) * X(m) + z(1);
	//   for i = 2:n
	//      z(i - 1) = num(i) * X(m) + z(i) - den(i) * Y(m);
	//   end
	//end
	std::vector<float> Y;
	for (size_t m = 0; m < signal.size(); m++) {
		Y.push_back(num[0] * signal.at(m) + z[0]);

		for (size_t i = 1; i < ARRAYSIZE(den); i++) {
			z[i-1] = num[i] * signal.at(m) + z[i] - den[i] * Y.at(m);
		}
	}

	//clear z
	//z(n) = 0;
	z[0] = 0;
	z[1] = 0;
	z[2] = 0;

	//flip vector
	std::reverse(Y.begin(),Y.end());
	signal.clear();
	signal = Y;

	//Y    = zeros(size(X));
	// second round filtering (backward)
	//for m = 1:length(Y)
	//   Y(m) = b(1) * X(m) + z(1);
	//   for i = 2:n
	//      z(i - 1) = b(i) * X(m) + z(i) - a(i) * Y(m);
	//   end
	//end
	Y.clear();
	for (size_t m = 0; m < signal.size(); m++) {
		Y.push_back(num[0] * signal.at(m) + z[0]);

		for (size_t i = 1; i < ARRAYSIZE(den); i++) {
			z[i-1] = num[i] * signal.at(m) + z[i] - den[i] * Y.at(m);
		}
	}

	//flip again
	std::reverse(Y.begin(),Y.end());
	signal = Y;
	Y.clear();

	//============================================================
	//end of zero phase digital filter function
	//============================================================
}

GADGET_FACTORY_DECLARE(CS_Retro_NavigatorGadget)
