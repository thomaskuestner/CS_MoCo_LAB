#include "CS_Retro_NavigatorGadget.h"

namespace Gadgetron{
// class constructor
CS_Retro_NavigatorGadget::CS_Retro_NavigatorGadget() {}
 
// class destructor - delete temporal buffer/memory
CS_Retro_NavigatorGadget::~CS_Retro_NavigatorGadget(){
	
}

// read flexible data header
int CS_Retro_NavigatorGadget::process_config(ACE_Message_Block* mb)
{
	
	return GADGET_OK;
}

int CS_Retro_NavigatorGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, GadgetContainerMessage< hoNDArray <std::complex<float> > >* m3){

	iNoChannels_ = m1->getObjectPtr()->channels;
	iNoNav_		 = m1->getObjectPtr()->user_int[5];
	lNoScans_	 = m3->getObjectPtr()->get_size(1);

	if (getNav2D(*m2->getObjectPtr())){
		if (bMatlab_){
			//mexPrintf("Error in getNav2D\n");mexEvalString("drawnow;");
		}
		else{
			GADGET_DEBUG1("Error in getNav2D\n");
		}
	}

	GadgetContainerMessage< hoNDArray<float> >* tmp_m2 = new GadgetContainerMessage< hoNDArray<float> >();

	tmp_m2->getObjectPtr()->create(vNavInt_.size());

	// convert vector to array
	float* fPtr = tmp_m2->getObjectPtr()->get_data_ptr();
	for (long iI = 0; iI < vNavInt_.size(); iI++){
		fPtr[iI] = vNavInt_.at(iI);
	}

	//// debug output
	//if (bMatlab_){
	//	//mexPrintf("\nNavigator signal\n");mexEvalString("drawnow;");
	//}
	//else{
	//	GADGET_DEBUG1("\nNavigator signal\n");
	//	tmp_m2->getObjectPtr()->print(std::cout);
	//}

	m1->cont(tmp_m2);
	tmp_m2->cont(m3);

	if (this->next()->putq(m1) < 0) {
    	return GADGET_FAIL;
	}

	return GADGET_OK;
}

// get interpolated navigator signal
bool CS_Retro_NavigatorGadget::getNav2D(hoNDArray<std::complex<float>> &aNav){
	
	if (bMatlab_){
		//mexPrintf("\n\n**************************************\n********** get navigator 2D **********\n**************************************\n\n");mexEvalString("drawnow;");
		GlobalVar::instance()->vNavInd_ = vNavInd_;
	}
	else{
		GADGET_DEBUG1("\n\n**************************************\n********** get navigator 2D **********\n**************************************\n\n");
	}	
	
	// reconstruct the 1-D projections for all measurements and all channels
	if (bMatlab_){
		//mexPrintf("domain transformation - k-space to image\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("domain transformation - k-space to image\n");
	}

	/* MATLAB
	% Reconstruct the 1-D projections for all measurements and all channels
	dImg = fftshift(ifft(ifftshift(dKSpace)));
	dImg = fftshift(ifft(ifftshift(dImg, 3), [], 3), 3);
	dImg = flipdim(dImg, 1); % Invert the RO direction: 1-N -> H-F
	dImg = dImg(iNSamples/4:iNSamples.*3/4 - 1, :, :, :); % RO x t x PE x CH
	*/
	hoNDArray<std::complex<float>> aImg = aNav;
	hoNDFFT_CS<float>::instance()->ifftshift3D(aImg);
	hoNDFFT_CS<float>::instance()->ifft1(aImg);
	hoNDFFT_CS<float>::instance()->fftshift3D(aImg);
	hoNDFFT_CS<float>::instance()->ifft(&aImg, 2, true);
	flip_array(aImg, 0);

	// crop center part due to twofold oversampling
	if (bMatlab_){
		//mexPrintf("crop center part of two-fold oversampled data..\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("crop center part of two-fold oversampled data..\n");
	}

	std::vector<size_t> vStart, vSize;
	vStart.push_back((size_t)(std::floor(0.25*aNav.get_size(0)))-1);  vStart.push_back(0); vStart.push_back(0); vStart.push_back(0);
	vSize.push_back(std::floor(0.5*aNav.get_size(0))); vSize.push_back(aNav.get_size(1)); vSize.push_back(aNav.get_size(2)); vSize.push_back(aNav.get_size(3)); 
	get_subarray(aImg, vStart, vSize, aImg);

	// get channel power and normalize channel images
	if (bMatlab_){
		//mexPrintf("calculate channel power and normalize channel images..\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("calculate channel power and normalize channel images..\n");
	}

	std::vector<float> fPower(iNoChannels_);
	hoNDArray<std::complex<float>> aPower = aImg; aPower.fill(std::complex<float>(0,0));
	for (int c = 0; c < iNoChannels_; c++){
		size_t offset = aImg.get_size(0)*aImg.get_size(1)*aImg.get_size(2)*c;
		hoNDArray<std::complex<float>> SubArray(aImg.get_size(0), aImg.get_size(1), aImg.get_size(2), aImg.get_data_ptr()+offset, false);
		fPower.at(c) = asum(&SubArray);
		// fill part of the 3D array					
		#pragma  omp parallel for
		for (long i = 0; i < aImg.get_size(0)*aImg.get_size(1)*aImg.get_size(2); i++)
			aPower.at(i+offset) = (fPower.at(c), fPower.at(c));
	}
	divide(aImg, aPower, aImg);
	if (bMatlab_){
		//mexPrintf("channel images normalized..\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("channel images normalized..\n");
	}
	aPower.clear();

	// get x range of respiratory motion & FFT without scrambling
	hoNDArray<std::complex<float>> aFreq = aImg;
	hoNDFFT_CS<float>::instance()->fft(&aFreq, 1, false);
	multiplyConj(aFreq,aFreq,aPower);
	
	// conversion from complex float to float
	hoNDArray<float> afPower(aPower.get_dimensions());
	for (long i = 0; i < afPower.get_number_of_elements(); i++){
		afPower[i] = aPower[i].real();
	}

	/*
	dIMGres = 1./(double(dNavPeriod)./1000.*double(iNMeasurements)); % The frequency resolution of dIMG in Hz
	dPower = squeeze(sum(dPower(:, round(1./(5.*dIMGres)):round(1./(3.*dIMGres)), :, :), 2)); % RO x PE x CH
	*/
	float fIMGRes = 1.0/(((float)iNavPeriod_/(float)1000)*(float)iNoNav_); // frequency resolution of aImg in Hz
	hoNDArray<float> aPowerInChan, aPowerAcrossChan;
	vStart.clear(); vStart.push_back(0); vStart.push_back(std::floor(1.0/(5*fIMGRes)-.5));vStart.push_back(0);vStart.push_back(0);
	vSize.clear(); vSize.push_back(afPower.get_size(0)); vSize.push_back(std::ceil(1.0/(3*fIMGRes)-.5)-std::ceil(1.0/(5*fIMGRes)-.5)+1); vSize.push_back(afPower.get_size(2)); vSize.push_back(afPower.get_size(3));
	get_subarray(afPower, vStart, vSize, afPower);				
	sum_dim(afPower, 1, afPower); // RO x PE x CH
	sum_dim(afPower, 1, aPowerInChan); // RO x CH
	sum_dim(aPowerInChan, 1, aPowerAcrossChan); // RO x 1

	// Prevent detection of regions in the abdomen
	for (long i = aPowerAcrossChan.get_number_of_elements(); i > std::floor(aPowerAcrossChan.get_number_of_elements()*.75); i--)
		aPowerAcrossChan.at(i) = 0;

	// get Gaussian filter kernel
	std::vector<float> vGaussian;
	filter1DGaussian(vGaussian, 20);
	arrayConv(aPowerAcrossChan, vGaussian, 0);	

	// find index of maximum
	int iMaxIndex = amax(&aPowerAcrossChan);
	if (bMatlab_){
		//mexPrintf("data filtered and maximum determined.. iMaxIndex: %i\n", iMaxIndex);mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG2("data filtered and maximum determined.. iMaxIndex: %i\n", iMaxIndex);
	}
	if ( (iMaxIndex<20) || (iMaxIndex > aPowerInChan.get_size(1)-20)){		
		if (bMatlab_){
			//mexPrintf("error occured in CS_Retro - iMaxIndex out of bounds..\n");mexEvalString("drawnow;");
		}
		else{
			GADGET_DEBUG1("error occured in CS_Retro - iMaxIndex out of bounds..\n");
		}
		return GADGET_FAIL;
	}

	//-------------------------------------------------------------------------
	// sort out channels with no relevant information in target area		
	if (bMatlab_){
		//mexPrintf("get channels which contain most information..\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("get channels which contain most information..\n");
	}
	vStart.clear(); vStart.push_back(iMaxIndex-20); vStart.push_back(0);
	vSize.clear(); vSize.push_back(41); vSize.push_back(iNoChannels_);
	get_subarray(aPowerInChan, vStart, vSize,aPowerInChan);
	if (bMatlab_){
		//mexPrintf("41 elements around maximum cropped..\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("41 elements around maximum cropped..\n");
	}

	std::vector<float> vGoodChannels;
	hoNDArray<float> aTmp;
	// loop over channels
	for (int c = 0; c < iNoChannels_; c++){
		size_t offset = aPowerInChan.get_size(0)*c;
		aTmp.create(aPowerInChan.get_size(0), aPowerInChan.get_data_ptr()+offset, false);
		int iTmpInd = amax(&aTmp);
		vGoodChannels.push_back(aTmp.at(iTmpInd));
	}
	int iIndex = std::max_element(vGoodChannels.begin(), vGoodChannels.end())- vGoodChannels.begin();
	float fMax = vGoodChannels.at(iIndex);
	int iNumGood = 0;
	for (int i = 0; i < vGoodChannels.size(); i++)
		if (vGoodChannels.at(i) > .2*fMax){
			vGoodChannels.at(i) = 1;
			iNumGood++;
		}
		else
			vGoodChannels.at(i) = 0;
	
	// debug output
	for (size_t i = 0; i < vGoodChannels.size(); i++){
		if (bMatlab_){
			//mexPrintf("vGoodChannels[%i]: %f\n", i,vGoodChannels.at(i));mexEvalString("drawnow;");
		}
		else{
			GADGET_DEBUG2("vGoodChannels[%i]: %f\n", i,vGoodChannels.at(i));
		}
	}

	//-------------------------------------------------------------------------
	// get the best PE line
	// get sub array of good channels	
	if (bMatlab_){
		//mexPrintf("get best PE line..\n afPower:\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("get best PE line..\n afPower:\n");
	}
	hoNDArray<float> aPowerInPE(41,afPower.get_size(1), iNumGood); aTmp.clear();
	vStart.clear(); vStart.push_back(iMaxIndex-20); vStart.push_back(0); vStart.push_back(0);
	vSize.clear(); vSize.push_back(41); vSize.push_back(afPower.get_size(1)); vSize.push_back(1);
	if (bMatlab_){
		//mexPrintf("vStart: %i, %i, %i, vSize: %i, %i, %i, afPower size: %i, %i, %i\n",  vStart.at(0), vStart.at(1), vStart.at(2), vSize.at(0), vSize.at(1), vSize.at(2), afPower.get_size(0), afPower.get_size(1), afPower.get_size(2));mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG2("vStart: %i, %i, %i, vSize: %i, %i, %i, afPower size: %i, %i, %i\n",  vStart.at(0), vStart.at(1), vStart.at(2), vSize.at(0), vSize.at(1), vSize.at(2), afPower.get_size(0), afPower.get_size(1), afPower.get_size(2));
	}
	size_t o = 0; //helper
	for (int c = 0; c < iNoChannels_; c++){
		if (vGoodChannels.at(c) == 1){
			vStart.at(2) = c;
			get_subarray(afPower, vStart, vSize, aTmp);
			// fill array
			size_t offset = aPowerInPE.get_size(0)*aPowerInPE.get_size(1)*o;
			for (long i = 0; i < aTmp.get_number_of_elements(); i++){
				aPowerInPE.at(i + offset) = aTmp.at(i);
			}
			o++;
		}
	}

	sum_dim(aPowerInPE, 2, aPowerInPE);
	sum_dim(aPowerInPE, 0, aPowerInPE);
	if (bMatlab_){
		//mexPrintf("\n aPowerInPE\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("\n aPowerInPE\n");
		aPowerInPE.print(std::cout);
	}
	iMaxIndex = amax(&aPowerInPE);
	int iMaxChan = iMaxIndex;	
	if (bMatlab_){
		//mexPrintf("found at %i\n", iMaxIndex);mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG2("found at %i\n", iMaxIndex);
	}	

	//-------------------------------------------------------------------------
	// find best corresponding channels according to best phase encoding position
	if (bMatlab_){
		//mexPrintf("searching for channels according to best PE line..\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("searching for channels according to best PE line..\n");
	}
	// get sub array
	vStart.clear(); vStart.push_back(0); vStart.push_back(0); vStart.push_back(iMaxIndex); vStart.push_back(0);
	vSize.clear(); vSize.push_back(aImg.get_size(0)); vSize.push_back(aImg.get_size(1)); vSize.push_back(1); vSize.push_back(aImg.get_size(3));	
	aFreq.clear();
	get_subarray(aImg, vStart, vSize, aFreq);

	hoNDFFT_CS<float>::instance()->fft(&aFreq, 1, false);
	multiplyConj(aFreq,aFreq,aPower);
	
	// conversion from complex float to float
	afPower.clear(); afPower.create(*aPower.get_dimensions());
	for (long i = 0; i < afPower.get_number_of_elements(); i++){
		afPower[i] = aPower[i].real();
	}

	fIMGRes = 1.0/(((float)iNavPeriod_/(float)1000)*(float)iNoNav_); // frequency resolution of aImg in Hz
	aPowerInChan.clear(); aPowerAcrossChan.clear();
	vStart.clear(); vStart.push_back(0); vStart.push_back(std::ceil(1.0/(5*fIMGRes)-.5));vStart.push_back(0);vStart.push_back(0);
	vSize.clear(); vSize.push_back(afPower.get_size(0)); vSize.push_back(std::ceil(1.0/(3*fIMGRes)-.5)-std::ceil(1.0/(5*fIMGRes)-.5)+1); vSize.push_back(afPower.get_size(2)); vSize.push_back(afPower.get_size(3));
	if (bMatlab_){
//		mexPrintf("vSize: %i, %i, %i, %i - fIMGRes: %f\n", vSize.at(0), vSize.at(1), vSize.at(2), vSize.at(3), fIMGRes);mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG2("vSize: %i, %i, %i, %i - fIMGRes: %f\n", vSize.at(0), vSize.at(1), vSize.at(2), vSize.at(3), fIMGRes);	
	}
	get_subarray(afPower, vStart, vSize, afPower);	
	sum_dim(afPower, 1, afPower); // RO x PE x CH
	sum_dim(afPower, 1, aPowerInChan); // RO x CH
	sum_dim(aPowerInChan, 1, aPowerAcrossChan); // RO x 1

	// Prevent detection of regions in the abdomen
	for (long i = aPowerAcrossChan.get_number_of_elements(); i > std::floor(aPowerAcrossChan.get_number_of_elements()*.75); i--)
		aPowerAcrossChan.at(i) = 0;
	
	if (bMatlab_){
//		mexPrintf("filter data with Gaussian kernel..\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("filter data with Gaussian kernel..\n");	
	}
	// get Gaussian filter kernel
	vGaussian.clear();
	filter1DGaussian(vGaussian, 20);
	arrayConv(aPowerAcrossChan, vGaussian);

	// find index of maximum
	iMaxIndex = amax(&aPowerAcrossChan);
	int dX = iMaxIndex;
	
	if (bMatlab_){
//		mexPrintf("found at %i\n", iMaxIndex);mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG2("found at %i\n", iMaxIndex);
	}

	// get good channels
	vStart.clear(); vStart.push_back(iMaxIndex-20); vStart.push_back(0);
	vSize.clear(); vSize.push_back(41); vSize.push_back(aPowerInChan.get_size(1));
	get_subarray(aPowerInChan, vStart, vSize,aPowerInChan);
	vGoodChannels.clear();aTmp.clear();
	// loop over channels
	for (int c = 0; c < iNoChannels_; c++){
		size_t offset = aPowerInChan.get_size(0)*c;
		aTmp.create(aPowerInChan.get_size(0), aPowerInChan.get_data_ptr()+offset, false);
		int iTmpInd = amax(&aTmp);
		vGoodChannels.push_back(aTmp.at(iTmpInd));
	}
	iIndex = std::max_element(vGoodChannels.begin(), vGoodChannels.end())-vGoodChannels.begin();
	fMax =vGoodChannels.at(iIndex);
	iNumGood = 0;
	for (int i = 0; i < vGoodChannels.size(); i++){
		if (vGoodChannels.at(i) > .2*fMax){
			vGoodChannels.at(i) = 1;
			iNumGood++;
		}
		else
			vGoodChannels.at(i) = 0;
	}
	// debug output
	for (size_t i = 0; i < vGoodChannels.size(); i++){
		if (bMatlab_){
//			mexPrintf("vGoodChannels[%i]: %f\n", i,vGoodChannels.at(i));mexEvalString("drawnow;");
		}
		else{
			GADGET_DEBUG2("vGoodChannels[%i]: %f\n", i,vGoodChannels.at(i));
		}
	}

	// get relevant image - dRelevantImg = squeeze(dImg(:, :, dPos, lGoodChannels));
	vStart.clear(); vStart.push_back(0); vStart.push_back(0); vStart.push_back(iMaxChan); vStart.push_back(0);
	vSize.clear(); vSize.push_back(aImg.get_size(0)); vSize.push_back(aImg.get_size(1)); vSize.push_back(1); vSize.push_back(1);
	hoNDArray<std::complex<float>> acRelevantImg(aImg.get_size(0), aImg.get_size(1), iNumGood);
	o = 0; //helper - only non-zero entries of vGoodChannels are of interest
	hoNDArray<std::complex<float>> afTmp;	
	for (int c = 0; c < iNoChannels_; c++){
		if (vGoodChannels.at(c) == 1){
			vStart.at(3) = c;
			get_subarray(aImg, vStart, vSize, afTmp);
			// fill array
			size_t offset = acRelevantImg.get_size(0)*acRelevantImg.get_size(1)*o;
			for (long i = 0; i < afTmp.get_number_of_elements(); i++){
				acRelevantImg.at(i + offset) = afTmp.at(i);
			}
			o++;
		}
	}
	
	//-------------------------------------------------------------------------
	if (bMatlab_){
//		mexPrintf("get SOSImg\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("get SOSImg\n");
	}
	hoNDArray<float> aRelevantImg(*acRelevantImg.get_dimensions());
	multiplyConj(acRelevantImg,acRelevantImg,acRelevantImg);
	// complex float to float datatype
	for (long i = 0; i < aRelevantImg.get_number_of_elements(); i++){
		aRelevantImg[i] = acRelevantImg[i].real();
	}
	
	// loop over channels
	size_t tOffset = aRelevantImg.get_size(0)*aRelevantImg.get_size(1);
	hoNDArray<float> hafMax(aRelevantImg.get_dimensions());
	float* fPtr = hafMax.get_data_ptr();
	for (int iI = 0; iI < aRelevantImg.get_size(2); iI++){
		hoNDArray<float> tmp(aRelevantImg.get_size(0), aRelevantImg.get_size(1), aRelevantImg.get_data_ptr()+tOffset*iI);
		iMaxIndex = amax(&tmp);
		fMax = tmp.at(iMaxIndex);
		for (int iL = 0; iL < tOffset; iL++){
			fPtr[iL+tOffset*iI] = fMax;
		}
	}
	divide(aRelevantImg, hafMax, aRelevantImg);
	
	hoNDArray<float> aSOSImg;
	sum_dim(aRelevantImg, 2, aSOSImg);
	sqrt(aSOSImg, aSOSImg);

	// convert array
	hoNDArray<std::complex<float>> cfaSOSImgTest(aSOSImg.get_dimensions()); cfaSOSImgTest.fill(std::complex<float>(0.0, 0.0));
	std::complex<float> *cfPointer = cfaSOSImgTest.get_data_ptr();
	for (long iI = 0; iI < cfaSOSImgTest.get_number_of_elements(); iI++){
		cfPointer[iI] = std::complex<float>(aSOSImg.at(iI), 0.0);
	}

	//-------------------------------------------------------------------------
	// get navigator
	int iDisplacementMax  = 80; // [mm] diaphragm displacement (max. +/- 80mm)
	int iDisplacement = (int)std::ceil((float)iDisplacementMax/((float)field_of_view_[0]/(float)aImg.get_size(0))-.5);

	// fill last line in ref image
	hoNDArray<float> aRefImg = aSOSImg; aRefImg.fill(0.0);
	for (int iI = 0; iI < aSOSImg.get_size(0); iI++)
		aRefImg.at(iI + aSOSImg.get_size(0)*(aSOSImg.get_size(1)-1)) = aSOSImg.at(iI + aSOSImg.get_size(0)*(aSOSImg.get_size(1)-1));
	
	// create index vector and navigator vector (filled with zeros)
	std::vector<int> vIdx;
	for (int i = 0; i < aRefImg.get_size(1); i++) vIdx.push_back(i);
	vNav_.clear();
	for (int i = 0; i < aRefImg.get_size(1); i++) vNav_.push_back(0);
		
	// performance measurement
	double fClockSumDim = 0.0, fClockCirc = 0.0, fClockMultiply = 0.0, fClockFirst = 0.0, fClockSec = 0.0, fClockMem = 0.0, fClockSumDim2 = 0.0, fClockTh = 0.0;
	clock_t b;
	hoNDArray<float> aRMSImg;
	for (int i = aRefImg.get_size(1)-2; i >1;i--){// -1 ; i--){
		aRMSImg.create(aRefImg.get_size(0), 2*iDisplacement+1); aRMSImg.fill(0.0);	

		//MATLAB: tmp = dRefImg(:,idx(i+1:end))
		hoNDArray<float> tmp;
		vStart.clear(); vStart.push_back(0); vStart.push_back(i+1); vSize.clear(); vSize.push_back(aRefImg.get_size(0)); vSize.push_back(aRefImg.get_size(1)-i-1);
		get_subarray(aRefImg, vStart, vSize, tmp);
		hoNDArray<float> repTmp(tmp.get_dimensions());

		//MATLAB: dSOSImg(:,idx(i)
		hoNDArray<float> aTmp;
		vStart.clear(); vStart.push_back(0); vStart.push_back(vIdx.at(i)); vSize.clear(); vSize.push_back(aSOSImg.get_size(0)); vSize.push_back(1);		
		get_subarray(aSOSImg, vStart, vSize, aTmp);

		hoNDArray<float> aTmp2 = aTmp;
		circshift(aTmp2, -iDisplacement-1, 0);		
		for (int l = -iDisplacement; l <= iDisplacement; l++){
			
			//MATLAB: circshift(dSOSImg(:,idx(i)), iD)				
			circshift(aTmp2, 1, 0);	

			//MATLAB: (tmp - repmat(circshift(dSOSImg(:,idx(i)), iD),[1 size(tmp,2)]))
			hoNDArray<float> tmp2(tmp.get_dimensions());tmp2.fill(0.0); // result of subtraction	
			int N = tmp.get_size(0), LE = tmp.get_size(1); float *pA = tmp.begin(), *pB = aTmp2.begin(), *pR = tmp2.begin();				
			#pragma omp parallel for default(none) schedule(static) shared(N, pA, pB, pR, LE)
			for (int iL = 0; iL < LE; iL++)
				for (int iE = 0; iE < N; iE++)
					pR[iE + N*iL] = pA[iE + N*iL] - pB[iE];	

			//MATLAB: (tmp - repmat(circshift(dSOSImg(:,idx(i)), iD),[1 size(tmp,2)])).^2			
			multiply(tmp2, tmp2, tmp2);

			//MATLAB: sum((tmp - repmat(circshift(dSOSImg(:,idx(i)), iD),[1 size(tmp,2)])).^2,2)
			hoNDArray<float> hTmp2;
			if (tmp2.get_number_of_dimensions() > 1){
				std::vector<size_t> vD = *tmp2.get_dimensions();
				vD.pop_back();
				hTmp2.create(&vD); hTmp2.fill(0.0);
				float *pNewArray	= hTmp2.get_data_ptr(), *pOldArray	= tmp2.get_data_ptr();	
				int N = tmp2.get_size(1), L = hTmp2.get_number_of_elements();

				#pragma omp parallel for default(none) schedule(static) shared(N, L, pNewArray, pOldArray)
				for (int sum_dim = 0; sum_dim < N; sum_dim++)
					for (int i = 0; i < L; i++)
						pNewArray[i] += pOldArray[i + L*sum_dim];
			}
			
			//MATLAB: dRMSImg(:,dDisplacement-iD+1) = sum((tmp - repmat(circshift(dSOSImg(:,idx(i)), iD),[1 size(tmp,2)])).^2,2)			
			int iOffset = (iDisplacement-l)*hTmp2.get_size(0); 			
			memcpy(aRMSImg.begin() + iOffset, hTmp2.begin(), sizeof(float)*hTmp2.get_size(0));
		}
	

		//MATLAB: sum(dRMSImg(dX-round(dDisplacement/2):dX+round(dDisplacement/2),:))
		vStart.clear(); vStart.push_back(dX-(int)(((float)iDisplacement/2)+.5)); vStart.push_back(0); vSize.clear(); vSize.push_back(iDisplacement+1); vSize.push_back(aRMSImg.get_size(1));	
		get_subarray(aRMSImg, vStart, vSize, aRMSImg);

		hoNDArray<std::complex<float>> cfaRMSImgTest(aRMSImg.get_dimensions()); cfaRMSImgTest.fill(std::complex<float>(0.0, 0.0));
		std::complex<float> *cfPointer = cfaRMSImgTest.get_data_ptr();
		for (long iI = 0; iI < cfaRMSImgTest.get_number_of_elements(); iI++){
			cfPointer[iI] = std::complex<float>(aRMSImg.at(iI), 0.0);
		}

		sum_dim(aRMSImg, 0, aRMSImg);

		cfaRMSImgTest.clear();
		cfaRMSImgTest.create(aRMSImg.get_dimensions()); cfaRMSImgTest.fill(std::complex<float>(0.0, 0.0));
		cfPointer = cfaRMSImgTest.get_data_ptr();
		for (long iI = 0; iI < cfaRMSImgTest.get_number_of_elements(); iI++){
			cfPointer[iI] = std::complex<float>(aRMSImg.at(iI), 0.0);
		}

		//MATLAB: min(sum(dRMSImg(dX-round(dDisplacement/2):dX+round(dDisplacement/2),:)))
		int iMinVal = amin(&aRMSImg);		

		//MATLAB: dDisplacement + 1 - dNav(i)
		vNav_.at(i) =iDisplacement - iMinVal;

		//MATLAB: circshift(dSOSImg(:,idx(i)), dNav(i))
		hoNDArray<float> aTmp3 = aTmp;
		circshift(aTmp3, vNav_.at(i), 0);
		if (i%20 == 0){
			if (bMatlab_){
//				mexPrintf("Getting Navigator - %.1f %%\n", (float)(aRefImg.get_size(1)-2-i)/(float)(aRefImg.get_size(1)-2)*100);mexEvalString("drawnow;");
			}
			else{
				GADGET_DEBUG2("Getting Navigator - %.1f %%\n", (float)(aRefImg.get_size(1)-2-i)/(float)(aRefImg.get_size(1)-2)*100);		
			}
		}	

		//MATLAB: dRefImg(:,idx(i)) = circshift(dSOSImg(:,idx(i)), dNav(i))
		memcpy(aRefImg.get_data_ptr()+vIdx.at(i)*aRefImg.get_size(0), aTmp3.get_data_ptr(), sizeof(float)*aTmp3.get_size(0));
	}

	for (long i = 0; i < vNav_.size(); i++) vNav_.at(i) *= -1;

	// get Gaussian filter kernel and calculate convolution with navigator data
	vGaussian.clear();
	filter1DGaussian(vGaussian, 5);
	vectorConv(vNav_, vGaussian, 0);
	
	//-------------------------------------------------------------------------
	// interpolate navigator data signal to TR intervals
	if (bMatlab_){
//		mexPrintf("interpolation of navigator data to TR intervals..\n");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("interpolation of navigator data to TR intervals..\n");
	}
	for (long i = 0; i < vNav_.size(); i++)
		vNav_.at(i) = -vNav_.at(i);
	
	int iMin = std::min_element(vNav_.begin(), vNav_.end())-vNav_.begin();
	float fMin = vNav_.at(iMin);
	for (long i = 0; i < vNav_.size(); i++)
		vNav_.at(i) -= fMin;

	//interpolation - vNavInt_
	std::vector<float> vNavIndNew;
	for (long i = 0; i < lNoScans_; i++) vNavIndNew.push_back(i);
		
	if (bMatlab_){
//		mexPrintf("vNavInd size: %i, vNav_ size: %i, vNavIndNew size: %i\n", GlobalVar::instance()->vNavInd_.size(), vNav_.size(), vNavIndNew.size());mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG2("vNavInd size: %i, vNav_ size: %i, vNavIndNew size: %i\n", GlobalVar::instance()->vNavInd_.size(), vNav_.size(), vNavIndNew.size());
	}
	std::vector<float> vNavInd = GlobalVar::instance()->vNavInd_;
	vNavInt_ = interp1(vNavInd, vNav_, vNavIndNew);

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_NavigatorGadget)
}
