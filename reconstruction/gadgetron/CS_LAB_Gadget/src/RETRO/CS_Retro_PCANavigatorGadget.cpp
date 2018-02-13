#include "CS_Retro_PCANavigatorGadget.h"

namespace Gadgetron{
	// class constructor
	CS_Retro_PCANavigatorGadget::CS_Retro_PCANavigatorGadget() {

	}

	// class destructor - delete temporal buffer/memory
	CS_Retro_PCANavigatorGadget::~CS_Retro_PCANavigatorGadget(){

	}

	// read flexible data header
	int CS_Retro_PCANavigatorGadget::process_config(ACE_Message_Block* mb) {
		return GADGET_OK;
	}

	int CS_Retro_PCANavigatorGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, GadgetContainerMessage< hoNDArray <std::complex<float> > >* m3) {
		iNoChannels_	= m1->getObjectPtr()->channels;
		iNoNav_			= m1->getObjectPtr()->user_int[5];
		lNoScans_		= m3->getObjectPtr()->get_size(1);

		if (getNav2DPCA(*m2->getObjectPtr())) {
			if (bMatlab_) {
				//mexPrintf("Error in getNav2D\n");mexEvalString("drawnow;");
			} else {
				GADGET_DEBUG1("Error in getNav2D\n");
			}
		}

		GadgetContainerMessage< hoNDArray<float> >* tmp_m2 = new GadgetContainerMessage< hoNDArray<float> >();

		// dumb: set vNavInt_ to GlobalVar::vNavInd_. Maybe this is the same?
		// TODO: Delete vNavInt_ on success!
		vNavInt_ = GlobalVar::instance()->vNavInd_;

		tmp_m2->getObjectPtr()->create(vNavInt_.size());

		// convert vector to array
		float* fPtr = tmp_m2->getObjectPtr()->get_data_ptr();
		for (long iI = 0; iI < vNavInt_.size(); iI++) {
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
	bool CS_Retro_PCANavigatorGadget::getNav2DPCA(hoNDArray<std::complex<float>> &aNav) {
		if (bMatlab_) {
			//mexPrintf("\n\n**************************************\n********** get navigator 2D **********\n**************************************\n\n");mexEvalString("drawnow;");
		} else {
			GADGET_DEBUG1("\n\n**************************************\n********** get navigator 2D **********\n**************************************\n\n");
		}

		// reconstruct the 1-D projections for all measurements and all channels
		if (bMatlab_) {
			//mexPrintf("domain transformation - k-space to image\n");mexEvalString("drawnow;");
		} else {
			GADGET_DEBUG1("domain transformation - k-space to image\n");
		}

		/*assumptions:
		* iMeasurementTime_ is the total scan time in seconds
		* lNoScans_ is the same as length(iLC(:,15)) in Matlab
		* */

		int iNSamples =			aNav.get_size(0);
		int iNMeasurment =		aNav.get_size(1);
		int iNavRes	 =			aNav.get_size(2);
		int iNChannels	 =		aNav.get_size(3);
		double pi = 3.14159265358979323846;

		/* MATLAB
		% Reconstruct the 1-D projections for all measurements and all channels
		dImg = fftshift(ifft(ifftshift(dKSpace)));
		*/
		hoNDArray<std::complex<float>> aImg = aNav;
		hoNDFFT_CS<float>::instance()->ifftshift3D(aImg);
		hoNDFFT_CS<float>::instance()->ifft1(aImg);
		hoNDFFT_CS<float>::instance()->fftshift3D(aImg);

		//Concatenate arrays
		std::vector<size_t> A_dims;
		A_dims.push_back(iNSamples*iNavRes*iNChannels);
		A_dims.push_back(iNMeasurment);

		hoNDArray< std::complex<float> > A;
		A.create(&A_dims);

		for (unsigned long int i = 0; i < (iNSamples*iNavRes*iNChannels*iNMeasurment); i++) {
			A.at(i) = aImg.at(i);
		}

		hoNDKLT< std::complex<float> > * VT = new hoNDKLT < std::complex<float> >;
		std::vector<size_t> coeff_dims;
		coeff_dims.push_back(iNMeasurment);
		coeff_dims.push_back(iNMeasurment);

		hoNDArray< std::complex<float> > coeff;
		coeff.create(&coeff_dims);

		//Compute PCA based on KLT principal components are saved in coeff in descending order
		VT->prepare(A, static_cast<size_t>(1), static_cast<size_t>(0), true);
		VT->eigen_vector(coeff);

		double Fs = static_cast<float>(GlobalVar::instance()->iMeasurementTime_)/(static_cast<float>(iNMeasurment)*1000.0); // Get the sampling frequency

		//fft(result, ...,1);
		//fft only 1 dimensional(first dimension)
		hoNDFFT<float>::instance()->fft1(coeff);

		//nfft2 = 2.^nextpow2(nfft);
		int nfft2 = std::pow(2, std::ceil(log(coeff.get_size(0))/log(2)));
		hoNDArray< std::complex<float> > absresult;
		absresult.create(&coeff_dims);

		//fy = fft(y, nfft2);
		//fy = abs(fy);
		for(int i =0; i < coeff.get_number_of_elements(); i++) {
			absresult.at(i) = abs(coeff.at(i));
		}

		int inspectednr = 15; // only search in the first 15 principal components
		std::vector<size_t> fy_dims;
		fy_dims.push_back(inspectednr);
		fy_dims.push_back(nfft2/2-1);
		hoNDArray< std::complex<float> > fy;
		fy.create(&fy_dims);

		//fy = fy(1:nfft2/2);
		for(int x = 0; x < fy.get_size(0); x++) {
			for(int i = 0; i < fy.get_size(1); i++) {
				fy.at((i)+(x*fy.get_size(1))) = absresult.at((i)+(x*fy.get_size(1)));
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
		int frequency = 0, searcharea = std::floor(Fu)-std::floor(Fl), colmnnr =0;
		for(int x = 0; x < inspectednr; x++) {
			for(int i = 0; i < searcharea; i++) {
				if(compare_complex_values(maxvalue, fy.at(((i+std::floor(Fl)))+(x*iNMeasurment))) < 0) {
					maxvalue = fy.at(((i+std::floor(Fl)))+(x*iNMeasurment));
					frequency = i+1;
					colmnnr = x;
				}
			}
		}

		//frequency =  ((Fl)+frequency-2)*Fs/nfft2;
		float realfrequency =  (Fl + frequency - 2) * Fs / nfft2;

		//dECG = real(coeff(:,coeffnumber)) - imag(coeff(:,coeffnumber));
		// Fs= 1/dTR  changing sampling rate because signal is going to be interpolated
		Fs = static_cast<float>(GlobalVar::instance()->iMeasurementTime_)/(static_cast<float>(lNoScans_)*1000.0);
		std::vector<size_t> dECG_dims;
		dECG_dims.push_back(iNMeasurment);
		dECG_dims.push_back(1);

		hoNDArray< std::complex<float> > dECGtemp;
		dECGtemp.create(&dECG_dims);

		for(int i = 0; i < iNMeasurment; i++) {
			dECGtemp.at(i) = coeff.at(i+(colmnnr * iNMeasurment));
		}

		//get the real and the imag part of the signal and subtract them.
		hoNDArray< std::complex<float> > dECGhoNDArray;
		dECGhoNDArray.create(&dECG_dims);
		hoNDArray< std::complex<float> > realpart;
		realpart.create(&dECG_dims);

		// extract real part, for newer compilers also consider:
		//realpart = real(&dECGtemp);
		for(int i = 0; i < iNMeasurment; i++) {
			realpart.at(i) = dECGtemp.at(i).real();
		}

		hoNDArray< std::complex<float> > imaginarypart;
		imaginarypart.create(&dECG_dims);

		// extract imaginary part, for newer compilers also consider:
		//imaginarypart = imag(&dECGtemp);
		for(int i = 0; i < iNMeasurment; i++) {
			imaginarypart.at(i) = dECGtemp.at(i).imag();
		}

		subtract(&realpart, &imaginarypart, &dECGhoNDArray);
		std::vector<std::complex<float> > dECG;
		for (long i = 0; i < iNMeasurment; i++) {
			dECG.push_back(0);
		}

		//type conversion from complex to float and to vector
		for(int i = 0; i<iNMeasurment; i++) {
			dECG.at(i) = real(dECGhoNDArray.at(i));
		}

		//factor = length(iLC)/size(coeff,1);
		//dECG = fScale(dECG , factor);
		std::vector<std::complex<float> > dECGIndtemp;
		for (long i = 1; i <= lNoScans_; i++) {
			dECGIndtemp.push_back(i);
		}

		std::vector<std::complex<float> > dECGInt;
		for (long i = 0; i < iNMeasurment; i++) {
			dECGInt.push_back(i*lNoScans_/iNMeasurment);
		}

		std::vector<std::complex<float> > dECGInd;
		dECGInd = interp1<std::complex<float> >(dECGInt, dECG, dECGIndtemp);

		// Filter the Signal with a first order butterworth filter

		//===============================================================
		//calculate the numerator and denominator of a first order butterworth filter. (End of calulation is indicated by =======)
		//===============================================================

		//ul = 4*tan(pi*fl/2);
		//uh = 4*tan(pi*fh/2);
		//den = [1 0 0];
		float ul = 4*tan(pi*(realfrequency-0.1)/Fs/2/2);
		float uh = 4*tan(pi*(realfrequency+0.1)/Fs/2/2);
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
		hoNDArray< std::complex<float> > t1;
		t1.create(&t1_dims);
		hoNDArray< std::complex<float> > t2;
		t2.create(&t1_dims);
		hoNDArray< std::complex<float> > ad;
		ad.create(&t1_dims);

		t1.at(0) = 1+(Wn*(-Bw/Wn)/4);
		t1.at(1) = Wn/4;
		t1.at(2) = -Wn/4;
		t1.at(3) = 1;

		// also transpose
		t2.at(0)= 1-(Wn*(-Bw/Wn)/4);
		t2.at(2)= -Wn/4;//(switched indices because of transpose)
		t2.at(1)= Wn/4;
		t2.at(3)= 1;

		multiply(t2, t1, ad);

		//%den = poly(ad);
		//e = eig(ad);
		std::vector<size_t> e_dims;
		e_dims.push_back(2);
		e_dims.push_back(1);
		hoNDArray< std::complex<float> > e;
		e.create(&e_dims);

		std::vector<size_t> kern_dims;
		kern_dims.push_back(3);
		kern_dims.push_back(1);
		hoNDArray< std::complex<float> > kern;
		kern.create(&kern_dims);
		std::vector<std::complex<float> > num;
		num.push_back(0);
		num.push_back(0);
		num.push_back(0);

		//eig
		hoNDKLT< std::complex<float> > * Vtemp = new hoNDKLT < std::complex<float> >;

		//get the eigen value
		Vtemp->prepare(ad, static_cast<size_t>(1), static_cast<size_t>(0), false);
		Vtemp->eigen_value(e);

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
		for(int k = 0; k<3; k++) {
			std::complex<float> i = std::complex<float>(0.0, 1.0);		// define imaginary constant: make code compliant to older compiler versions
			kern.at(k) = std::exp(-i*static_cast<std::complex<float> >(Wn)*static_cast<std::complex<float> >(k));
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
		std::vector<std::complex<float> > Y;
		for (long i = 0; i < dECG.size(); i++) {
			Y.push_back(0);
		}

		//for m = 1:length(Y)
		//  Y(m) = num(1) * X(m) + z(1);
		//   for i = 2:n
		//      z(i - 1) = num(i) * X(m) + z(i) - den(i) * Y(m);
		//   end
		//end
		for(int m = 0; m < Y.size(); m++) {
			Y.at(m) = num.at(0) * dECG.at(m) + z.at(0);

			for(int i = 1; i < den.size(); i++) {
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
		Y.clear();
		for (long i = 0; i < dECG.size(); i++) {
			Y.push_back(0);
		}

		// second round filtering (backward)
		//for m = 1:length(Y)
		//   Y(m) = b(1) * X(m) + z(1);
		//   for i = 2:n
		//      z(i - 1) = b(i) * X(m) + z(i) - a(i) * Y(m);
		//   end
		//end
		for(int m = 0; m < Y.size(); m++) {
			Y.at(m) = num.at(0) * dECG.at(m) + z.at(0);

			for(int i = 1; i < den.size(); i++) {
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
		for(int i = 0; i<dECG.size()-1; i++) {
			dECGdiff.push_back(dECG.at(i+1)-dECG.at(i));
		}

		return GADGET_OK;
	}

	GADGET_FACTORY_DECLARE(CS_Retro_PCANavigatorGadget)
} // close namespace Gadgetron
