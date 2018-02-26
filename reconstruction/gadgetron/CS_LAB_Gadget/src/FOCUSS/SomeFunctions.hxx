/*
file name	: 	SomeFunctions.hxx

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	implementation of the prototypes in SomeFunctions.h

reference	:	hoNDArray_math_util.h from the Gadgetron implementation
*/

#include "SomeFunctions.h"

namespace Gadgetron
{
	// r = a-b*c - all arrays same size
	template <typename T>
	bool fAminusBmultC(hoNDArray<T> &a, hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r)
	{
		// check if number of elements are equal
		GADGET_DEBUG_CHECK_RETURN_FALSE(a.get_number_of_elements()==b.get_number_of_elements());
		GADGET_DEBUG_CHECK_RETURN_FALSE(b.get_number_of_elements()==c.get_number_of_elements());

		const T *pA = a.begin();
		const T *pB = b.begin();
		const T *pC = c.begin();
		T *pR = r.begin();

		try {
			if (r.get_number_of_elements()!=a.get_number_of_elements()) {
				r = a;
			}

			long long N = (long long)a.get_number_of_elements();
			long long n;

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pR)
			for (n=0; n < (long long)N; n++) {
				pR[n] = pA[n];
			}

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pB, pC, pR)
			for (n=0; n < (long long)N; n++) {
				pR[n] -= pB[n]*pC[n];
			}
		} catch (...) {
			GERROR("Error occurred in AminusBmultC(const hoNDArray<T> &a, const hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r) ... ");
			return false;
		}

		return true;
	}

	// r = a*b-c
	template <typename T>
	bool fAmultBminusC(const hoNDArray<T> &a, const hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r)
	{
		// check if number of elements are equal
		GADGET_DEBUG_CHECK_RETURN_FALSE(a.get_number_of_elements()==b.get_number_of_elements());
		GADGET_DEBUG_CHECK_RETURN_FALSE(b.get_number_of_elements()==c.get_number_of_elements());

		const T*pA = a.begin();
		const T*pB = b.begin();
		const T*pC = c.begin();
		T *pR = r.begin();

		long long N = a.get_number_of_elements();

		try {
			if (r.get_number_of_elements()!=a.get_number_of_elements()) {
				r.create(a.get_dimensions());
			}

			long long n;

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pB, pR)
			for (n=0; n < (long long)N; n++) {
				pR[n] = pA[n]*pB[n];
			}

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pC, pR)
			for (n=0; n < (long long)N; n++) {
				pR[n] -= pC[n];
			}

		} catch (...) {
			GERROR("Error occurred in AmultBminusC(const hoNDArray<T> &a, const hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r) ... ");
			return false;
		}

		return true;
	}

	// r = a+b*c
	template <typename T>
	bool fAplusBmultC(const hoNDArray<T> &a, const hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r)
	{
		// check if number of elements are equal
		GADGET_DEBUG_CHECK_RETURN_FALSE(a.get_number_of_elements()==b.get_number_of_elements());
		GADGET_DEBUG_CHECK_RETURN_FALSE(b.get_number_of_elements()==c.get_number_of_elements());

		const T*pA = a.begin();
		const T*pB = b.begin();
		const T*pC = c.begin();
		T *pR = r.begin();

		long long N = a.get_number_of_elements();

		try {
			if (r.get_number_of_elements() != a.get_number_of_elements()) {
				r.create(a.get_dimensions());
			}

			long long n;

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pB, pC, pR)
			for (n=0; n < (long long)N; n++) {
				pR[n] = pA[n];
			}

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pB, pC, pR)
			for (n=0; n < (long long)N; n++) {
				pR[n] += pB[n]*pC[n];
			}
		} catch (...) {
			GERROR("Error occurred in AplusBmultC(const hoNDArray<T> &a, const hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r) ... ");
			return false;
		}

		return true;
	}

	// r = a + b; a: constant value, b: array
	template <typename T>
	bool fAplusB(T a, const hoNDArray<T> &b, hoNDArray<T> &r)
	{
		const T*pA = a.begin();
		const T*pB = b.begin();
		T *pR = r.begin();

		long long N = b.get_number_of_elements();

		try {
			if (r.get_number_of_elements()!=b.get_number_of_elements()) {
				r.create(b.get_dimensions());
			}

			long long n;

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, a, pB, pR)
			for (n=0; n < (long long)N; n++) {
				pR[n] = a+pB[n];
			}
		} catch (...) {
			GERROR("Error occurred in AplusB(T a, const hoNDArray<T> &b, hoNDArray<T> &r) ... ");
			return false;
		}

		return true;
	}

	// calculate abs and pow to b
	template <typename T>
	bool fAbsPow(hoNDArray<T> &a, float b)
	{
		T*pA = a.begin();

		long long N = a.get_number_of_elements();

		try {
			long long n;

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, b)
			for (n=0; n < (long long)N; n++) {
				pA[n] = std::pow(std::abs(pA[n]), b);
			}
		} catch (...) {
			GERROR("Error occurred in pow_abs(hoNDArray<T> &a, const T b) ... ");
			return false;
		}

		return true;
	}

	// calculate abs, pow to b and divide by signal energy
	template <typename T>
	bool fAbsPowDivide(hoNDArray<T> &a, float b, const hoNDArray<T> &c)
	{
		T*pA = a.begin();
		const T* pC = c.begin();

		long long N = a.get_number_of_elements();

		try {
			long long n;

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, b)
			for (n=0; n < (long long)N; n++) {
				pA[n] = std::pow(std::abs(pA[n]), b);
			}

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pC)
			for (n=0; n < (long long)N; n++) {
				pA[n] /= pC[n];
			}
		} catch (...) {
			GERROR("Error occurred in pow_abs(hoNDArray<T> &a, const T b) ... ");
			return false;
		}

		return true;
	}

	// multiply: a = a.*b
	template <typename T>
	bool fMultiply(hoNDArray<T> &a, const hoNDArray<T> &b)
	{
		T*pA = a.begin();
		const T* pB = b.begin();

		long long N = a.get_number_of_elements();

		try {
			long long n;

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pB)
			for (n=0; n < (long long)N; n++) {
				pA[n] *= pB[n];
			}
		} catch (...) {
			GERROR("Error occurred in pow_abs(hoNDArray<T> &a, const T b) ... ");
			return false;
		}

		return true;
	}

	// G: r = -conj(a)*b+c*d+e*f
	template <typename T>
	bool fCalcGradient(const hoNDArray<T> &a, const hoNDArray<T> &b, T c, const hoNDArray<T> &d, T e, const hoNDArray<T> &f, hoNDArray<T> &r)
	{
		// check if number of elements are equal
		GADGET_DEBUG_CHECK_RETURN_FALSE(a.get_number_of_elements()==b.get_number_of_elements());
		GADGET_DEBUG_CHECK_RETURN_FALSE(b.get_number_of_elements()==d.get_number_of_elements());
		GADGET_DEBUG_CHECK_RETURN_FALSE(d.get_number_of_elements()==f.get_number_of_elements());

		const T *pA = a.begin();
		const T *pD = d.begin();
		const T *pF = f.begin();
		T *pR = r.begin();

		long long N = b.get_number_of_elements();

		try {
			if (r.get_number_of_elements()!=a.get_number_of_elements()) {
				r.create(a.get_dimensions());
			}

			long long n;
			r=b;

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pR)
			for (n=0; n < (long long)N; n++) {
				pR[n] *= -std::conj(pA[n]);
			}

			if (c != std::complex<float>(0.0)) {
				//GDEBUG("calc gradient- lambda: %f + i %f..\n", c.real(), c.imag());
				#pragma omp parallel for default(none) schedule(static) private(n) shared(N, c, pD, pR)
				for (n=0; n < (long long)N; n++) {
					pR[n] += c*pD[n];
				}
			}

			if (e != std::complex<float>(0.0)) {
				#pragma omp parallel for default(none) schedule(static) private(n) shared(N, e, pF, pR)
				for (n=0; n < (long long)N; n++) {
					pR[n] += e*pF[n];
				}
			}
		} catch (...) {
			GERROR("Error occurred in calc_gradient(...) .. ");
			return false;
		}

		return true;
	}

	//ESPReSSo: 0.5 * Lambda * (conj(W) .* q .* (1+phase) + conj(W).*W.*conj(q).*phase - fTrafkSpaceCombi.*(1+phase));
	template <typename T>
	bool fESPReSSoOut(const hoNDArray<T> &W, const hoNDArray<T> &q, const hoNDArray<T> &phase, const hoNDArray<T> &ifft_kSpaceCombi, hoNDArray<T> &out)
	{
		// check if number of elements are equal
		GADGET_DEBUG_CHECK_RETURN_FALSE(W.get_number_of_elements()==q.get_number_of_elements());
		GADGET_DEBUG_CHECK_RETURN_FALSE(q.get_number_of_elements()==phase.get_number_of_elements());
		GADGET_DEBUG_CHECK_RETURN_FALSE(phase.get_number_of_elements()==ifft_kSpaceCombi.get_number_of_elements());

		const T* pW = W.begin();
		const T* pQ = q.begin();
		const T* pPh = phase.begin();
		const T* pkS = ifft_kSpaceCombi.begin();
		T* pO = out.begin();

		std::complex<float> cVal(1.0);
		std::complex<float> cVal2(0.5);

		long long N = W.get_number_of_elements();

		try {
			if (out.get_number_of_elements() != W.get_number_of_elements()) {
				out.create(W.get_dimensions());
			}

			long long n;

			#pragma omp parallel for default(none) private(n) shared(N, pW, pQ, pO)
			for (n = 0; n < N; n++) {
				pO[n] = std::conj(pW[n])*pQ[n];
			}

			#pragma omp parallel for default(none) private(n) shared(N, pPh, pO, cVal)
			for (n = 0; n < N; n++) {
				pO[n] *= (cVal + pPh[n]);
			}

			#pragma omp parallel for default(none) private(n) shared(N, pW, pQ, pPh, pO)
			for (n = 0; n < N; n++) {
				pO[n] += std::conj(pW[n])*pW[n]*std::conj(pQ[n])*pPh[n];
			}

			#pragma omp parallel for default(none) private(n) shared(N, pPh, pkS, pO, cVal)
			for (n = 0; n < N; n++) {
				pO[n] -= pkS[n]*(cVal + pPh[n]);
			}

			#pragma omp parallel for default(none) private(n) shared(N, pO, cVal2)
			for (n = 0; n < N; n++) {
				pO[n] *= cVal2;
			}
		} catch (...) {
			GERROR("Error occurred in ESPReSSoOut(const hoNDArray<T> &W, const hoNDArray<T> &q, const hoNDArray<T> &phase, const hoNDArray<T> &ifft_kSpaceCombi, hoNDArray<T> &out) ... ");
			return false;
		}

		return true;
	}

	// crop center part with knowledge that x isn't PF - dim: y-z-x-c
	template <typename T>
	bool fCropArrYZ(const hoNDArray<T> &Array, int a, int b, hoNDArray<T> &result)
	{
		int y = a;
		int z = b;
		// get dims of incoming array
		std::vector<size_t> dims = *Array.get_dimensions();

		// get center indices
		size_t iCenterX = std::floor((float)(dims[2]-1)/2);
		size_t iCenterY = std::floor((float)(dims[0]-1)/2);
		size_t iCenterZ = std::floor((float)(dims[1]-1)/2);

		// get border indices
		//GDEBUG("y: %i, z: %i\n", y,z);
		size_t iBorderY_Upper, iBorderY_Lower, iBorderZ_Upper, iBorderZ_Lower;
		if (y%2 == 0) {
			iBorderY_Upper = iCenterY+std::ceil((float)y/2);
			iBorderY_Lower = iCenterY-std::floor((float)y/2)+1;
		} else {
			iBorderY_Upper = iCenterY+std::ceil((float)y/2)-1;
			iBorderY_Lower = iCenterY-std::floor((float)y/2);
		}

		if (z%2 == 0) {
			iBorderZ_Upper = iCenterZ+std::ceil((float)z/2);
			iBorderZ_Lower = iCenterZ-std::floor((float)z/2)+1;
		} else{
			iBorderZ_Upper = iCenterZ+std::ceil((float)z/2)-1;
			iBorderZ_Lower = iCenterZ-std::floor((float)z/2);
		}

		//GDEBUG("Border: y: %i,%i z: %i,%i\n", iBorderY_Upper, iBorderY_Lower, iBorderZ_Upper, iBorderZ_Lower);
		// indices out of bounds?
		if (iBorderY_Upper > dims[0]-1) {
			iBorderY_Upper = dims[0]-1;
			y = dims[0];
		}

		if (iBorderZ_Upper > dims[1]-1) {
			iBorderZ_Upper = dims[1]-1;
			z = dims[1];
		}

		if (iBorderY_Lower < 0) {
			iBorderY_Lower = 0;
		}

		if (iBorderZ_Lower < 0) {
			iBorderZ_Lower = 0;
		}

		// create output
		result.create(y,z);

		// fill output
		T* PtrIn = Array.get_data_ptr();
		T* PtrOut= result.get_data_ptr();
		for (size_t iy = iBorderY_Lower; iy <= iBorderY_Upper; iy++) {
			for(size_t iz = iBorderZ_Lower; iz <= iBorderZ_Upper; iz++) {
				PtrOut[(iy-iBorderY_Lower)+(iz-iBorderZ_Lower)*y] = PtrIn[iy + iz*dims[0] + iCenterX*dims[0]*dims[1]];
			}
		}

		return true;
	}

	// calculate signal energy of array
	template <typename T>
	float fCalcEnergy(hoNDArray<T> a)
	{
		const T *pA = a.begin();
		hoNDArray<float> r(a.get_dimensions());

		float * pR = r.begin();
		float fResult = 0.0;

		try {
			long long N = (long long)a.get_number_of_elements();
			long long n;

			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pR)
			for (n=0; n < (long long)N; n++) {
				pR[n] = std::pow(abs(pA[n]),2);
			}

			for (n=0; n < (long long)N; n++) {
				fResult += pR[n];
			}
		} catch (...) {
			GERROR("Error occurred in AminusBmultC(const hoNDArray<T> &a, const hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r) ... ");
			return false;
		}

		return fResult;
	}

	inline bool fAllZero(const hoNDArray<bool> &Array)
	{
		bool isAllZero = true;
		bool *PtrIn = Array.get_data_ptr();

		for (size_t i = 0; i < Array.get_number_of_elements(); i++) {
			if (PtrIn[i] != false) {
				isAllZero = false;
				break;
			}
		}

		return isAllZero;
	}

	inline bool fAllZero(const hoNDArray<std::complex<float> > &Array)
	{
		bool isAllZero = true;
		std::complex<float> *PtrIn = Array.get_data_ptr();

		for (size_t i = 0; i < Array.get_number_of_elements(); i++) {
			if (PtrIn[i] != std::complex<float>(0, 0)) {
				isAllZero = false;
				break;
			}
		}

		return isAllZero;
	}

	// function, which checks if all elements are one (bool)
	inline bool fAllOne(const hoNDArray<bool> &Array)
	{
		bool *PtrIn = Array.get_data_ptr();

		for (size_t i = 0; i < Array.get_number_of_elements(); i++) {
			if (PtrIn[i] == false) {
				return false;
			}
		}

		return true;
	}

	// calculate L-point hanning window with L = elements-1
	/*

	implementation reference:


	%von-Hann (Hanning) window
	%
	% w = hann(L) returns an L-point symmetric Hann window in the column vector w.
	% L must be a positive integer.
	%
	% The coefficients of a Hann window are computed from the following equation:
	%
	%      w(n) = 0.5 * (1 + cos(2*pi*n/N)),   0 <= n <= N
	%
	% The window length is L = N+1.
	%
	% w = hann(L,'sflag') returns an L-point Hann window using the window sampling specified by 'sflag',
	% which can be either 'periodic' or 'symmetric' (the default). The 'periodic' flag is useful for DFT/FFT purposes,
	% such as in spectral analysis.
	% The DFT/FFT contains an implicit periodic extension and the periodic flag enables a signal windowed
	% with a periodic window to have perfect periodic extension.
	% When 'periodic' is specified, hann computes a length L+1 window and returns the first L points.
	% When using windows for filter design, the 'symmetric' flag should be used.
	%
	% --> http://www.mathworks.de/de/help/signal/ref/hann.html
	% --> https://de.wikipedia.org/wiki/Hann-Fenster

	% implemented by Michael.Voelker@mr-bavaria.de, 2012
	*/
	inline std::vector<float>* fGetHanningWindow(int iElements)
	{
		const float pi = 3.141592653589793238463;
		std::vector<float>* vfHanning = new std::vector<float>(iElements);

		// L = N - 1
		float fL = float(iElements-1);

		// calc n
		for (std::vector<int>::size_type i = 0; i != vfHanning->size(); i++) {
			vfHanning->at(i) = float(i)-fL/2;
		}

		// calc w(n)
		for (std::vector<int>::size_type i = 0; i != vfHanning->size(); i++) {
			vfHanning->at(i) = .5*(1+std::cos(2*pi*vfHanning->at(i)/fL));
		}

		return vfHanning;
	}

	// calculate L-point hamming window with L = N (periodic case)
	/*

	implementation reference:

	%Hamming window
	%
	% w = hamming(L) returns an L-point symmetric Hamming window in the column vector w.
	% L should be a positive integer.
	%
	%  The coefficients of a Hamming window are computed from the following equation:
	%
	%       w(n) = 0.54  +  0.46 * cos(2*pi*n/N),   0 <= n <= N
	%
	%
	% w = hamming(L, 'symFlag') returns an L-point Hamming window using the window sampling
	% specified by 'symFlag', which can be either 'periodic'  or 'symmetric' (the default).
	% The 'periodic' flag is useful for DFT/FFT purposes, such as in spectral analysis.
	% The DFT/FFT contains an implicit periodic extension and the periodic flag enables a signal
	% windowed with a periodic window to have perfect periodic extension.
	% When 'periodic' is specified, hamming computes a length L+1 window and returns the first L points.
	% When using windows for filter design, the 'symmetric' flag should be used.
	%
	% --> http://www.mathworks.de/de/help/signal/ref/hamming.html
	% --> https://de.wikipedia.org/wiki/Hamming-Fenster

	% implemented by Michael.Voelker@mr-bavaria.de, 2012*/
	inline std::vector<float>* fGetHammingWindow(int iElements)
	{
		const float pi = 3.141592653589793238463;
		std::vector<float> *vfHamming = new std::vector<float>(iElements);

		// L = N - periodic case
		float fL = float(iElements);

		// calc n
		for (std::vector<int>::size_type i = 0; i != vfHamming->size(); i++) {
			vfHamming->at(i) = float(i)-fL/2;
		}

		// calc w(n)
		for (std::vector<int>::size_type i = 0; i != vfHamming->size(); i++) {
			vfHamming->at(i) = 0.54 + 0.46*std::cos(2*pi*vfHamming->at(i)/fL);
		}

		return vfHamming;
	}

	// sum array in specified dimension and squeeze the result - compare to MATLAB sum(array, dim)
	template <typename T>
	bool sum_dim(hoNDArray<T> &Array, unsigned int dimension, hoNDArray<T> &result)
	{
		hoNDArray<T> NewTmpArray = Array;

		// get number of dimensions
		size_t iNoDim = NewTmpArray.get_number_of_dimensions();

		// get dimensions
		std::vector<size_t> vDims = *NewTmpArray.get_dimensions();

		// check dimensions
		if (Array.get_number_of_dimensions() <= dimension) {
			BOOST_THROW_EXCEPTION(runtime_error("sum_dim: Error occurred - array is smaller than expected by the user!\n"));
		} else {
			if (iNoDim > 1) {
				// permute array - summation dimension will be last dimension
				std::vector<size_t> vPermute;
				for (size_t i = 0; i < iNoDim; i++) {
					if (i != dimension) {
						vPermute.push_back(i);
					}
				}

				vPermute.push_back(dimension);
				NewTmpArray = *permute(&NewTmpArray, &vPermute,false);

				// create output array
				std::vector<size_t> vNewDims = *NewTmpArray.get_dimensions();
				vNewDims.pop_back();
				hoNDArray<T> tmp(&vNewDims); tmp.fill(T(0.0));

				// summation in specified dimension
				T* pNewArray = tmp.get_data_ptr();
				T* pOldArray = NewTmpArray.get_data_ptr();
				for (size_t sum_dim = 0; sum_dim < vDims.at(dimension); sum_dim++) {
					for (size_t i = 0; i < tmp.get_number_of_elements(); i++) {
						pNewArray[i] += pOldArray[i + tmp.get_number_of_elements()*sum_dim];
					}
				}

				// output
				result = tmp;
			} else {
				// create output
				hoNDArray<T> tmp(1);

				// summation
				T* pOldArray = NewTmpArray.get_data_ptr();
				T* pNewArray = tmp.get_data_ptr();
				pNewArray[0] = 0;
				for (size_t sum_dim = 0; sum_dim < vDims.at(0); sum_dim++) {
					pNewArray[0] += pOldArray[sum_dim];
				}

				// output
				result = tmp;
			}
		}

		return GADGET_OK;
	}

	// copy all AcquisitionHeader values
	//inline int fCopyHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1, GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1_new) {
	//	GC_acq_m1_new->getObjectPtr()->acquisition_time_stamp		= GC_acq_m1->getObjectPtr()->acquisition_time_stamp;
	//	GC_acq_m1_new->getObjectPtr()->active_channels				= GC_acq_m1->getObjectPtr()->active_channels;
	//	GC_acq_m1_new->getObjectPtr()->available_channels			= GC_acq_m1->getObjectPtr()->available_channels;
	//	GC_acq_m1_new->getObjectPtr()->center_sample				= GC_acq_m1->getObjectPtr()->center_sample;
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[0]				= GC_acq_m1->getObjectPtr()->channel_mask[0];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[1]				= GC_acq_m1->getObjectPtr()->channel_mask[1];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[2]				= GC_acq_m1->getObjectPtr()->channel_mask[2];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[3]				= GC_acq_m1->getObjectPtr()->channel_mask[3];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[4]				= GC_acq_m1->getObjectPtr()->channel_mask[4];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[5]				= GC_acq_m1->getObjectPtr()->channel_mask[5];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[6]				= GC_acq_m1->getObjectPtr()->channel_mask[6];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[7]				= GC_acq_m1->getObjectPtr()->channel_mask[7];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[8]				= GC_acq_m1->getObjectPtr()->channel_mask[8];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[9]				= GC_acq_m1->getObjectPtr()->channel_mask[9];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[10]				= GC_acq_m1->getObjectPtr()->channel_mask[10];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[11]				= GC_acq_m1->getObjectPtr()->channel_mask[11];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[12]				= GC_acq_m1->getObjectPtr()->channel_mask[12];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[13]				= GC_acq_m1->getObjectPtr()->channel_mask[13];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[14]				= GC_acq_m1->getObjectPtr()->channel_mask[14];
	//	GC_acq_m1_new->getObjectPtr()->channel_mask[15]				= GC_acq_m1->getObjectPtr()->channel_mask[15];
	//	GC_acq_m1_new->getObjectPtr()->discard_post					= GC_acq_m1->getObjectPtr()->discard_post;
	//	GC_acq_m1_new->getObjectPtr()->discard_pre					= GC_acq_m1->getObjectPtr()->discard_pre;
	//	GC_acq_m1_new->getObjectPtr()->encoding_space_ref			= GC_acq_m1->getObjectPtr()->encoding_space_ref;
	//	GC_acq_m1_new->getObjectPtr()->flags						= GC_acq_m1->getObjectPtr()->flags;
	//	GC_acq_m1_new->getObjectPtr()->idx.average					= GC_acq_m1->getObjectPtr()->idx.average;
	//	GC_acq_m1_new->getObjectPtr()->idx.contrast					= GC_acq_m1->getObjectPtr()->idx.contrast;
	//	GC_acq_m1_new->getObjectPtr()->idx.kspace_encode_step_1		= GC_acq_m1->getObjectPtr()->idx.kspace_encode_step_1;
	//	GC_acq_m1_new->getObjectPtr()->idx.kspace_encode_step_2		= GC_acq_m1->getObjectPtr()->idx.kspace_encode_step_2;
	//	GC_acq_m1_new->getObjectPtr()->idx.phase					= GC_acq_m1->getObjectPtr()->idx.phase;
	//	GC_acq_m1_new->getObjectPtr()->idx.repetition				= GC_acq_m1->getObjectPtr()->idx.repetition;
	//	GC_acq_m1_new->getObjectPtr()->idx.segment					= GC_acq_m1->getObjectPtr()->idx.segment;
	//	GC_acq_m1_new->getObjectPtr()->idx.set						= GC_acq_m1->getObjectPtr()->idx.set;
	//	GC_acq_m1_new->getObjectPtr()->idx.slice					= GC_acq_m1->getObjectPtr()->idx.slice;
	//	GC_acq_m1_new->getObjectPtr()->idx.user[0]					= GC_acq_m1->getObjectPtr()->idx.user[0];
	//	GC_acq_m1_new->getObjectPtr()->idx.user[1]					= GC_acq_m1->getObjectPtr()->idx.user[1];
	//	GC_acq_m1_new->getObjectPtr()->idx.user[2]					= GC_acq_m1->getObjectPtr()->idx.user[2];
	//	GC_acq_m1_new->getObjectPtr()->idx.user[3]					= GC_acq_m1->getObjectPtr()->idx.user[3];
	//	GC_acq_m1_new->getObjectPtr()->idx.user[4]					= GC_acq_m1->getObjectPtr()->idx.user[4];
	//	GC_acq_m1_new->getObjectPtr()->idx.user[5]					= GC_acq_m1->getObjectPtr()->idx.user[5];
	//	GC_acq_m1_new->getObjectPtr()->idx.user[6]					= GC_acq_m1->getObjectPtr()->idx.user[6];
	//	GC_acq_m1_new->getObjectPtr()->idx.user[7]					= GC_acq_m1->getObjectPtr()->idx.user[7];
	//	GC_acq_m1_new->getObjectPtr()->measurement_uid				= GC_acq_m1->getObjectPtr()->measurement_uid;
	//	GC_acq_m1_new->getObjectPtr()->number_of_samples			= GC_acq_m1->getObjectPtr()->number_of_samples;
	//	GC_acq_m1_new->getObjectPtr()->patient_table_position[0]	= GC_acq_m1->getObjectPtr()->patient_table_position[0];
	//	GC_acq_m1_new->getObjectPtr()->patient_table_position[1]	= GC_acq_m1->getObjectPtr()->patient_table_position[1];
	//	GC_acq_m1_new->getObjectPtr()->patient_table_position[2]	= GC_acq_m1->getObjectPtr()->patient_table_position[2];
	//	GC_acq_m1_new->getObjectPtr()->phase_dir[0]					= GC_acq_m1->getObjectPtr()->phase_dir[0];
	//	GC_acq_m1_new->getObjectPtr()->phase_dir[1]					= GC_acq_m1->getObjectPtr()->phase_dir[1];
	//	GC_acq_m1_new->getObjectPtr()->phase_dir[2]					= GC_acq_m1->getObjectPtr()->phase_dir[2];
	//	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[0]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[0];
	//	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[1]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[1];
	//	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[2]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[2];
	//	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[3]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[3];
	//	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[4]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[4];
	//	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[5]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[5];
	//	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[6]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[6];
	//	GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[7]		= GC_acq_m1->getObjectPtr()->physiology_time_stamp[7];
	//	GC_acq_m1_new->getObjectPtr()->position[0]					= GC_acq_m1->getObjectPtr()->position[0];
	//	GC_acq_m1_new->getObjectPtr()->position[1]					= GC_acq_m1->getObjectPtr()->position[1];
	//	GC_acq_m1_new->getObjectPtr()->position[2]					= GC_acq_m1->getObjectPtr()->position[2];
	//	GC_acq_m1_new->getObjectPtr()->read_dir[0]					= GC_acq_m1->getObjectPtr()->read_dir[0];
	//	GC_acq_m1_new->getObjectPtr()->read_dir[1]					= GC_acq_m1->getObjectPtr()->read_dir[1];
	//	GC_acq_m1_new->getObjectPtr()->read_dir[2]					= GC_acq_m1->getObjectPtr()->read_dir[2];
	//	GC_acq_m1_new->getObjectPtr()->sample_time_us				= GC_acq_m1->getObjectPtr()->sample_time_us;
	//	GC_acq_m1_new->getObjectPtr()->scan_counter					= GC_acq_m1->getObjectPtr()->scan_counter;
	//	GC_acq_m1_new->getObjectPtr()->slice_dir[0]					= GC_acq_m1->getObjectPtr()->slice_dir[0];
	//	GC_acq_m1_new->getObjectPtr()->slice_dir[1]					= GC_acq_m1->getObjectPtr()->slice_dir[1];
	//	GC_acq_m1_new->getObjectPtr()->slice_dir[2]					= GC_acq_m1->getObjectPtr()->slice_dir[2];
	//	GC_acq_m1_new->getObjectPtr()->trajectory_dimensions		= GC_acq_m1->getObjectPtr()->trajectory_dimensions;
	//	GC_acq_m1_new->getObjectPtr()->user_float[0]				= GC_acq_m1->getObjectPtr()->user_float[0];
	//	GC_acq_m1_new->getObjectPtr()->user_float[1]				= GC_acq_m1->getObjectPtr()->user_float[1];
	//	GC_acq_m1_new->getObjectPtr()->user_float[2]				= GC_acq_m1->getObjectPtr()->user_float[2];
	//	GC_acq_m1_new->getObjectPtr()->user_float[3]				= GC_acq_m1->getObjectPtr()->user_float[3];
	//	GC_acq_m1_new->getObjectPtr()->user_float[4]				= GC_acq_m1->getObjectPtr()->user_float[4];
	//	GC_acq_m1_new->getObjectPtr()->user_float[5]				= GC_acq_m1->getObjectPtr()->user_float[5];
	//	GC_acq_m1_new->getObjectPtr()->user_float[6]				= GC_acq_m1->getObjectPtr()->user_float[6];
	//	GC_acq_m1_new->getObjectPtr()->user_float[7]				= GC_acq_m1->getObjectPtr()->user_float[7];
	//	GC_acq_m1_new->getObjectPtr()->user_int[0]					= GC_acq_m1->getObjectPtr()->user_int[0];
	//	GC_acq_m1_new->getObjectPtr()->user_int[1]					= GC_acq_m1->getObjectPtr()->user_int[1];
	//	GC_acq_m1_new->getObjectPtr()->user_int[2]					= GC_acq_m1->getObjectPtr()->user_int[2];
	//	GC_acq_m1_new->getObjectPtr()->user_int[3]					= GC_acq_m1->getObjectPtr()->user_int[3];
	//	GC_acq_m1_new->getObjectPtr()->user_int[4]					= GC_acq_m1->getObjectPtr()->user_int[4];
	//	GC_acq_m1_new->getObjectPtr()->user_int[5]					= GC_acq_m1->getObjectPtr()->user_int[5];
	//	GC_acq_m1_new->getObjectPtr()->user_int[6]					= GC_acq_m1->getObjectPtr()->user_int[6];
	//	GC_acq_m1_new->getObjectPtr()->user_int[7]					= GC_acq_m1->getObjectPtr()->user_int[7];
	//	GC_acq_m1_new->getObjectPtr()->version						= GC_acq_m1->getObjectPtr()->version;
	//
	//	return GADGET_OK;
	//}


#ifndef __GADGETRON_VERSION_HIGHER_3_6__
	// save array
	inline bool save_array(hoNDArray< std::complex<float> > &Array, std::string file_prefix)
	{
		std::string file_path_;
		std::string ismrmrd_file_name_;
		boost::shared_ptr<ISMRMRD::IsmrmrdDataset> ismrmrd_dataset_;

		//Generate filename
		ismrmrd_file_name_ = file_prefix + std::string(".h5");
		ISMRMRD::HDF5Exclusive lock; //This will ensure threadsafe access to HDF5
		ismrmrd_dataset_ = boost::shared_ptr<ISMRMRD::IsmrmrdDataset>(new ISMRMRD::IsmrmrdDataset(ismrmrd_file_name_.c_str(), "dataset"));

		std::vector<size_t> dim = *Array.get_dimensions();

		GINFO("Try to write array on disk..\n");

		try {
			//ismrmrd_dataset_->appendImageHeader(*m1->getObjectPtr(),"image.head");
			hoNDArray< std::complex<float> >* buffer_ = new hoNDArray< std::complex<float> >(dim, Array.get_data_ptr(),false);
			std::vector<unsigned int> dims(dim.size());

			size_t i;
			for (i = 0; i< dim.size(); i++) {
				dims[i] = dim[i];
			}

			if (ismrmrd_dataset_->appendArray(dims, Array.get_data_ptr(), "image_0.img") < 0) {
				GERROR("Failed to write image data\n");
				return GADGET_FAIL;
			}
		} catch (...) {
			GERROR("Error attempting to append images to HDF5 file\n");
			return GADGET_FAIL;
		}

		return GADGET_OK;
	}
#endif

	// flip array in specified dimension - reference to: hoNDFFT.cpp (cpufft - original Gadgetron)
	template <typename T>
	bool flip_array(hoNDArray<T> &Array, unsigned int dimension)
	{
		// check dimensions
		if (Array.get_number_of_dimensions() <= dimension) {
			BOOST_THROW_EXCEPTION(runtime_error("flip_array: Error occurred - array is smaller than expected by the user!\n"));
			return false;
		}

		int stride     = 1;           //Distance between points in transform
		int dist       = 1;           //Distance between vectors
		int trafos     = 1;           //Transformations per chunk
		int chunks     = 1;           //Number of chunks
		int chunk_size = 1;           //Points per chunk
		int length     = 1;           //Length of each transform
		int total_dist = 1;

		hoNDArray<T> buffer;
		const T* data_ptr = Array.begin();
		buffer.create(Array.get_dimensions());
		T* buffer_ptr = buffer.begin();

		//Set sizes
		length = (int)Array.get_size(dimension);

		if (dimension != 0) {
			for (size_t i = 0; i < dimension; i++) {
				chunk_size *= (int)Array.get_size(i);
			}

			stride = chunk_size;
			trafos = chunk_size;
			chunk_size *= length;

			for (size_t i = dimension+1; i < Array.get_number_of_dimensions(); i++) {
				chunks *= (int)Array.get_size(i);
			}
		} else {
			for (size_t i = 1; i < Array.get_number_of_dimensions(); i++) {
				trafos *= (int)Array.get_size(i);
			}

			chunk_size = trafos*length;

			dist = length;
		}

		total_dist = trafos*dist;

		int idx1_max = chunks*chunk_size;
		int idx1, idx2;       //Index variables
		int idx2_limit;
		int length2 = length;
		int stride2 = stride;

		//Loop over all chunks
		for (idx1 = 0; idx1 < idx1_max; idx1+=chunk_size) {
			idx2_limit = idx1+total_dist;

			//Loop over all transformations
			for (idx2 = idx1; idx2 < idx2_limit; idx2+=dist) {
				///Copy data to buffer.
				int idx3 = idx2;
				for (int j = 0; j < length2; idx3+=stride2) {
					buffer_ptr[idx3] = data_ptr[idx2 + stride2*(length2 - j - 1)];
					j++;
				}
			} //Loop over transformations
		} //Loop over chunks

		// output
		Array = buffer;

		return true;
	}

	// flip line in array by specified offset
	template <typename T>
	bool flip_line(hoNDArray<T> &Array, size_t offset)
	{
		//temporal array
		hoNDArray<T> aTmp(Array.get_size(0));

		// fill temp array
		memcpy(aTmp.begin(), Array.begin() + offset, sizeof(T)*Array.get_size(0));

		// flip array
		flip_array(aTmp, 0);

		// copy back
		memcpy(Array.begin()+offset, aTmp.begin(), sizeof(T)*Array.get_size(0));

		return GADGET_OK;
	}

	// sum array in specified dimension and squeeze the result - compare to MATLAB sum(array, dim)
	/*template <typename T>
	bool sum_dim(hoNDArray<T> &Array, int dimension, hoNDArray<T> &result) {

		hoNDArray<T> NewTmpArray = Array;

		// get number of dimensions
		int iNoDim = NewTmpArray.get_number_of_dimensions();

		// get dimensions
		std::vector<size_t> vDims = *NewTmpArray.get_dimensions();

		// check dimensions
		if (Array.get_number_of_dimensions() <= dimension) {
			BOOST_THROW_EXCEPTION(runtime_error("sum_dim: Error occurred - array is smaller than expected by the user!\n"));
		}
		else{
			if (iNoDim > 1) {
				// permute array - summation dimension will be last dimension
				std::vector<size_t> vPermute;
				for (int i = 0; i < iNoDim; i++)
					if (i != dimension)
						vPermute.push_back(i);

				vPermute.push_back(dimension);
				NewTmpArray = *permute(&NewTmpArray, &vPermute,false);

				// number of elements of all dimension except the summation dimension
				long lMax = NewTmpArray.get_number_of_elements()/vDims.at(dimension);
				//GDEBUG("lMax: %i\n", lMax);
				// create output array
				std::vector<size_t> vNewDims = *NewTmpArray.get_dimensions();
				vNewDims.pop_back();
				hoNDArray<T> tmp(&vNewDims); tmp.fill(T(0.0));

				// summation in specified dimension
				T* pNewArray = tmp.get_data_ptr();
				T* pOldArray = NewTmpArray.get_data_ptr();
				for (int sum_dim = 0; sum_dim < vDims.at(dimension); sum_dim++)
					for (int i = 0; i < tmp.get_number_of_elements(); i++)
						pNewArray[i] += pOldArray[i + tmp.get_number_of_elements()*sum_dim];

				// output
				result = tmp;
			}
			else{
				// create output
				hoNDArray<T> tmp(1);
				
				// summation
				T* pOldArray = NewTmpArray.get_data_ptr();
				T* pNewArray = tmp.get_data_ptr(); pNewArray[0] = 0;
				for (long sum_dim = 0; sum_dim < vDims.at(0); sum_dim++)
					pNewArray[0] += pOldArray[sum_dim];

				// output
				result = tmp;
			}
		}

		return GADGET_OK;
	}*/

	// create symmetric Gaussian 1D kernel
	template <typename T>
	void filter1DGaussian(std::vector<T> &result, int length)
	{
		// get corner points
		int iStart	= -std::floor((float)length/2);
		int iEnd	= +std::floor((float)length/2);

		//GDEBUG("iStart = %i, iEnd = %i\n", iStart, iEnd);

		// get Gaussian kernel values and push on result vector
		for (int i = iStart; i <= iEnd; i++) {
			result.push_back(exp(-.5*std::pow((float)i*6/length,2)));
		}

		// get sum
		T sum = 0.0;
		for (size_t i = 0; i < result.size(); i++) {
			sum += result.at(i);
		}

		// norm output
		for (size_t i = 0; i < result.size(); i++) {
			result.at(i) /= sum;
		}
	}

	// vector convolution - output same size (option - 0: same size, 1: full convolution)
	template <typename T>
	void vectorConv(std::vector<T> &vector, std::vector<T> &kernelVec, int option)
	{
		// get vector size
		int	N = vector.size();
		int K = kernelVec.size();

		// output size - full convolution
		int L = K + N - 1;

		// temporal vector
		std::vector<T> vTmp(L);

		// full convolution
		for (int l = 0; l < L; l++) {
			int kMin, kMax, k;

			// initialize output element
			vTmp[l] = 0;

			// determine min and max indices
			kMin = (l >= K - 1) ? l - (K - 1) : 0;
			kMax = (l < N - 1) ? l : N - 1;

			// loop over kernel
			for (k = kMin; k <= kMax; k++) {
				vTmp[l] += vector[k] * kernelVec[l-k];
			}
		}

		vector = vTmp;

		// if option 'same size' - cut out from center
		if (option == 0) {
			std::vector<T> vTmp2(N);

			// offset
			int iOffset = std::floor(((float)L/2 - (float)N/2));

			memcpy(&vTmp2[0], &vTmp[0]+iOffset, sizeof(T)*N);

			// fill output
			vector = vTmp2;
		}
	}

	// 1D array convolution - first dimension - same size
	template <typename T>
	void arrayConv(hoNDArray<T> &Array, std::vector<T> &kernelVec)
	{
		// temporal Array
		hoNDArray<T> aOutput(*Array.get_dimensions());

		// temporal vector
		std::vector<T> vTmp(aOutput.get_size(0));

		// loop over lines
		for (size_t i = 0; i < aOutput.get_size(1); i++) {
			// offset of current line
			int offset = i*aOutput.get_size(0);

			// copy data into temp vec
			memcpy(&vTmp[0], Array.begin()+offset, sizeof(T)*(aOutput.get_size(0)));

			// do convolution
			vectorConv(vTmp, kernelVec, 0);

			// copy result into output array
			memcpy(aOutput.begin()+offset, &vTmp[0], sizeof(T)*(aOutput.get_size(0)));
		}

		// fill output array
		Array = aOutput;

		//aOutput.clear();
	}

	// 1D array convolution - same size - arbitrary dimension - ref to: hoNDFFT.cpp from original Gadgetron implementation
	template <typename T>
	void arrayConv(hoNDArray<T> &Array, std::vector<T> &kernelVec, unsigned int dimension)
	{
		// check dimensions
		if (Array.get_number_of_dimensions() <= dimension) {
			BOOST_THROW_EXCEPTION(runtime_error("arrayConv(...): error occurred - array is smaller than expected by the user!\n"));
		}

		int stride     = 1;           //Distance between points in transform
		int dist       = 1;           //Distance between vectors
		int trafos     = 1;           //Transformations per chunk
		int chunks     = 1;           //Number of chunks
		int chunk_size = 1;           //Points per chunk
		int length     = 1;           //Length of each transform
		int total_dist = 1;

		hoNDArray<T> buffer;
		const T* data_ptr = Array.begin();
		buffer.create(Array.get_dimensions());
		T* buffer_ptr = buffer.begin();

		//Set sizes
		length = (int)Array.get_size(dimension);

		if (dimension != 0) {
			for (size_t i = 0; i < dimension; i++) {
				chunk_size *= (int)Array.get_size(i);
			}

			stride = chunk_size;
			trafos = chunk_size;
			chunk_size *= length;

			for (size_t i = dimension+1; i < Array.get_number_of_dimensions(); i++) {
				chunks *= (int)Array.get_size(i);
			}
		} else {
			for (size_t i = 1; i < Array.get_number_of_dimensions(); i++) {
				trafos *= (int)Array.get_size(i);
			}

			chunk_size = trafos*length;

			dist = length;
		}

		total_dist = trafos*dist;

		int idx1_max = chunks*chunk_size;
		int idx1, idx2;       //Index variables
		int idx2_limit;
		int length2 = length;
		int stride2 = stride;

		std::vector<T> vTmp;

		//Loop over all chunks
		for (idx1 = 0; idx1 < idx1_max; idx1+=chunk_size) {
			idx2_limit = idx1+total_dist;

			//Loop over all transformations
			for (idx2 = idx1; idx2 < idx2_limit; idx2+=dist) {
				///Copy data to buffer.
				int idx3 = idx2;
				for (int j = 0; j < length2; idx3+=stride2) {
					vTmp.push_back(data_ptr[idx3]);
					j++;
				}

				// do convolution
				vectorConv(vTmp, kernelVec, 0);

				// copy data to temporal output
				idx3 = idx2;
				for (int j = 0; j < length2; idx3+=stride2) {
					buffer_ptr[idx3] = vTmp.at(j);
					j++;
				}

				// clear temporal vector
				vTmp.clear();

			} //Loop over transformations
		} //Loop over chunks

		// output
		Array = buffer;
	}

	// circluar 1D shift - only positive shifts
	template <typename T>
	void circshift(std::vector<T> &Array, int shift)
	{
		// out of boundary?
		if (static_cast<unsigned int>(std::abs(shift)) > Array.size()) {
			if (shift >= 0) {
				shift = shift%Array.size();
			} else {
				shift = -std::abs(shift)%Array.size();
			}
		}

		// temporal vector
		std::vector<T> vOutput(Array.size());

		// copy data into new vector
		// positive shift
		if (shift >= 0) {
			memcpy(&vOutput[0]+shift, &Array[0], sizeof(T)*(Array.size()-shift));
			memcpy(&vOutput[0], &Array[0]+(Array.size()-shift), sizeof(T)*shift);
		} else {
			// negative shift
			memcpy(&vOutput[0], &Array[0]-shift, sizeof(T)*(Array.size()+shift));
			memcpy(&vOutput[0]+(Array.size()+shift), &Array[0], sizeof(T)*(-shift));
		}

		// copy to input vector
		Array = vOutput;

		vOutput.clear();
	}

	// circular 1D array shift (arbitrary dimension) - ref to: hoNDFFT.cpp from original Gadgetron implementation
	template <typename T>
	void circshift(hoNDArray<T> &Array, int shift, unsigned int dimension)
	{
		// check dimensions
		if (Array.get_number_of_dimensions() <= dimension) {
			BOOST_THROW_EXCEPTION(runtime_error("circshift(...): error occurred - array is smaller than expected by the user!\n"));
		}

		// out of boundary?
		if (static_cast<unsigned int>(std::abs(shift)) > Array.get_size(dimension)) {
			if (shift >= 0) {
				shift = shift%Array.get_size(dimension);
			} else {
				shift = -std::abs(shift)%Array.get_size(dimension);
			}
		}

		int stride     = 1;           //Distance between points in transform
		int dist       = 1;           //Distance between vectors
		int trafos     = 1;           //Transformations per chunk
		int chunks     = 1;           //Number of chunks
		int chunk_size = 1;           //Points per chunk
		int length     = 1;           //Length of each transform
		int total_dist = 1;

		hoNDArray<T> buffer;
		const T* data_ptr = Array.begin();
		buffer.create(Array.get_dimensions());
		T* buffer_ptr = buffer.begin();

		//Set sizes
		length = (int)Array.get_size(dimension);

		if (dimension != 0) {
			for (size_t i = 0; i < dimension; i++) {
				chunk_size *= (int)Array.get_size(i);
			}

			stride = chunk_size;
			trafos = chunk_size;
			chunk_size *= length;

			for (size_t i = dimension+1; i < Array.get_number_of_dimensions(); i++) {
				chunks *= (int)Array.get_size(i);
			}
		} else {
			for (size_t i = 1; i < Array.get_number_of_dimensions(); i++) {
				trafos *= (int)Array.get_size(i);
			}

			chunk_size = trafos*length;
			dist = length;
		}

		total_dist = trafos*dist;

		int idx1_max = chunks*chunk_size;
		int idx1, idx2;       //Index variables
		int idx2_limit;
		int length2 = length;
		int stride2 = stride;

		std::vector<T> vTmp;

		//Loop over all chunks
		for (idx1 = 0; idx1 < idx1_max; idx1+=chunk_size) {
			idx2_limit = idx1+total_dist;

			//Loop over all transformations
			for (idx2 = idx1; idx2 < idx2_limit; idx2+=dist) {
				///Copy data to buffer.
				int idx3 = idx2;
				for (int j = 0; j < length2; idx3+=stride2) {
					vTmp.push_back(data_ptr[idx3]);
					j++;
				}

				// circular shift
				circshift(vTmp, shift);

				// copy data to temporal output
				idx3 = idx2;
				for (int j = 0; j < length2; idx3+=stride2) {
					buffer_ptr[idx3] = vTmp.at(j);
					j++;
				}

				vTmp.clear();

			} //Loop over transformations
		} //Loop over chunks

		// output
		Array = buffer;
	}

	// output a linear equally spaced vector
	template <typename T>
	std::vector<T> &linspace(T fStart, T fEnd, int iElements)
	{
		// empty output vector
		std::vector<T> result;

		// distance between adjacent elements
		T fDistance = (fEnd-fStart)/(iElements-1);

		// fill output vector
		for (long i = 0; i < iElements; i++) {
			result.push_back(fStart + i*fDistance);
		}

		// check last element
		if (result.back() != fEnd) {
			BOOST_THROW_EXCEPTION(runtime_error("linspace: error occurred - last element not as expected!\n"));
		}

		return result;
	}

	// interp1 from http://stackoverflow.com/questions/9394867/c-implementation-of-matlab-interp1-function-linear-interpolation
	template <typename T>
	inline std::vector<T> interp1(std::vector<T> &x, std::vector<T> &y, std::vector<T> &x_new)
	{
		std::vector<T> y_new;
		y_new.reserve(x_new.size());

		std::vector<T> dx, dy, slope, intercept;
		dx.reserve(x.size());
		dy.reserve(x.size());
		slope.reserve(x.size());
		intercept.reserve(x.size());

		for (size_t i = 0; i < x.size(); ++i) {
			if (i < x.size()-1) {
				dx.push_back(x.at(i+1) - x.at(i));
				dy.push_back(y.at(i+1) - y.at(i));
				slope.push_back(dy.at(i) / dx.at(i));
				intercept.push_back(y.at(i) - x.at(i) * slope.at(i));
			} else {
				dx.push_back(dx.at(i-1));
				dy.push_back(dy.at(i-1));
				slope.push_back(slope.at(i-1));
				intercept.push_back(intercept.at(i-1));
			}
		}

		for (size_t i = 0; i < x_new.size(); ++i) {
			int idx = findNearestNeighbourIndex<T>(x_new.at(i), x);

			if (idx < 0) {
				GWARN("idx<0: Something went wrong!\n");
				continue;
			}

			y_new.push_back(slope.at(idx) * x_new.at(i) + intercept.at(idx));
		}

		return y_new;
	}

	template <typename T>
	inline int findNearestNeighbourIndex(T value, std::vector< T > &x)
	{
		float dist = FLT_MAX;
		int idx = -1;

		for (size_t i = 0; i < x.size(); ++i) {
			float newDist = std::abs(value - x[i]);
			if (newDist > 0 && newDist < dist) {
				dist = newDist;
				idx = static_cast<int>(i);
			}
		}

		return idx;
	}

	// interpolation - options: 0: 'linear' or 1: 'pchip'
	//template <typename T>
	//std::vector<T>& interp1(std::vector<T> &vInX, std::vector<T> &vInY, std::vector<T> vOutX, int option) {
	//
	//	// same size of input vectors
	//	if (vInX.size() != vInY.size()) {
	//		BOOST_THROW_EXCEPTION(runtime_error("interp1: error occurred - vector with different size!\n"));
	//	}
	//
	//	// create temporal output vector
	//	std::vector<T> vTmpOut(vOutX.size());
	//
	//	GINFO("inside interpolation 1 - no. elements: %i, %i, %i - option: %i\n", vInX.size(), vInY.size(), vOutX.size(), option);
	//	// linear interpolation
	//	if (option == 0) {
	//		// helper variables
	//		float fLower, fHigher;
	//	
	//		// loop over new coord vector
	//		for (long i = 0; i < vOutX.size(); i++) {
	//
	//			float fCoord = vOutX.at(i);
	//
	//			// find interval in old coord system
	//			std::vector<T> vTmp = vInX;
	//
	//			for (long n = 0; n < vTmp.size(); n++) {
	//				vTmp.at(n) = std::abs(vTmp.at(n)-fCoord);
	//			}
	//			int iIndexMin = std::min_element(vTmp.begin(), vTmp.end())-vTmp.begin();
	//			float fMin    = vInX.at(iIndexMin);
	//
	//			// descent
	//			float m;
	//
	//			// lower, same or higher than the given element
	//			// same
	//			if (fCoord == vInX.at(iIndexMin)) {
	//				vTmpOut.at(i) = vInY.at(iIndexMin);
	//			}
	//
	//			// lower
	//			else if (fCoord < vInX.at(iIndexMin)) {
	//				// out of boundaries?
	//				if (iIndexMin != 0) m = (vInY.at(iIndexMin)-vInY.at(iIndexMin-1))/(vInX.at(iIndexMin)-vInX.at(iIndexMin-1));
	//				else				m = (vInY.at(iIndexMin+1)-vInY.at(iIndexMin))/(vInX.at(iIndexMin+1)-vInX.at(iIndexMin));
	//
	//				// new point
	//				vTmpOut.at(i) = m*(fCoord - vInX.at(iIndexMin))+vInY.at(iIndexMin);
	//
	//			}
	//
	//			// higher
	//			else if (fCoord > vInX.at(iIndexMin)) {
	//				// out of boundaries
	//				if (iIndexMin != (vInX.size()-1)) m = (vInY.at(iIndexMin+1)-vInY.at(iIndexMin))/(vInX.at(iIndexMin+1)-vInX.at(iIndexMin));
	//				else 					      m = (vInY.at(iIndexMin)-vInY.at(iIndexMin-1))/(vInX.at(iIndexMin)-vInX.at(iIndexMin-1));
	//
	//				// new point
	//				vTmpOut.at(i) = m*(fCoord - vInX.at(iIndexMin))+vInY.at(iIndexMin);
	//			}
	//			else{
	//				GDEBUG("fCoord: %f, vInxAt: %f, iIndexMin: %i\n", fCoord, vInX.at(iIndexMin), iIndexMin);
	//				BOOST_THROW_EXCEPTION(runtime_error("interp1: Error occurred - element cannot be found\n"));
	//			}
	//		}
	//
	//	}
	//	// piecewise constant interpolation
	//	else if(option == 1) {
	//		GERROR("interp1: Error occurred - option is not implemented in current version!\n");
	//	}
	//	else{
	//		BOOST_THROW_EXCEPTION(runtime_error("interp1: Error occurred - option not specified or unknown!\n"));
	//	}
	//
	//	return vTmpOut;
	//}



	// all elements true
	inline bool allTrue(hoNDArray<bool> &Array)
	{
		bool bFound = true;

		// get data pointer
		bool *pData = Array.get_data_ptr();

		// loop over array
		for (unsigned int i = 0; i < Array.get_number_of_elements(); i++) {
			if (pData[i] != true) {
				bFound = false;
				return bFound;
			}
		}

		return bFound;
	}

	inline bool allTrue(std::vector<bool> &Vector)
	{
		bool bFound = true;

		// loop over array
		for (unsigned int i = 0; i < Vector.size(); i++) {
			if (Vector[i] != true) {
				bFound = false;
				return bFound;
			}
		}

		return bFound;
	}

	// get sub array of 2D, 3D, 4D array
	template <typename T>
	void get_subarray(hoNDArray<T> &input, std::vector<size_t> vStart, std::vector<size_t> vSize, hoNDArray<T> &out)
	{
		// check vector size
		if (vStart.size() != vSize.size()) {
			BOOST_THROW_EXCEPTION(runtime_error("SomeFunctions::get_sub_array failed - number of dimensions unequal!\n"));
		}

		// out of bounds?
		for (size_t i = 0; i < vStart.size(); i++) {
			if ((vStart.at(i)+vSize.at(i)) > input.get_size(i)+1) {
				vStart.at(i) = input.get_size(i)-vSize.at(i);
				if (vStart.at(i) < 0) {
					vStart.at(i) = 0;
				}
	//			mexPrintf("SomeFunctions::get_sub_array failed - out of bounds - sub-array dimension size adjusted\n");mexEvalString("drawnow;");
			}
		}

		// create new array
		hoNDArray<T> aOut(&vSize);
		//input.print(std::cout);
		//aOut.print(std::cout);
		//for(size_t i = 0; i < aOut.get_number_of_dimensions(); i++) GDEBUG("dim %i: %i\n", i, aOut.get_size(i));
		try {
			// get sub array
			int iDims = vSize.size();
			switch (iDims) {
			case 2:
				for (size_t y = vStart.at(1), i = 0; y < vStart.at(1)+vSize.at(1); y++, i++) {
					size_t offset_new = i*vSize.at(0);
					size_t offset_old = y*input.get_size(0);
					memcpy(aOut.begin()+offset_new, input.begin()+offset_old+vStart.at(0), sizeof(T)*vSize.at(0));
				}
				break;
			case 3:
				for (size_t z = vStart.at(2), j = 0; z < vStart.at(2)+vSize.at(2); z++, j++) {
					for (size_t y = vStart.at(1), i = 0; y < vStart.at(1)+vSize.at(1); y++, i++) {
						size_t offset_new = i*vSize.at(0)+j*vSize.at(0)*vSize.at(1);
						size_t offset_old = y*input.get_size(0)+z*input.get_size(0)*input.get_size(1);
						memcpy(aOut.begin()+offset_new, input.begin()+offset_old+vStart.at(0), sizeof(T)*vSize.at(0));
					}
				}
				break;
			case 4:
				//GDEBUG("inside case 4\n");
				for (size_t c = vStart.at(3), k = 0; c < vStart.at(3) + vSize.at(3); c++, k++) {
					for (size_t z = vStart.at(2), j = 0; z < vStart.at(2)+vSize.at(2); z++, j++) {
						for (size_t y = vStart.at(1), i = 0; y < vStart.at(1)+vSize.at(1); y++, i++) {
							size_t offset_new = i*vSize.at(0)+j*vSize.at(0)*vSize.at(1)+k*vSize.at(0)*vSize.at(1)*vSize.at(2);
							size_t offset_old = y*input.get_size(0)+z*input.get_size(0)*input.get_size(1)+c*input.get_size(0)*input.get_size(1)*input.get_size(2);
							memcpy(aOut.begin()+offset_new, input.begin()+offset_old+vStart.at(0), sizeof(T)*vSize.at(0));
						}
					}
				}
				break;
			default:
				BOOST_THROW_EXCEPTION(runtime_error("SomeFunctions::get_sub_array failed - specified dimension not supported!\n"));
				break;
			}

			// return output
			out = aOut;
		} catch(...) {
			BOOST_THROW_EXCEPTION(runtime_error("SomeFunctions::get_subarray failed - catched error!\n"));
		}
	}

	template <typename T>
	bool sum_dim_g(hoNDArray<T> &Array, int dimension)
	{
		/*
		int stride     = 1;           //Distance between points in transform
		int dist       = 1;           //Distance between vectors
		int trafos     = 1;           //Transformations per chunk
		int chunks     = 1;           //Number of chunks
		int chunk_size = 1;           //Points per chunk
		int length     = 1;           //Length of each transform
		int total_dist = 1;

		std::vector<size_t> vDimsNew = *Array.get_dimensions();
		vDimsNew.at(dimension)=1;

		hoNDArray<T> buffer;
		const T* data_ptr = Array.begin();

		buffer.create(vDimsNew);
		T* buffer_ptr = buffer.begin();

		//Set sizes
		length = (int)Array.get_size(dimension);

		if (dimension != 0)
		{
			for (size_t i = 0; i < dimension; i++)
			{
				chunk_size *= (int)Array.get_size(i);
			}
			stride = chunk_size;
			trafos = chunk_size;
			chunk_size *= length;

			for (size_t i = dimension+1; i < Array.get_number_of_dimensions(); i++)
			{
				chunks *= (int)Array.get_size(i);
			}
		}
		else
		{
			for (size_t i = 1; i < Array.get_number_of_dimensions(); i++)
			{
				trafos *= (int)Array.get_size(i);
			}
			chunk_size = trafos*length;

			dist = length;
		}
		total_dist = trafos*dist;

		int idx1_max = chunks*chunk_size;
		int idx1, idx2;       //Index variables
		int idx2_limit;
		int length2 = length;
		int stride2 = stride;

		std::vector<T> vTmp;

		for (idx1 = 0; idx1 < idx1_max; idx1+=chunk_size) //Loop over all chunks
		{
			idx2_limit = idx1+total_dist;
			for (idx2 = idx1; idx2 < idx2_limit; idx2+=dist) //Loop over all transformations
			{
				///Copy data to buffer.
				{
					int j, idx3 = idx2;
					for (j = 0; j < length2; idx3+=stride2)
					{
						vTmp.push_back(data_ptr[idx3]);
						j++;
					}
				}

				T sum_elem = T(0.0);
				for (int iN = 0; iN < vTmp.size(); iN++)
					sum_elem += vTmp[iN];

				buffer_ptr[idx2] = sum_elem;

				vTmp.clear();

			} //Loop over transformations
		} //Loop over chunks

		// output
		Array = buffer;*/

		// TODO: Implement something
		return false;
	}

	inline int fCopyAcqHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1, GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1_new)
	{
		GC_acq_m1_new->getObjectPtr()->acquisition_time_stamp		= GC_acq_m1->getObjectPtr()->acquisition_time_stamp;
		GC_acq_m1_new->getObjectPtr()->active_channels				= GC_acq_m1->getObjectPtr()->active_channels;
		GC_acq_m1_new->getObjectPtr()->available_channels			= GC_acq_m1->getObjectPtr()->available_channels;
		GC_acq_m1_new->getObjectPtr()->center_sample				= GC_acq_m1->getObjectPtr()->center_sample;

		for (unsigned int i = 0; i < ISMRMRD::ISMRMRD_CHANNEL_MASKS; i++) {
			GC_acq_m1_new->getObjectPtr()->channel_mask[i]			= GC_acq_m1->getObjectPtr()->channel_mask[i];
		}

		GC_acq_m1_new->getObjectPtr()->discard_post					= GC_acq_m1->getObjectPtr()->discard_post;
		GC_acq_m1_new->getObjectPtr()->discard_pre					= GC_acq_m1->getObjectPtr()->discard_pre;
		GC_acq_m1_new->getObjectPtr()->encoding_space_ref			= GC_acq_m1->getObjectPtr()->encoding_space_ref;
		GC_acq_m1_new->getObjectPtr()->flags						= GC_acq_m1->getObjectPtr()->flags;
		GC_acq_m1_new->getObjectPtr()->idx.average					= GC_acq_m1->getObjectPtr()->idx.average;
		GC_acq_m1_new->getObjectPtr()->idx.contrast					= GC_acq_m1->getObjectPtr()->idx.contrast;
		GC_acq_m1_new->getObjectPtr()->idx.kspace_encode_step_1		= GC_acq_m1->getObjectPtr()->idx.kspace_encode_step_1;
		GC_acq_m1_new->getObjectPtr()->idx.kspace_encode_step_2		= GC_acq_m1->getObjectPtr()->idx.kspace_encode_step_2;
		GC_acq_m1_new->getObjectPtr()->idx.phase					= GC_acq_m1->getObjectPtr()->idx.phase;
		GC_acq_m1_new->getObjectPtr()->idx.repetition				= GC_acq_m1->getObjectPtr()->idx.repetition;
		GC_acq_m1_new->getObjectPtr()->idx.segment					= GC_acq_m1->getObjectPtr()->idx.segment;
		GC_acq_m1_new->getObjectPtr()->idx.set						= GC_acq_m1->getObjectPtr()->idx.set;
		GC_acq_m1_new->getObjectPtr()->idx.slice					= GC_acq_m1->getObjectPtr()->idx.slice;

		for (unsigned int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++) {
			GC_acq_m1_new->getObjectPtr()->idx.user[i]				= GC_acq_m1->getObjectPtr()->idx.user[i];
		}

		GC_acq_m1_new->getObjectPtr()->measurement_uid				= GC_acq_m1->getObjectPtr()->measurement_uid;
		GC_acq_m1_new->getObjectPtr()->number_of_samples			= GC_acq_m1->getObjectPtr()->number_of_samples;

		for (unsigned int i = 0; i < 3; i++) {
			GC_acq_m1_new->getObjectPtr()->patient_table_position[i]= GC_acq_m1->getObjectPtr()->patient_table_position[i];
		}

		for (unsigned int i = 0; i < 3; i++) {
			GC_acq_m1_new->getObjectPtr()->phase_dir[i]				= GC_acq_m1->getObjectPtr()->phase_dir[i];
		}

		for (unsigned int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++) {
			GC_acq_m1_new->getObjectPtr()->physiology_time_stamp[i]	= GC_acq_m1->getObjectPtr()->physiology_time_stamp[i];
		}

		for (unsigned int i = 0; i < 3; i++) {
			GC_acq_m1_new->getObjectPtr()->position[i]				= GC_acq_m1->getObjectPtr()->position[i];
		}

		for (unsigned int i = 0; i < 3; i++) {
			GC_acq_m1_new->getObjectPtr()->read_dir[i]				= GC_acq_m1->getObjectPtr()->read_dir[i];
		}

		GC_acq_m1_new->getObjectPtr()->sample_time_us				= GC_acq_m1->getObjectPtr()->sample_time_us;
		GC_acq_m1_new->getObjectPtr()->scan_counter					= GC_acq_m1->getObjectPtr()->scan_counter;

		for (unsigned int i = 0; i < 3; i++) {
			GC_acq_m1_new->getObjectPtr()->slice_dir[i]				= GC_acq_m1->getObjectPtr()->slice_dir[i];
		}

		GC_acq_m1_new->getObjectPtr()->trajectory_dimensions		= GC_acq_m1->getObjectPtr()->trajectory_dimensions;

		for (unsigned int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++) {
			GC_acq_m1_new->getObjectPtr()->user_float[i]			= GC_acq_m1->getObjectPtr()->user_float[i];
		}

		for (unsigned int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++) {
			GC_acq_m1_new->getObjectPtr()->user_int[i]				= GC_acq_m1->getObjectPtr()->user_int[i];
		}

		GC_acq_m1_new->getObjectPtr()->version						= GC_acq_m1->getObjectPtr()->version;

		return GADGET_OK;
	}

	inline int fCopyImageHeader(GadgetContainerMessage<ISMRMRD::ImageHeader> *tmp_m1, GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
	{
		tmp_m1->getObjectPtr()->average						= m1->getObjectPtr()->average;
		tmp_m1->getObjectPtr()->channels					= m1->getObjectPtr()->channels;
		tmp_m1->getObjectPtr()->contrast					= m1->getObjectPtr()->contrast;

		for (unsigned int i = 0; i < 3; i++) {
			tmp_m1->getObjectPtr()->field_of_view[i]		= m1->getObjectPtr()->field_of_view[i];
		}

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		tmp_m1->getObjectPtr()->data_type					= m1->getObjectPtr()->data_type;
#else
		tmp_m1->getObjectPtr()->image_data_type				= m1->getObjectPtr()->image_data_type;
#endif

		tmp_m1->getObjectPtr()->image_series_index			= m1->getObjectPtr()->image_series_index;
		tmp_m1->getObjectPtr()->image_index					= m1->getObjectPtr()->image_index;
		tmp_m1->getObjectPtr()->image_type					= m1->getObjectPtr()->image_type;

		for (unsigned int i = 0; i < 3; i++) {
			tmp_m1->getObjectPtr()->matrix_size[i]			= m1->getObjectPtr()->matrix_size[i];
		}

		tmp_m1->getObjectPtr()->phase						= m1->getObjectPtr()->phase;
		tmp_m1->getObjectPtr()->repetition					= m1->getObjectPtr()->repetition;
		tmp_m1->getObjectPtr()->set							= m1->getObjectPtr()->set;
		tmp_m1->getObjectPtr()->slice						= m1->getObjectPtr()->slice;
		tmp_m1->getObjectPtr()->acquisition_time_stamp		= m1->getObjectPtr()->acquisition_time_stamp;
		tmp_m1->getObjectPtr()->flags						= m1->getObjectPtr()->flags;
		tmp_m1->getObjectPtr()->measurement_uid				= m1->getObjectPtr()->measurement_uid;

		for (unsigned int i = 0; i < 3; i++) {
			tmp_m1->getObjectPtr()->patient_table_position[i]= m1->getObjectPtr()->patient_table_position[i];
		}

		for (unsigned int i = 0; i < 3; i++) {
			tmp_m1->getObjectPtr()->phase_dir[i]			= m1->getObjectPtr()->phase_dir[i];
		}

		for (unsigned int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++) {
			tmp_m1->getObjectPtr()->physiology_time_stamp[i]= m1->getObjectPtr()->physiology_time_stamp[i];
		}

		for (unsigned int i = 0; i < 3; i++) {
			tmp_m1->getObjectPtr()->position[i]				= m1->getObjectPtr()->position[i];
		}

		for (unsigned int i = 0; i < 3; i++) {
			tmp_m1->getObjectPtr()->read_dir[i]				= m1->getObjectPtr()->read_dir[i];
		}

		for (unsigned int i = 0; i < 3; i++) {
			tmp_m1->getObjectPtr()->slice_dir[i]			= m1->getObjectPtr()->slice_dir[i];
		}

		for (unsigned int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++) {
			tmp_m1->getObjectPtr()->user_float[i]			= m1->getObjectPtr()->user_float[i];
		}

		for (unsigned int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++) {
			tmp_m1->getObjectPtr()->user_int[i]				= m1->getObjectPtr()->user_int[i];
		}

		tmp_m1->getObjectPtr()->version						= m1->getObjectPtr()->version;

		return GADGET_OK;
	}

	template <typename T>
	inline int compare_complex_values(std::complex<T> c1, std::complex<T> c2)
	{
		T abs_c1 = std::abs(c1);
		T abs_c2 = std::abs(c2);

		if (abs_c1 > abs_c2) {
			// c1 is larger than c2 -> return +1
			return 1;
		} else if (abs_c1 < abs_c2) {
			// c1 is smaller than c2 -> return -1
			return -1;
		} else {
			// both magnitudes are equal, compare phases
			T arg_c1 = std::arg(c1);
			T arg_c2 = std::arg(c2);

			if (arg_c1 > arg_c2) {
				// c1 is larger than c2 -> return +1
				return 1;
			} else if (arg_c1 < arg_c2) {
				// c1 is smaller than c2 -> return -1
				return -1;
			} else {
				// both are equal -> return 0
				return 0;
			}
		}
	}

	template <typename T>
	void print_not_enough_ram_msg(const std::vector<T> &sizes, const int bytes_per_element)
	{
		// multiply all element values of dimkSpace_ vector
		long m = 1;
		for (const auto &e: sizes) {
			m *= e;
		}

		// calculate needed RAM in GiB
		int needed_ram = (m * bytes_per_element)/std::pow(1024, 3);
		GERROR("Not enough RAM to calculate! You need at least %dGiB of RAM.\n", needed_ram);
	}
}
