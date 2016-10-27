/*
file name	: 	SomeFunctions.hxx

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	implementation of the prototypes in SomeFunctions.h

reference	:	hoNDArray_math_util.h from the Gadgetron implementation
*/

#include "SomeFunctions.h"
#include "Gadget.h"

namespace Gadgetron{

// r = a-b*c - all arrays same size
template <typename T>
bool fAminusBmultC(hoNDArray<T>& a, hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r)
{
	// check if number of elements are equal
	GADGET_DEBUG_CHECK_RETURN_FALSE(a.get_number_of_elements()==b.get_number_of_elements());
	GADGET_DEBUG_CHECK_RETURN_FALSE(b.get_number_of_elements()==c.get_number_of_elements());

	const T *pA = a.begin();
	const T *pB = b.begin();
	const T *pC = c.begin();
	T *pR = r.begin();

    try
    {
		if ( r.get_number_of_elements()!=a.get_number_of_elements())
        {
            r = a;
        }

		long long N = (long long)a.get_number_of_elements();
        long long n;

        #pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pR)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] = pA[n];
        }

		#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pB, pC, pR)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] -= pB[n]*pC[n];
        }
    }
    catch(...)
    {
				#ifdef __GADGETRON_VERSION_HIGHER_3_6__
					GERROR_STREAM("Error occurred in AminusBmultC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r) ... ");
				#else
					GADGET_ERROR_MSG("Error occurred in AminusBmultC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r) ... ");
				#endif
        return false;
    }

    return true;
}

// r = a*b-c
template <typename T>
bool fAmultBminusC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r)
{
	// check if number of elements are equal
	GADGET_DEBUG_CHECK_RETURN_FALSE(a.get_number_of_elements()==b.get_number_of_elements());
	GADGET_DEBUG_CHECK_RETURN_FALSE(b.get_number_of_elements()==c.get_number_of_elements());

	const T*pA = a.begin();
	const T*pB = b.begin();
	const T*pC = c.begin();
	T *pR = r.begin();

	long long N = a.get_number_of_elements();

    try
    {
		if ( r.get_number_of_elements()!=a.get_number_of_elements())
        {
            r.create(a.get_dimensions());
        }

        long long n;

        #pragma omp parallel for default(none)  schedule(static) private(n) shared(N, pA, pB, pR)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] = pA[n]*pB[n];
        }

		#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pC, pR)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] -= pC[n];
        }

    }
    catch(...)
    {
			#ifdef __GADGETRON_VERSION_HIGHER_3_6__
				GERROR_STREAM("Error occurred in AmultBminusC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r) ... ");
			#else
				GADGET_ERROR_MSG("Error occurred in AmultBminusC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r) ... ");
			#endif
        return false;
    }

    return true;
}

// r = a+b*c
template <typename T>
bool fAplusBmultC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r)
{
	// check if number of elements are equal
	GADGET_DEBUG_CHECK_RETURN_FALSE(a.get_number_of_elements()==b.get_number_of_elements());
	GADGET_DEBUG_CHECK_RETURN_FALSE(b.get_number_of_elements()==c.get_number_of_elements());

	const T*pA = a.begin();
	const T*pB = b.begin();
	const T*pC = c.begin();
	T *pR = r.begin();

	long long N = a.get_number_of_elements();

    try
    {
		if ( r.get_number_of_elements()!=a.get_number_of_elements())
        {
            r.create(a.get_dimensions());
        }

        long long n;

        #pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pB, pC, pR)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] = pA[n];
        }

		#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pB, pC, pR)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] += pB[n]*pC[n];
        }
    }
    catch(...)
    {
				#ifdef __GADGETRON_VERSION_HIGHER_3_6__
					GERROR_STREAM("Error occurred in AplusBmultC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r) ... ");
				#else
					GADGET_ERROR_MSG("Error occurred in AplusBmultC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r) ... ");
				#endif
        return false;
    }

    return true;
}

// r = a + b; a: constant value, b: array
template <typename T>
bool fAplusB(T a, const hoNDArray<T>& b, hoNDArray<T>& r)
{
	const T*pA = a.begin();
	const T*pB = b.begin();
	T *pR = r.begin();

	long long N = b.get_number_of_elements();

    try
    {
		if ( r.get_number_of_elements()!=b.get_number_of_elements())
        {
            r.create(b.get_dimensions());
        }

        long long n;

        #pragma omp parallel for default(none) schedule(static) private(n) shared(N, a, pB, pR)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] = a+pB[n];
        }

    }
    catch(...)
    {
			#ifdef __GADGETRON_VERSION_HIGHER_3_6__
				GERROR_STREAM("Error occurred in AplusB(T a, const hoNDArray<T>& b, hoNDArray<T>& r) ... ");
			#else
				GADGET_ERROR_MSG("Error occurred in AplusB(T a, const hoNDArray<T>& b, hoNDArray<T>& r) ... ");
			#endif
        return false;
    }

    return true;
}

// calculate abs and pow to b
template <typename T>
bool fAbsPow(hoNDArray<T> &a, float b){
	T*pA = a.begin();

	long long N = a.get_number_of_elements();

    try
    {
        long long n;

        #pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, b)
        for ( n=0; n<(long long)N; n++ )
        {
            pA[n] = std::pow(std::abs(pA[n]), b);
        }
    }
    catch(...)
    {
			#ifdef __GADGETRON_VERSION_HIGHER_3_6__
				GERROR_STREAM("Error occurred in pow_abs(hoNDArray<T> &a, const T b) ... ");
			#else
				GADGET_ERROR_MSG("Error occurred in pow_abs(hoNDArray<T> &a, const T b) ... ");
			#endif
        return false;
    }

	return true;
}

// calculate abs, pow to b and divide by signal energy
template <typename T>
bool fAbsPowDivide(hoNDArray<T> &a, float b, const hoNDArray<T> &c){
	T*pA = a.begin();
	const T* pC = c.begin();

	long long N = a.get_number_of_elements();

    try
    {
        long long n;

        #pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, b)
        for ( n=0; n<(long long)N; n++ )
        {
            pA[n] = std::pow(std::abs(pA[n]), b);
        }

		#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pC)
        for ( n=0; n<(long long)N; n++ )
        {
            pA[n] /= pC[n];
        }
    }
    catch(...)
    {
			#ifdef __GADGETRON_VERSION_HIGHER_3_6__
				GERROR_STREAM("Error occurred in pow_abs(hoNDArray<T> &a, const T b) ... ");
			#else
				GADGET_ERROR_MSG("Error occurred in pow_abs(hoNDArray<T> &a, const T b) ... ");
			#endif
        return false;
    }

	return true;
}

// multiply: a = a.*b
template <typename T>
bool fMultiply(hoNDArray<T> &a, const hoNDArray<T> &b){
	T*pA = a.begin();
	const T* pB = b.begin();

	long long N = a.get_number_of_elements();

    try
    {
        long long n;

        #pragma omp parallel for default(none)  schedule(static) private(n) shared(N, pA, pB)
        for ( n=0; n<(long long)N; n++ )
        {
            pA[n] *= pB[n];
        }
    }
    catch(...)
    {
				#ifdef __GADGETRON_VERSION_HIGHER_3_6__
					GERROR_STREAM("Error occurred in pow_abs(hoNDArray<T> &a, const T b) ... ");
				#else
					GADGET_ERROR_MSG("Error occurred in pow_abs(hoNDArray<T> &a, const T b) ... ");
				#endif
        return false;
    }

	return true;
}

// G: r = -conj(a)*b+c*d+e*f
template <typename T>
bool fCalcGradient(const hoNDArray<T>& a, const hoNDArray<T>& b, T c, const hoNDArray<T>& d, T e, const hoNDArray<T>& f, hoNDArray<T>& r)
{
	// check if number of elements are equal
	GADGET_DEBUG_CHECK_RETURN_FALSE(a.get_number_of_elements()==b.get_number_of_elements());
	GADGET_DEBUG_CHECK_RETURN_FALSE(b.get_number_of_elements()==d.get_number_of_elements());
	GADGET_DEBUG_CHECK_RETURN_FALSE(d.get_number_of_elements()==f.get_number_of_elements());

	const T*pA = a.begin();
	const T*pB = b.begin();
	const T*pD = d.begin();
	const T*pF = f.begin();
	T *pR = r.begin();

	long long N = b.get_number_of_elements();

    try
    {
		if ( r.get_number_of_elements()!=a.get_number_of_elements()){
            r.create(a.get_dimensions());

        }

        long long n;
		r=b;

		#pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pR)
        for ( n=0; n<(long long)N; n++ ){
            pR[n] *= -std::conj(pA[n]);
        }
		if (c != std::complex<float>(0.0)){
			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, c, pD, pR)
			for ( n=0; n<(long long)N; n++ ){
				pR[n] += c*pD[n];
			}
		}
		if (e != std::complex<float>(0.0)){
			#pragma omp parallel for default(none) schedule(static) private(n) shared(N, e, pF, pR)
			for ( n=0; n<(long long)N; n++ ){
				pR[n] += e*pF[n];
			}
		}
    }
    catch(...)
    {
			#ifdef __GADGETRON_VERSION_HIGHER_3_6__
				GERROR_STREAM("Error occurred in calc_gradient(...) .. ");
			#else
				GADGET_ERROR_MSG("Error occurred in calc_gradient(...) .. ");
			#endif
        return false;
    }

    return true;
}

//ESPReSSo: 0.5 * Lambda * (conj(W) .* q .* (1+phase) + conj(W).*W.*conj(q).*phase - fTrafkSpaceCombi.*(1+phase));
template <typename T>
bool fESPReSSoOut(const hoNDArray<T>& W, const hoNDArray<T>& q, const hoNDArray<T>& phase, const hoNDArray<T>& ifft_kSpaceCombi, hoNDArray<T>& out)
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

	try
	{
		if (out.get_number_of_elements() != W.get_number_of_elements()){
			out.create(W.get_dimensions());
		}

		long long n;

		#pragma omp parallel for default(none) private(n) shared(N, pW, pQ, pO)
		for (n = 0; n < N; n++){
			pO[n] = std::conj(pW[n])*pQ[n];
		}

		#pragma omp parallel for default(none) private(n) shared(N, pPh, pO, cVal)
		for (n = 0; n < N; n++){
			pO[n] *= (cVal + pPh[n]);
		}

		#pragma omp parallel for default(none) private(n) shared(N, pW, pQ, pPh, pO)
		for (n = 0; n < N; n++){
			pO[n] += std::conj(pW[n])*pW[n]*std::conj(pQ[n])*pPh[n];
		}

		#pragma omp parallel for default(none) private(n) shared(N, pPh, pkS, pO, cVal)
		for (n = 0; n < N; n++){
			pO[n] -= pkS[n]*(cVal + pPh[n]);
		}

		#pragma omp parallel for default(none) private(n) shared(N, pO, cVal2)
		for (n = 0; n < N; n++){
			pO[n] *= cVal2;
		}
	}
	catch(...)
	{
		#ifdef __GADGETRON_VERSION_HIGHER_3_6__
			GERROR_STREAM("Error occurred in ESPReSSoOut(const hoNDArray<T>& W, const hoNDArray<T>& q, const hoNDArray<T>& phase, const hoNDArray<T>& ifft_kSpaceCombi, hoNDArray<T>& out) ... ");
		#else
			GADGET_ERROR_MSG("Error occurred in ESPReSSoOut(const hoNDArray<T>& W, const hoNDArray<T>& q, const hoNDArray<T>& phase, const hoNDArray<T>& ifft_kSpaceCombi, hoNDArray<T>& out) ... ");
		#endif
		return false;
	}

	return true;
}

// crop center part with knowledge that x isn't PF - dim: y-z-x-c
template <typename T>
bool fCropArrYZ(const hoNDArray<T> &Array, int a, int b, hoNDArray<T> &result){
	int y = a, z = b;
	// get dims of incoming array
	std::vector<size_t> dims = *Array.get_dimensions();

	// get center indices
	int iCenterX = std::floor((float)(dims[2]-1)/2);
	int iCenterY = std::floor((float)(dims[0]-1)/2);
	int iCenterZ = std::floor((float)(dims[1]-1)/2);

	// get border indices
	//GADGET_DEBUG2("y: %i, z: %i\n", y,z);
	int iBorderY_Upper, iBorderY_Lower, iBorderZ_Upper, iBorderZ_Lower;
	if (y%2 == 0){
		iBorderY_Upper = iCenterY+std::ceil((float)y/2);
		iBorderY_Lower = iCenterY-std::floor((float)y/2)+1;
	}
	else{
		iBorderY_Upper = iCenterY+std::ceil((float)y/2)-1;
		iBorderY_Lower = iCenterY-std::floor((float)y/2);
	}
	if (z%2 == 0){
		iBorderZ_Upper = iCenterZ+std::ceil((float)z/2);
		iBorderZ_Lower = iCenterZ-std::floor((float)z/2)+1;
	}
	else{
		iBorderZ_Upper = iCenterZ+std::ceil((float)z/2)-1;
		iBorderZ_Lower = iCenterZ-std::floor((float)z/2);
	}

	//GADGET_DEBUG2("Border: y: %i,%i z: %i,%i\n", iBorderY_Upper, iBorderY_Lower, iBorderZ_Upper, iBorderZ_Lower);
	// indices out of bounds?
	if (iBorderY_Upper > dims[0]-1){
		iBorderY_Upper = dims[0]-1;
		y = dims[0];
	}
	if (iBorderZ_Upper > dims[1]-1){
		iBorderZ_Upper = dims[1]-1;
		z = dims[1];
	}
	if (iBorderY_Lower < 0)	iBorderY_Lower = 0;
	if (iBorderZ_Lower < 0)	iBorderZ_Lower = 0;

	// create output
	result.create(y,z);

	// fill output
	T* PtrIn = Array.get_data_ptr();
	T* PtrOut= result.get_data_ptr();
	for (int iy = iBorderY_Lower; iy <= iBorderY_Upper; iy++)
		for(int iz = iBorderZ_Lower; iz <= iBorderZ_Upper; iz++)
			PtrOut[(iy-iBorderY_Lower)+(iz-iBorderZ_Lower)*y] = PtrIn[iy + iz*dims[0] + iCenterX*dims[0]*dims[1]];

	return true;
}

// calculate signal energy of array
template <typename T>
float fCalcEnergy(hoNDArray<T> a){
	const T *pA = a.begin();
	hoNDArray<float> r(a.get_dimensions());

	float * pR = r.begin();
	float fResult;
    try
    {
		long long N = (long long)a.get_number_of_elements();
        long long n;

        #pragma omp parallel for default(none) schedule(static) private(n) shared(N, pA, pR)
        for ( n=0; n<(long long)N; n++ )
        {
            pR[n] = std::pow(abs(pA[n]),2);
        }

        for ( n=0; n<(long long)N; n++ )
        {
            fResult += pR[n];
        }
    }
    catch(...)
    {
			#ifdef __GADGETRON_VERSION_HIGHER_3_6__
				GERROR_STREAM("Error occurred in AminusBmultC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r) ... ");
			#else
				GADGET_ERROR_MSG("Error occurred in AminusBmultC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r) ... ");
			#endif
        return false;
    }

    return fResult;
}

// all zero?
inline bool fAllZero(const hoNDArray<bool> &Array){
	// output value
	bool isAllZero = true;
	
	bool *PtrIn = Array.get_data_ptr();
	for (long i = 0; i < Array.get_number_of_elements(); i++)
		if (PtrIn[i] != false){
			isAllZero = false;
			break;
		}

	return isAllZero;
}
inline bool fAllZero(const hoNDArray<std::complex<float> >  &Array){
	// output value
	bool isAllZero = true;
	
	std::complex<float> *PtrIn = Array.get_data_ptr();
	for (long i = 0; i < Array.get_number_of_elements(); i++)
		if (PtrIn[i] != std::complex<float>(0, 0)){
			isAllZero = false;
			break;
		}

	return isAllZero;
}

// function, which checks if all elements are one (bool)
inline bool fAllOne(const hoNDArray<bool> &Array){
	bool *PtrIn = Array.get_data_ptr();
	for (long i = 0; i < Array.get_number_of_elements(); i++)
		if (PtrIn[i] == false){		
			return false;
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
inline std::vector<float>& fGetHanningWindow(int iElements){
	const float pi = 3.141592653589793238463;
	std::vector<float> vfHanning(iElements);
	// L = N - 1
	float fL = float(iElements-1);
	// calc n
	for (std::vector<int>::size_type i = 0; i != vfHanning.size(); i++){
		vfHanning.at(i) = float(i)-fL/2;
	}
	// calc w(n)
	for (std::vector<int>::size_type i = 0; i != vfHanning.size(); i++){
		vfHanning.at(i) = .5*(1+std::cos(2*pi*vfHanning.at(i)/fL));
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
% w = hamming( L, 'symFlag') returns an L-point Hamming window using the window sampling
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
inline std::vector<float>& fGetHammingWindow(int iElements){
	const float pi = 3.141592653589793238463;
	std::vector<float> vfHamming(iElements);
	// L = N - periodic case
	float fL = float(iElements);
	// calc n
	for (std::vector<int>::size_type i = 0; i != vfHamming.size(); i++){
		vfHamming.at(i) = float(i)-fL/2;					
	}
	// calc w(n)
	for (std::vector<int>::size_type i = 0; i != vfHamming.size(); i++){
		vfHamming.at(i) = 0.54 + 0.46*std::cos(2*pi*vfHamming.at(i)/fL); 
	}

	return vfHamming;
}

// sum array in specified dimension and squeeze the result - compare to MATLAB sum(array, dim) 
template <typename T>
bool sum_dim(hoNDArray<T> &Array, int dimension, hoNDArray<T> &result){
	
	hoNDArray<T> NewTmpArray = Array;

	// get number of dimensions
	int iNoDim = NewTmpArray.get_number_of_dimensions();

	// get dimensions
	std::vector<size_t> vDims = *NewTmpArray.get_dimensions();

	// check dimensions
	if (Array.get_number_of_dimensions() <= dimension){
		BOOST_THROW_EXCEPTION( runtime_error("sum_dim: Error occurred - array is smaller than expected by the user!\n"));
	}
	else{
		if (iNoDim > 1){
			// permute array - summation dimension will be last dimension
			std::vector<size_t> vPermute;
			for (int i = 0; i < iNoDim; i++)
				if (i != dimension)
					vPermute.push_back(i);
		
			vPermute.push_back(dimension);
			NewTmpArray = *permute(&NewTmpArray, &vPermute,false);
	
			// number of elements of all dimension except the summation dimension
			long lMax = NewTmpArray.get_number_of_elements()/vDims.at(dimension);
			//GADGET_DEBUG2("lMax: %i\n", lMax);
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
}
}
