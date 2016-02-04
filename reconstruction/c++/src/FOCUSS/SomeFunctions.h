/*	
file name	: 	SomeFunctions.h

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	collection of mathematical operations on arrays
				
references	:	hoNDArray_math_util.h from the Gadgetron implementation
*/

#pragma once
#include "hoNDArray.h"
#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"
#include "ho6DArray.h"
#include "ho7DArray.h"
#include "hoNDImage.h"
#include "complext.h"
#include "cpucore_math_export.h"
#include "GadgetronCommon.h"
#include <complex>

namespace Gadgetron{
// r = a-b*c
template <typename T> 
bool fAminusBmultC(hoNDArray<T>& a, hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r);

// r = a*b-c
template <typename T> 
bool fAmultBminusC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r);
    
// r = a+b*c
template <typename T> 
bool fAplusBmultC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, hoNDArray<T>& r);
    
// r = a+b*c
template <typename T> 
bool fAplusBmultC(const hoNDArray<T>& a, const hoNDArray<T>& b, const hoNDArray<T>& c, const hoNDArray<T>& d, const hoNDArray<T>& e, const hoNDArray<T>& f, const hoNDArray<T>& g, const hoNDArray<T>& h, const hoNDArray<T>& i, hoNDArray<T>& r);

// r = a+b
template <typename T> 
bool fAplusB(T a, const hoNDArray<T>& b, hoNDArray<T>& r);

// calculate abs, pow to b and divide by signal energy
template <typename T>
bool fAbsPowDivide(hoNDArray<T> &a, const float b, const hoNDArray<T> &c);

// multiply: a = a.*b
template <typename T>
bool fMultiply(hoNDArray<T> &a, const hoNDArray<T> &b);

// G: r = -a*b+c*d+e*f+g*h+i*j 
template <typename T>
bool fCalcGradient(const hoNDArray<T>& a, const hoNDArray<T>& b, T c, const hoNDArray<T>& d, T e, const hoNDArray<T>& f, hoNDArray<T>& r);

//ESPReSSo: 0.5 * Lambda * (conj(W) .* q .* (1+phase) + conj(W).*W.*conj(q).*phase - fTrafkSpaceCombi.*(1+phase)); 
template <typename T>
bool fESPReSSoOut(const hoNDArray<T>& W, const hoNDArray<T>& q, const hoNDArray<T>& phase, const hoNDArray<T>& ifft_kSpaceCombi, hoNDArray<T>& out);

// crop array in y/z-plane
template <typename T>
bool fCropArrYZ(const hoNDArray<T> &Array, int a, int b, hoNDArray<T> &result);

// calculate signal energy of array
template <typename T>
float fCalcEnergy(hoNDArray<T> a);

// calculate abs and pow to b
template <typename T>
bool fAbsPow(hoNDArray<T> &a, const float b);

// all zero?
bool fAllZero(const hoNDArray<bool> &Array);
bool fAllZero(const hoNDArray<std::complex<float>> &Array);

// all non-zero?
bool fAllOne(const hoNDArray<bool> &Array);

// get hanning window values
std::vector<float>& fGetHanningWindow(int iElements);

// get hamming window values
std::vector<float>& fGetHammingWindow(int iElements);
}

