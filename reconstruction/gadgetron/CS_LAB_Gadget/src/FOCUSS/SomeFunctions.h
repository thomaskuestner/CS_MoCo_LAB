/*
file name	: 	SomeFunctions.h

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	collection of mathematical operations on arrays

references	:	hoNDArray_math_util.h from the Gadgetron implementation
*/

#ifndef SOMEFUNCTIONS_H
#define SOMEFUNCTIONS_H

#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
	//#define GET_MACRO(_1,_2,_3,NAME,...) NAME
	//#define GDEBUG(...) GET_MACRO(__VA_ARGS__, GADGET_DEBUG1, GADGET_DEBUG2)(__VA_ARGS__)
	#define GADGET_DEBUG1(__VA_ARGS__) GDEBUG(__VA_ARGS__)
	#define GADGET_DEBUG2(x, ...) GDEBUG(x, ##__VA_ARGS__)
	#define GADGET_DEBUG_EXCEPTION(x,y) GEXCEPTION(x,y)
#endif

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

#if __GADGETRON_VERSION_HIGHER_3_6__ == 0
	#include "GadgetronCommon.h"
	#include <fstream>
	#include <ismrmrd_hdf5.h>
	#include <Shlwapi.h>
#endif
#include <complex>

#include "hoNDArray_utils.h"
#include "Gadget.h"

#include <ismrmrd.h>


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
bool fAbsPowDivide(hoNDArray<T> &a, float b, const hoNDArray<T> &c);

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
bool fAbsPow(hoNDArray<T> &a, float b);

// all zero?
bool fAllZero(const hoNDArray<bool> &Array);
bool fAllZero(const hoNDArray<std::complex<float> >  &Array);

// all non-zero?
bool fAllOne(const hoNDArray<bool> &Array);

// get hanning window values
std::vector<float>* fGetHanningWindow(int iElements);

// get hamming window values
std::vector<float>* fGetHammingWindow(int iElements);

// sum array in specified dimension and squeeze the result - compare to MATLAB sum(array, dim)
template <typename T>
bool sum_dim(hoNDArray<T> &Array, int dimension, hoNDArray<T> &result);

//inline int fCopyHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1, GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1_new);  

#if __GADGETRON_VERSION_HIGHER_3_6__ & WIN32 == 0
	inline bool save_array(hoNDArray< std::complex<float> > &Array, std::string file_prefix);
#endif

// flip array in specified dimension
template <typename T>
bool flip_array(hoNDArray<T> &Array, int dimension);

// flip line in array by specified offset
template <typename T>
bool flip_line(hoNDArray<T> &Array, size_t offset);

// sum array in specified dimension and squeeze the result - compare to MATLAB sum(array, dim)
template <typename T>
bool sum_dim(hoNDArray<T> &Array, int dimension, hoNDArray<T> &result);

// Gaussian 1D filter
template <typename T>
void filter1DGaussian(std::vector<T> &result, int length);

// vector convolution
template <typename T>
void vectorConv(std::vector<T> &vector, std::vector<T> &kernelVec, int option);

// array convolution - first dimension - same size
template <typename T>
void arrayConv(hoNDArray<T> &Array, std::vector<T> &kernelVec);

// 1D array convolution - same size - arbitrary dimension
template <typename T>
void arrayConv(hoNDArray<T> &Array, std::vector<T> &kernelVec, int dimensions);

// circular 1D vector shift
template <typename T>
void circshift(std::vector<T> &Array, int shift);

// circular 1D array shift (first dimension)
template <typename T>
void circshift(hoNDArray<T> &Array, int shift);

// circular 1D array shift (arbitrary dimension)
template <typename T>
void circshift(hoNDArray<T> &Array, int shift, int dimension);

// output a linear equally spaced vector
template <typename T>

std::vector<T>& linspace(T fStart, T fEnd, int iElements);

// interpolation
inline std::vector< float > interp1( std::vector< float > &x, std::vector< float > &y, std::vector< float > &x_new );//, int option);

inline int findNearestNeighbourIndex( float value, std::vector< float > &x );

// all elements true
inline bool allTrue(hoNDArray<bool> &Array);
inline bool allTrue(std::vector<bool> &Vector);

// get sub array
template <typename T>
void get_subarray(hoNDArray<T> &input, std::vector<size_t> vStart, std::vector<size_t> vSize, hoNDArray<T> &out);


template <typename T>
bool sum_dim_g(hoNDArray<T> &Array, int dimension);

inline int fCopyAcqHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1, GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1_new);

inline int fCopyImageHeader(GadgetContainerMessage<ISMRMRD::ImageHeader> *tmp_m1, GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);
}
#endif //SOMEFUNCTIONS_H

#include "SomeFunctions.hxx"
