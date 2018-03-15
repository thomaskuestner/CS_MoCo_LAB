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

#include "gadgetron_messages.h"

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

#ifndef __GADGETRON_VERSION_HIGHER_3_6__
	#include "GadgetronCommon.h"
	#include <fstream>
	#include <ismrmrd_hdf5.h>
	#include <Shlwapi.h>
#endif
#include <complex>

#include "hoNDArray_utils.h"
#include "Gadget.h"

#include <ismrmrd.h>

#define ARRAYSIZE(a) sizeof(a)/sizeof(a[0])

namespace Gadgetron
{
	// r = a-b*c
	template <typename T>
	bool fAminusBmultC(hoNDArray<T> &a, hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r);

	// r = a*b-c
	template <typename T>
	bool fAmultBminusC(const hoNDArray<T> &a, const hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r);

	// r = a+b*c
	template <typename T>
	bool fAplusBmultC(const hoNDArray<T> &a, const hoNDArray<T> &b, const hoNDArray<T> &c, hoNDArray<T> &r);

	// r = a+b*c
	template <typename T>
	bool fAplusBmultC(const hoNDArray<T> &a, const hoNDArray<T> &b, const hoNDArray<T> &c, const hoNDArray<T> &d, const hoNDArray<T> &e, const hoNDArray<T> &f, const hoNDArray<T> &g, const hoNDArray<T> &h, const hoNDArray<T> &i, hoNDArray<T> &r);

	// r = a+b
	template <typename T>
	bool fAplusB(T a, const hoNDArray<T> &b, hoNDArray<T> &r);

	// calculate abs, pow to b and divide by signal energy
	template <typename T>
	bool fAbsPowDivide(hoNDArray<T> &a, float b, const hoNDArray<T> &c);

	// multiply: a = a.*b
	template <typename T>
	bool fMultiply(hoNDArray<T> &a, const hoNDArray<T> &b);

	// G: r = -a*b+c*d+e*f+g*h+i*j
	template <typename T>
	bool fCalcGradient(const hoNDArray<T> &a, const hoNDArray<T> &b, T c, const hoNDArray<T> &d, T e, const hoNDArray<T> &f, hoNDArray<T> &r);

	//ESPReSSo: 0.5 * Lambda * (conj(W) .* q .* (1+phase) + conj(W).*W.*conj(q).*phase - fTrafkSpaceCombi.*(1+phase));
	template <typename T>
	bool fESPReSSoOut(const hoNDArray<T> &W, const hoNDArray<T> &q, const hoNDArray<T> &phase, const hoNDArray<T> &ifft_kSpaceCombi, hoNDArray<T> &out);

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
	bool sum_dim(hoNDArray<T> &Array, unsigned int dimension, hoNDArray<T> &result);

	#if !defined(__GADGETRON_VERSION_HIGHER_3_6__) && WIN32 == 0
		inline bool save_array(hoNDArray<std::complex<float> > &Array, std::string file_prefix);
	#endif

	// flip array in specified dimension
	template <typename T>
	bool flip_array(hoNDArray<T> &Array, unsigned int dimension);

	// flip line in array by specified offset
	template <typename T>
	bool flip_line(hoNDArray<T> &Array, size_t offset);

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
	void arrayConv(hoNDArray<T> &Array, std::vector<T> &kernelVec, unsigned int dimensions);

	// circular 1D vector shift
	template <typename T>
	void circshift(std::vector<T> &Array, int shift);

	// circular 1D array shift (first dimension)
	template <typename T>
	void circshift(hoNDArray<T> &Array, int shift);

	// circular 1D array shift (arbitrary dimension)
	template <typename T>
	void circshift(hoNDArray<T> &Array, int shift, unsigned int dimension);

	// output a linear equally spaced vector
	template <typename T>
	std::vector<T>& linspace(T fStart, T fEnd, int iElements);

	// interpolation
	template <typename T>
	inline std::vector<T> interp1(std::vector<T> &x, std::vector<T> &y, std::vector<T> &x_new);//, int option);

	template <typename T>
	inline int findNearestNeighbourIndex(T value, std::vector<T> &x);

	// all elements true
	inline bool allTrue(hoNDArray<bool> &Array);
	inline bool allTrue(std::vector<bool> &Vector);

	// get sub array
	template <typename T>
	void get_subarray(hoNDArray<T> &input, std::vector<size_t> vStart, std::vector<size_t> vSize, hoNDArray<T> &out);

	template <typename T>
	bool sum_dim_g(hoNDArray<T> &Array, int dimension);

	inline int fCopyAcqHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_m1_new, const ISMRMRD::AcquisitionHeader *GC_acq_m1);

	inline int fCopyImageHeader(GadgetContainerMessage<ISMRMRD::ImageHeader> *tmp_m1, const ISMRMRD::ImageHeader *m1);

	/**
	* @brief Function to compare two complex values the MATLAB way.
	*
	* @param c1 complex value
	* @param c2 complex value
	* @return 1 if c1 > c2, -1 if c1 < c2 or 0 if c1 == c2
	*/
	template <typename T>
	inline int compare_complex_values(const std::complex<T> c1, const std::complex<T> c2);

	/**
	* @brief Prints message of too few RAM storage available.
	*
	* @param sizes a vector containing all element sizes to multiplicate
	* @param bytes_per_element amount of bytes each element is using
	*/
	template <typename T>
	void print_not_enough_ram_msg(const std::vector<T> &sizes, const int bytes_per_element);
}

#include "SomeFunctions.hxx"

#endif //SOMEFUNCTIONS_H
