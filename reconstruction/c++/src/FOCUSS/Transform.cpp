/*	
file name	: 	Transform.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	implementation of the class "Transform" (file Transform.h)
*/

#include "Transform.h"

namespace Gadgetron{
	Transform::Transform() : bIsActive_(false){};

	// array transformtion - inplace transformation without FFT wrapping
	bool Transform::KernelFTransform(hoNDArray<std::complex<float>> &Array){
		// image domain (Cartesian base) -> new base
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_sparsity_.size(); i++){
			if (!(dynamic_cast<FFTWrapper*>(TVec_.at(i)))) // dynamic cast result NULL if not specified object type
				TVec_.at(i)->KernelFTrafo(Array, dims_to_trans_sparsity_.at(i));
		}	
		return GADGET_OK;
	}
	bool Transform::KernelBTransform(hoNDArray<std::complex<float>> &Array){
		// new base -> image domain (Cartesian base)
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_sparsity_.size(); i++){
			if (!(dynamic_cast<FFTWrapper*>(TVec_.at(i))))
				TVec_.at(i)->KernelBTrafo(Array, dims_to_trans_sparsity_.at(i));	
		}
		return GADGET_OK;
	}

	// array forward transformation - in-place transformation
	// depending on the stored DFT dimensions and sparsifying dimensions
	bool Transform::FTransform(hoNDArray<std::complex<float>> &Array){
		// k-space -> image domain
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_FFT_.size(); i++){
			hoNDFFT<float>::instance()->ifft(&Array, (unsigned int)dims_to_trans_FFT_.at(i));
		}
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_sparsity_.size(); i++){
			if (!(dynamic_cast<FFTWrapper*>(TVec_.at(i)))) // dynamic cast result NULL if not specified object type
				TVec_.at(i)->KernelFTrafo(Array, dims_to_trans_sparsity_.at(i));
		}	
		return GADGET_OK;
	}
	bool Transform::BTransform(hoNDArray<std::complex<float>> &Array){
		// new base -> image domain (Cartesian base)
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_sparsity_.size(); i++){
			if (!(dynamic_cast<FFTWrapper*>(TVec_.at(i))))
				TVec_.at(i)->KernelBTrafo(Array, dims_to_trans_sparsity_.at(i));	
		}
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_FFT_.size(); i++){	
			hoNDFFT<float>::instance()->fft(&Array, (unsigned int)dims_to_trans_FFT_.at(i));
		}
		return GADGET_OK;
	}

	// do F/BTransformation (sparsifying with fft wrapping) in only one direction - input: array, dim to transform, associated function for this dimension will be called
	bool Transform::FTransform(hoNDArray<std::complex<float>> &Array, int dim_to_transform){
		// iterate vector dims to find the entry
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_sparsity_.size(); i++){
			if(dims_to_trans_sparsity_.at(i) == dim_to_transform){
				hoNDFFT<float>::instance()->ifft(&Array, (unsigned int)dims_to_trans_sparsity_.at(i));
				
				// image domain (cartesian base) -> new base
				if (!(dynamic_cast<FFTWrapper*>(TVec_.at(i))))
					TVec_.at(i)->KernelFTrafo(Array, dims_to_trans_sparsity_.at(i));
			}
		}
		return GADGET_OK;
	}
	bool Transform::BTransform(hoNDArray<std::complex<float>> &Array, int dim_to_transform){
		// iterate vector dims to find the entry
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_sparsity_.size(); i++){
			if(dims_to_trans_sparsity_.at(i) == dim_to_transform){
			
				// new base -> image domain (cartesian base)
				if (!(dynamic_cast<FFTWrapper*>(TVec_.at(i))))
					TVec_.at(i)->KernelBTrafo(Array, dims_to_trans_sparsity_.at(i));
				hoNDFFT<float>::instance()->fft(&Array, (unsigned int)dims_to_trans_sparsity_.at(i));
			}
		}
		return GADGET_OK;
	}

	// do F/BTransformation (sparsifying with specified transformation type and dim to transform) - fft wrapper is checked if active for specified dimension
	bool Transform::FTransform(hoNDArray<std::complex<float>> &Array, int transformation_type, int dim_to_transform){
		std::function<bool(hoNDArray<std::complex<float>>&,int)> fPtr;
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_FFT_.size(); i++)
			if(dims_to_trans_FFT_.at(i) == dim_to_transform)
				hoNDFFT<float>::instance()->ifft(&Array, (unsigned int)dims_to_trans_FFT_.at(i));

		TransformWrapper* f;
		switch(transformation_type){
			case 0:
				//new FFTWrapper()->KernelFTrafo(Array, dim_to_transform);
				break;
			case 1:
				f = new DCTWrapper();
				f->KernelFTrafo(Array, dim_to_transform);
				break;
			case 2:
				f = new PCAWrapper();
				f->KernelFTrafo(Array, dim_to_transform);
				break;
			default:
				break;
		};
		return GADGET_OK;
	}
	bool Transform::BTransform(hoNDArray<std::complex<float>> &Array, int transformation_type, int dim_to_transform){
		std::function<bool(hoNDArray<std::complex<float>>&,int)> fPtr;
		TransformWrapper* f;
		switch(transformation_type){
			case 0:
				//fPtr = &FFTWrapper::KernelBTrafo;
				break;
			case 1:
				f = new DCTWrapper();
				f->KernelFTrafo(Array, dim_to_transform);
				break;
			case 2:
				f = new PCAWrapper();
				f->KernelFTrafo(Array, dim_to_transform);
				break;
			default:
				break;
		};
		for (std::vector<int>::size_type i = 0; i != dims_to_trans_FFT_.size(); i++)
			if(dims_to_trans_FFT_.at(i) == dim_to_transform)
					hoNDFFT<float>::instance()->fft(&Array, (unsigned int) dims_to_trans_FFT_.at(i));	
		
		return GADGET_OK;
	}

	// set active
	void Transform::set_active(){
		bIsActive_ = true;
	}

	// get active - true if transformation is active
	bool Transform::get_active(){
		return bIsActive_;
	}

	// append new transformation with associated transformation dimension
	void Transform::set_transformation_sparsity(int transformation_name, int dim_to_transform){

		// push dimension on dim_to_trans vector
		dims_to_trans_sparsity_.push_back(dim_to_transform);
		
		// push transformation on object pointer
		switch(transformation_name){
			case 0:
				TVec_.push_back(new FFTWrapper());
				break;

			case 1:
				TVec_.push_back(new DCTWrapper());
				break;

			case 2:
				TVec_.push_back(new PCAWrapper());
				break;
		}
	}

	// append new fft transform dimension
	void Transform::set_transformation_fft(int dim_to_transform){
		// push dimension on dim vec for FFT
		dims_to_trans_FFT_.push_back(dim_to_transform);
	}
}