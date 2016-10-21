/*	
file name	: 	Transform.h

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	wrapper for several kernel transformations

input		:	Array				: 	data array, which will be transformed
				dim_to_transform	: 	direction of the transformation
				transformation_type	:	kind of transformation (DFT, DCT, PCA,..)
				
output		:	Array				: 	in-placed transformed array
				
functions	:	KernelB/FTransform	:	kernel backward or forward transformation
				B/FTransform(...)	:	backward and forward transformation in all/one dimension with stored or specified transformation
				set_active(...)		:	set transformation active
				get_active(...)		:	get boolean value if transformation is active or not
				set_transformation_sparsity(...)	:	append new sparsifying transformation in specified direction
				set_transformation_fft(...)			:	append new DFT transformation in specified direction
	
variables	:	dims_to_trans_sparsity_				:	vector, which stores the dimensions to be transformed with the sparsifying transformation
				dims_to_trans_FFT_					: 	vector, which stores the dimensions to be transformed with the DFT transformation
				bIsActive_							:	flag indicates whether the transformation is active or not
				TVec_								:	vector, which stores the pointer to the transformation objects				
*/

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "TransformWrapper.h"
#include "hoNDFFT.h"
#include "CS_LAB_export.h"

namespace Gadgetron{

	class EXPORTCSLAB Transform{

		public:
			// constructor
			Transform();

			// array forward transformtion - inplace transformation without FFT wrapping
			bool KernelFTransform(hoNDArray<std::complex<float> >  &Array);

			// array backward transformtion - inplace transformation without FFT wrapping
			bool KernelBTransform(hoNDArray<std::complex<float> >  &Array);

			// array forward transformation - inplace transformation
			bool FTransform(hoNDArray<std::complex<float> >  &Array);

			// array backward transformation - inplace transformation
			bool BTransform(hoNDArray<std::complex<float> >  &Array);

			// do transformation in only one direction - input: dim to transform, associated function for this dimension will be called
			bool FTransform(hoNDArray<std::complex<float> >  &Array, int dim_to_transform);
			bool BTransform(hoNDArray<std::complex<float> >  &Array, int dim_to_transform);

			// do transformation on array, input: array, type of transformation, dim_to_transform
			bool FTransform(hoNDArray<std::complex<float> >  &Array, int dim_to_transform, int transformation_type);
			bool BTransform(hoNDArray<std::complex<float> >  &Array, int dim_to_transform, int transformation_type);

			// do transformation in specified dimensions - input: array, vector with dimensions
			bool FTransform(hoNDArray<std::complex<float> >  &Array, std::vector<int> &dims);
			bool BTransform(hoNDArray<std::complex<float> >  &Array, std::vector<int> &dims);

			// do transformation in specified dimensions and transformation type - input: array, type, vector with dimensions
			bool FTransform(hoNDArray<std::complex<float> >  &Array, std::vector<int> &dims, int transformation_type);
			bool BTransform(hoNDArray<std::complex<float> >  &Array, std::vector<int> &dims, int transformation_type);
			
			// set active
			void set_active();

			// get active - true if transformation is active
			bool get_active();

			// append new sparsifying transform with associated transformation dimension
			void set_transformation_sparsity(int transformation_name, int dim_to_transform);

			// append new fft transform dimension
			void set_transformation_fft(int dim_to_transform);

	private:

		// vector containing the dimensions to transform (sparsity)
		std::vector<int> dims_to_trans_sparsity_;

		// vector containing the dimensions to transform (FFT)
		std::vector<int> dims_to_trans_FFT_;
	
		// activation flag
		bool bIsActive_;

		// vector for transformation objects
		std::vector<TransformWrapper*> TVec_;
	};
}