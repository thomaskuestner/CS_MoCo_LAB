/*	
file name	: 	hoNDDCT.h

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	in-place DCT transformation based on the FFTW library

input		:	input				: 	data array, which will be transformed
				dim_to_transform	: 	direction of the transformation
				
output		:	input				: 	in-placed transformed array
				
functions	:	dct_int(...)		:	 DCT transformation
	
variables	:	instance_			:	Singleton object of the class
				
references	:	- hoNDFFT.h from the original Gadgetron implementation
				- fftw.org
*/

#ifndef hoNDDCT_H
#define hoNDDCT_H

#include "hoNDArray.h"
#include "hoNDArray_blas.h"
#include "hoMatrix.h"
#include "hoNDArray_math_util.h"
#include "hoNDArray_math.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include <boost/thread/mutex.hpp>
#include <iostream>
#include <fftw3.h>
#include <complex>

namespace Gadgetron
{
	/**
	Generic class for Fast Fourier Transforms using FFTW on the hoNDArray class.
	This class is a singleton because the planning and memory allocation routines of FFTW are NOT threadsafe.
	The class' template type is a REAL, ie. float or double.

	Access using e.g.
	FFT<float>::instance()
	*/
	template <typename T> class hoNDDCT
	{
	public:
		static hoNDDCT<T>* instance();

		// DCT II transformation in one specific dimension
		void dct(hoNDArray<T> *input, unsigned int dim_to_transform)
		{
			//-1 refers to the sign of the transform, -1 for FFTW_FORWARD
			dct_int(input,dim_to_transform,-1);
		}

		// DCT III transformation in one specific dimension
		void idct(hoNDArray<T> *input, unsigned int dim_to_transform)
		{
			//1 refers to the sign of the transform, +1 for FFTW_BACKWARD
			dct_int(input,dim_to_transform,1);
		}

		// DCT II transformation in all array dimensions
		void dct(hoNDArray< T >* input)
		{
			for (size_t i = 0; i < input->get_number_of_dimensions(); i++) {
				//-1 refers to the sign of the transform, -1 for FFTW_FORWARD
				dct_int(input,i,-1);
			}
		}

		// DCT III transformation in all array dimensions
		void idct(hoNDArray< T >* input)
		{
			for (size_t i = 0; i < input->get_number_of_dimensions(); i++) {
				//1 refers to the sign of the transform, +1 for FFTW_BACKWARD
				dct_int(input,i,1);
			}
		}

	protected:
		// class constructor
		hoNDDCT()
		{
			set_function_pointers();
		}

		// class destructor
		virtual ~hoNDDCT()
		{
			fftw_cleanup_pcfPtr_();
		}

		// DCT transformation in one specific dimension
		void dct_int(hoNDArray<T> *input, size_t dim_to_transform, int sign);

		// setting pointers
		void set_function_pointers();

		// pointers
		int   (*fftw_import_wisdom_from_file_pcfPtr_)(FILE*);
		void  (*fftw_export_wisdom_to_file_pcfPtr_)(FILE*);
		void  (*fftw_cleanup_pcfPtr_)(void);
		void* (*fftw_malloc_pcfPtr_)(size_t);
		void  (*fftw_free_pcfPtr_)(void* p);
		void  (*fftw_execute_pcfPtr_)(void*);
		void* (*fftw_plan_r2r_1d_pcfPtr_)(int, void*, void*, fftw_r2r_kind, unsigned);
		void  (*fftw_destroy_plan_pcfPtr_)(void*);

		// object instance
		static hoNDDCT<T>* instance_;
		boost::mutex mutex_;
	};

	template class hoNDDCT<float>;
	template class hoNDDCT<double>;

	// Declaration of specializations of set_function_pointers() function
	template<> void hoNDDCT<double>::set_function_pointers();
	template<> void hoNDDCT<float>::set_function_pointers();
}

#endif //hoNDDCT_H
