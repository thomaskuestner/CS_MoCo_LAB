/*	
file name	: 	hoNDDCT.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	implementation of the class "hoNDDCT" (file hoNDDCT.h)

references	:	- hoNDFFT.h from the original Gadgetron implementation
				- fftw.org
				
notes		:	main parts are from the original hoNDFFT.h file (original Gadgetron)
*/

#include "hoNDDCT.h"

namespace Gadgetron{

	// Singleton implementation
    template<typename T> hoNDDCT<T>* hoNDDCT<T>::instance()
    {
        if (!instance_) instance_ = new hoNDDCT<T>();
        return instance_;
    }

    template<class T> hoNDDCT<T>* hoNDDCT<T>::instance_ = NULL;

	// transformation in one specified dimension
    template<class T> void hoNDDCT<T>::dct_int(hoNDArray< T >* input, size_t dim_to_transform, int sign)
    {
		//+1: DCT III, -1: DCT II
        if (sign != -1 && sign != 1) return;
        if (dim_to_transform >= input->get_number_of_dimensions()) return;

        int stride     = 1;           //Distance between points in transform
        int dist       = 1;           //Distance between vectors
        int trafos     = 1;           //Transformations per chunk
        int chunks     = 1;           //Number of chunks
        int chunk_size = 1;           //Points per chunk
        int length     = 1;           //Length of each transform
        int total_dist = 1;

        T scale = 0.0;

        void* fft_plan        = 0;
        T*    fft_storage     = 0;

        T* fft_buffer = 0;
        T* data_ptr = 0;

        //Set sizes
        length = (int)input->get_size(dim_to_transform);
		
		// scale DCT III (inverse multiplied by logical DFT size - N = 2n)
        if (sign == 1)
        {
            scale = (T)(1.0/(2.0*length));
        }
		
		// scale DCT II
        else
        {
            scale = (T)1.0;
        }

        if (dim_to_transform != 0)
        {
            for (size_t i = 0; i < dim_to_transform; i++)
            {
                chunk_size *= (int)input->get_size(i);
            }
            stride = chunk_size;
            trafos = chunk_size;
            chunk_size *= length;

            for (size_t i = dim_to_transform+1; i < input->get_number_of_dimensions(); i++)
            {
                chunks *= (int)input->get_size(i);
            }
        }
        else
        {
            for (size_t i = 1; i < input->get_number_of_dimensions(); i++)
            {
                trafos *= (int)input->get_size(i);
            }
            chunk_size = trafos*length;

            dist = length;
        }
        total_dist = trafos*dist;

        //Allocate storage and make plan
        {
            mutex_.lock();
            fft_storage = (T*)fftw_malloc_pcfPtr_(sizeof(T)*length);
            if (fft_storage == 0)
            {
                std::cout << "Failed to allocate buffer for DCT" << std::endl;
                return;
            }
            fft_buffer = (T*)fft_storage;
			
			// planning flags - in-place transformation (input can be destroyed) and plan: MEASURE
            unsigned planner_flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;
			
			// decide between REDFT10: DCT II (forward transform) and REDFT01: DCT III (backward transform)
			if (sign == +1){
				fft_plan = fftw_plan_r2r_1d_pcfPtr_(length, fft_storage, fft_storage, FFTW_REDFT01, planner_flags);
			}
			else{
				fft_plan = fftw_plan_r2r_1d_pcfPtr_(length, fft_storage, fft_storage, FFTW_REDFT10, planner_flags);
			}

            if (fft_plan == 0)
            {
                fftw_free_pcfPtr_(fft_storage);
                std::cout << "Failed to create plan for DCT" << std::endl;
                return;
            }
            mutex_.unlock();
        }
	
        //Grab address of data
        data_ptr = reinterpret_cast<T*>(input->get_data_ptr());
		
        register int idx1_max = chunks*chunk_size;
        register int idx1, idx2;       //Index variables
        register int idx2_limit;
 
        for (idx1 = 0; idx1 < idx1_max; idx1+=chunk_size) //Loop over all chunks
        {
            idx2_limit = idx1+total_dist;
            for (idx2 = idx1; idx2 < idx2_limit; idx2+=dist) //Loop over all transformations
            {
                ///Copy data to buffer.
                {
					register int j, idx3 = idx2;
					for (j=0; j < length; idx3+=stride){
						 fft_buffer[j++] = data_ptr[idx3  ];
					}
                }

				// scale input and output to get same result like in Matlab - see DCT definition of Matlab and FFTW
				if (sign == 1)
				{
					for (int iI = 0; iI < length; iI++)
						fft_buffer[iI] *= 1/std::sqrt(0.5/(float)length);
					fft_buffer[0] *= std::sqrt(2.0);
				}

				// execute transformation
				fftw_execute_pcfPtr_(fft_plan);

				if (sign == -1)
				{
					for (int iI = 0; iI < length; iI++)
						fft_buffer[iI] *= std::sqrt(0.5/(float)length);
					fft_buffer[0] *= 1/std::sqrt(2.0);
				}

				
                {
                    register int j, idx3 = idx2;
					for (j=0; j < length; idx3+=stride){
						 data_ptr[idx3  ] = fft_buffer[j++]*scale;
					}
                }

            } //Loop over transformations
        } //Loop over chunks
	
        //clean up
        {
            mutex_.lock();
            if (fft_plan != 0)
            {
                fftw_destroy_plan_pcfPtr_(fft_plan);
            }

            if (fft_storage != 0)
            {
                fftw_free_pcfPtr_(fft_storage);
            }
            mutex_.unlock();
        }
    }

    template<> void hoNDDCT<float>::set_function_pointers()
    {
        fftw_import_wisdom_from_file_pcfPtr_ = &fftwf_import_wisdom_from_file;
        fftw_export_wisdom_to_file_pcfPtr_ = &fftwf_export_wisdom_to_file;
        fftw_cleanup_pcfPtr_ = &fftwf_cleanup;
        fftw_malloc_pcfPtr_ = &fftwf_malloc;
        fftw_free_pcfPtr_ = &fftwf_free;
        fftw_execute_pcfPtr_ = (void (*)(void*))(&fftwf_execute);
        fftw_plan_r2r_1d_pcfPtr_ = (void* (*)(int, void*, void*, fftw_r2r_kind, unsigned))(&fftwf_plan_r2r_1d);
        fftw_destroy_plan_pcfPtr_ = (void (*)(void*))(&fftwf_destroy_plan);
    }

    template<> void hoNDDCT<double>::set_function_pointers()
    {
        fftw_import_wisdom_from_file_pcfPtr_ = &fftw_import_wisdom_from_file;
        fftw_export_wisdom_to_file_pcfPtr_ = &fftw_export_wisdom_to_file;
        fftw_cleanup_pcfPtr_ = &fftw_cleanup;
        fftw_malloc_pcfPtr_ = &fftw_malloc;
        fftw_free_pcfPtr_ = &fftw_free;
        fftw_execute_pcfPtr_ = (void (*)(void*))(&fftw_execute);
        fftw_plan_r2r_1d_pcfPtr_ = (void* (*)(int, void*, void*, fftw_r2r_kind, unsigned))(&fftw_plan_r2r_1d);
        fftw_destroy_plan_pcfPtr_ = (void (*)(void*))(&fftw_destroy_plan);
    }
	
    template class hoNDDCT<float>;
    template class hoNDDCT<double>;
}
