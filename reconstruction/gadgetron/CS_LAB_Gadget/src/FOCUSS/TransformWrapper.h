/*

file name	: 	TransformWrapper.h



author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)



version		: 	1.0



date		: 	03.01.2015



description	: 	interface for the base transformations. For a detailed description, it is suggested to read p. 57 of the thesis.

*/

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoMatrix_util.h"
#include "hoNDArray_utils.h"

#include "hoNDFFT_CS.h"
#include "hoNDDCT.h"
#include "hoNDKLT_CS.h"
#include "GlobalVar.h"

namespace Gadgetron{
	class TransformWrapper{

		public:
			virtual bool KernelFTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform) = 0;
			virtual bool KernelBTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform) = 0;
	};

	// interface for the DFT transformation
	class FFTWrapper : public TransformWrapper
	{
		public:
			// forward transformation
			virtual bool KernelFTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){
				hoNDFFT_CS<float>::instance()->ifft(&Array, (unsigned int)dim_to_transform);
				return GADGET_OK;
			};

			// backward transformation
			virtual bool KernelBTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){
				hoNDFFT_CS<float>::instance()->fft(&Array, (unsigned int)dim_to_transform);
				return GADGET_OK;
			};
	};

	// interface for the DCT transformation
	class DCTWrapper : public TransformWrapper
	{
		public:

			// forward transformation - DCT II
			virtual bool KernelFTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){
				// split real and imaginary part for split DCT
				hoNDArray<float> RealArray(Array.get_dimensions()), ImagArray(Array.get_dimensions());
				float *fRealPtr = RealArray.get_data_ptr(), *fImagPtr = ImagArray.get_data_ptr();
				std::complex<float> *cArrayPtr = Array.get_data_ptr();				

				for (int i = 0; i < RealArray.get_number_of_elements(); i++){
					fRealPtr[i] = cArrayPtr[i].real();
					fImagPtr[i] = cArrayPtr[i].imag();
				}

				// transform real part
				hoNDDCT<float>::instance()->dct(&RealArray, dim_to_transform);

				// transform imaginary part
				hoNDDCT<float>::instance()->dct(&ImagArray, dim_to_transform);

				// put real and imaginary parts into original array
				#pragma  omp parallel for
				for (int i = 0; i < RealArray.get_number_of_elements(); i++)
					cArrayPtr[i] = std::complex<float>(fRealPtr[i],fImagPtr[i]);

				return GADGET_OK;
			};

			// backward transformation - DCT III
			virtual bool KernelBTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){
				// split real and imaginary part for split DCT
				hoNDArray<float> RealArray(Array.get_dimensions()), ImagArray(Array.get_dimensions());
				float *fRealPtr = RealArray.get_data_ptr(), *fImagPtr = ImagArray.get_data_ptr();
				std::complex<float> *cArrayPtr = Array.get_data_ptr();

				for (int i = 0; i < RealArray.get_number_of_elements(); i++){
					fRealPtr[i] = cArrayPtr[i].real();
					fImagPtr[i] = cArrayPtr[i].imag();
				}

				// transform real part
				hoNDDCT<float>::instance()->idct(&RealArray, dim_to_transform);

				// transform imaginary part
				hoNDDCT<float>::instance()->idct(&ImagArray, dim_to_transform);

				// put real and imaginary parts into original array
				#pragma  omp parallel for
				for (int i = 0; i < RealArray.get_number_of_elements(); i++)
					cArrayPtr[i] = std::complex<float>(fRealPtr[i],fImagPtr[i]);

				return GADGET_OK;
			};
	};

	// interface for the PCA transformation
	class KLTWrapper : public TransformWrapper
	{
		public:			
			// forward transformation
			virtual bool KernelFTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){

				// check if transformation was prepared
				if (GlobalVar::instance()->vbStatPrinc_.at(dim_to_transform) == false){					
					GWARN("KLTWrapper: KLT transformation not prepared!\n");
					GlobalVar::instance()->KLTVec_.at(dim_to_transform)->prepare(Array, (size_t)dim_to_transform, (size_t)0, true);
				}

				// transform
				GlobalVar::instance()->KLTVec_.at(dim_to_transform)->ftransform(Array, Array, dim_to_transform);

				return GADGET_OK;
			};
			
			// backward transformation
			virtual bool KernelBTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){
				
				// check if transformation was prepared
				if (GlobalVar::instance()->vbStatPrinc_.at(dim_to_transform) == false){					
					GWARN("KLTWrapper: KLT transformation not prepared!\n");
				}

				// transform
				GlobalVar::instance()->KLTVec_.at(dim_to_transform)->btransform(Array, Array, dim_to_transform);

				return GADGET_OK;
			};
	};

	// interface for the Wavelet transformation - not further implemented in this version
	class WaveletWrapper : public TransformWrapper
	{
		public:
			virtual bool KernelFTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){
				GERROR("not implemented in this version...\n");
				return GADGET_OK;
			};

			virtual bool KernelBTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){
				GERROR("not implemented in this version...\n");
				return GADGET_OK;
			};
	};
}
