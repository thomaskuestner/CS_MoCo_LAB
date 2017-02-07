/*

file name	: 	TransformWrapper.h



author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)



version		: 	1.0



date		: 	03.01.2015



description	: 	interface for the base transformations. For a detailed description, it is suggested to read p. 57 of the thesis.

*/



#include "Gadget.h"

#include "GadgetMRIHeaders.h"

#include "hoNDFFT_CS.h"

#include "hoNDDCT.h"

//#include "hoNDPCA.h"

#include "GlobalVar_FOCUSS.h"

#include "hoMatrix_util.h"

#include "hoNDArray_utils.h"



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

	class PCAWrapper : public TransformWrapper

	{

		public:



			// forward transformation

			virtual bool KernelFTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){

				/*// vector for permutation

				std::vector<size_t> vtDimOrder;



				// get number of dims and get permutation vector

				int iNoDims = Array.get_number_of_dimensions();



				// 2D

				if (iNoDims == 3){

					vtDimOrder.push_back(0); vtDimOrder.push_back(1); vtDimOrder.push_back(2);

				}

				// 3D

				else if (iNoDims == 4){

					switch(dim_to_transform){

						case 0:

							// no permutation y-z-x-cha

							vtDimOrder.push_back(0); vtDimOrder.push_back(1); vtDimOrder.push_back(2); vtDimOrder.push_back(3);

							break;



						case 1: // permutation z-y-x-cha

							vtDimOrder.push_back(1); vtDimOrder.push_back(0); vtDimOrder.push_back(2); vtDimOrder.push_back(3);

							break;



						default:

							GADGET_DEBUG1("Not implemented in this version!\n");

							break;

					}

				}

				else

					GADGET_DEBUG1("Not implemented in this version - only 2D/3D support\n");



				// permute array

				hoNDArray<std::complex<float> >  TmpArray = *permute(&Array, &vtDimOrder, false);

				std::vector<size_t> dims = *TmpArray.get_dimensions();



				// check if principle components are stored for this dimension an ifno principle components found - initialize

				if (GlobalVar_FOCUSS::instance()->iStatPrinc_.at(dim_to_transform) == false){

					hoNDPCA::instance()->init_pc(TmpArray, dim_to_transform);

				}



				// transform

				hoNDPCA::instance()->ipca(*GlobalVar_FOCUSS::instance()->vPrincipleComponents_.at(dim_to_transform), TmpArray);

				*/

				return GADGET_OK;

			};



			// backward transformation

			virtual bool KernelBTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){

				/*// vector for permutation

				std::vector<size_t> vtDimOrder;



				// get number of dims and get permutation vector

				int iNoDims = Array.get_number_of_dimensions();



				// 2D

				if (iNoDims == 3){

					vtDimOrder.push_back(0); vtDimOrder.push_back(1); vtDimOrder.push_back(2);

				}

				// 3D

				else if (iNoDims == 4){

					switch(dim_to_transform){

						case 0:

							// no permutation y-z-x-cha

							vtDimOrder.push_back(0); vtDimOrder.push_back(1); vtDimOrder.push_back(2); vtDimOrder.push_back(3);

							break;



						case 1: // permutation z-y-x-cha

							vtDimOrder.push_back(1); vtDimOrder.push_back(0); vtDimOrder.push_back(2); vtDimOrder.push_back(3);

							break;



						default:

							GADGET_DEBUG1("Not implemented in this version!\n");

							break;

					}

				}

				else

					GADGET_DEBUG1("Not implemented in this version - only 2D/3D support\n");



				// permute array

				hoNDArray<std::complex<float> >  TmpArray = *permute(&Array, &vtDimOrder, false);

				std::vector<size_t> dims = *TmpArray.get_dimensions();



				// check if principle components are stored for this dimension an if no principle components found - initialize

				if (GlobalVar_FOCUSS::instance()->iStatPrinc_.at(dim_to_transform) == false){

					hoNDPCA::instance()->init_pc(TmpArray, dim_to_transform);

				}



				// transform

				hoNDPCA::instance()->pca(*GlobalVar_FOCUSS::instance()->vPrincipleComponents_.at(dim_to_transform), TmpArray);

				*/

				return GADGET_OK;

			};

	};



	// interface for the Wavelet transformation - not further implemented in this version

	class WaveletWrapper : public TransformWrapper

	{

		public:

			virtual bool KernelFTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){

				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1

					GDEBUG("not implemented in this version...\n");

				#else

					GADGET_DEBUG1("not implemented in this version...\n");

				#endif

				return GADGET_OK;

			};



			virtual bool KernelBTrafo(hoNDArray<std::complex<float> >  &Array, int dim_to_transform){

				#if __GADGETRON_VERSION_HIGHER_3_6__ == 1

					GDEBUG("not implemented in this version...\n");

				#else

					GADGET_DEBUG1("not implemented in this version...\n");

				#endif

				return GADGET_OK;

			};

	};

}
