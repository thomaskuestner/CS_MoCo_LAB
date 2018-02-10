/*
file name	: 	CS_FOCUSS.h
author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
version		: 	1.2
date		: 	17.12.2015
description	: 	abstract class for the member variables and member function prototypes. Base class for the dimension specific FOCUSS algorithm classes. For a detailed description, it is suggested to read ch. 4 of the thesis.
input		:	-
output		:	-
functions	:	process(...)		:	function prototype - reconstruct the Compressed Sensing k-space data with the FOCUSS algorithm and additional constraints for the Conjugate Gradient method (Total Variation, ESPReSSo)
				process_config(...)	:	read the flexible header information and XML configuration file - this function is the same for all sub-classes
				fInitVal(...)		:	map CS specific values (ESPReSSo, fully sampled, VDMap,..) between user values in the header fields and respective member variables - this function is the same for all sub-classes
				fGradESPReSSo(...)	:	function prototype - calculates the gradient of the ESPReSSo constraint
				fInitESPReSSo(...)	:	function prototype - calculates the symmetrical and conjugate sampling pattern as well as the 3D filter (Hanning/Hamming based) for the ESPReSSo reconstruction
				fWindowing(...)		:	function prototype - windows the incoming data set for the initial estimate for the FOCUSS reconstruction
				SetupTransformation(): creates and initializes the "transform"-object. The initialization is based on the XML configuration parameters.
				fSetupTransformationMatlab(...): see above..called from MATLAB
				fGetHanningWindow(...):calculates the Hann filter coefficients for the ESPReSSo filter
				fGetHammingWindow(...):calculates the Hamming filter coefficients for the ESPReSSo filter
				fRecon(...)			:	start FOCUSS reconstruction
variables	:	pbPtrN_				:	data pointer for boolean variables
				bESPReSSoIsLower_	:	upper or lower Partial Fourier data (true: data is lower Partial Fourier sampled)
				bControl_			:	flag, which indicates if the class is used as standalone version or in use with the "CS_Control"-class (the data set is put on the stream or passed to the CS_Control-class).
				habMaskConj_		:	boolean mask for the ESPReSSo algorithm
				habMaskConj2_		:	boolean mask for the ESPReSSo algorithm
				habMaskRight_		:	boolean mask for the ESPReSSo algorithm for pure CS data sets
				habMaskLeft_		:	boolean mask for the ESPReSSo algorithm for pure CS data sets
				habKSpaceCenter_	:	boolean mask, which stores region of the fully sampled k-space center
				iNChannels_			:	number of active channels
				iNOuter_			:	number of maximum FOCUSS iterations
				iNInner_			:	number of maximum CG iterations
				iESPReSSoDirection_	:	direction the Partial Fourier sampling (1: phase encoding direction, 2: partition encoding direction)
				vtDim_				:	dimension of the data object
				fP_					:	p-norm value (default: p = .5 for FOCUSS algorithm)
				fEpsilon_			:	tolerance - boundary for the CG algorithm
				fCSAccel_			:	acceleration value of the CS acquisitions
				pcfPtrN_			:	complex float data pointer
				hacfFilter_			:	3D filter for the ESPReSSo algorithm
				Transform_KernelTransform_: transformation object for the kernel transformation
				Transform_fftBA_	: 	transformation object, which controls the transformation before all (like the Fourier transformation in x-direction in the "normal" FOCUSS algorithm)
				Transform_fftAA_	:	transformation object, which controls the transformation after all, e. g. for outputting k-space data instead of image data
references	:	ESPReSSo: Küstner, T. et al. (2014):"ESPReSSo: A Compressed Sensing Partial k-Space Acquisition and Reconstruction"
				FOCUSS:
					- Algorithm				:	Gorodnitsky, I. and Rao, B. (1997): "Sparse Signal Reconstruction from Limited Data Using FOCUSS: A Re-weighted Minimum Norm Algorithm"
					- FOCUSS in MRI			:	Jung, H. et al. (2009): "k-t FOCUSS: A General Compressed Sensing Framework for High Resolution Dynamic MRI"
*/

#ifndef CS_FOCUSS_H
#define CS_FOCUSS_H

#pragma once
#include "CS_LAB_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_math_util.h"
#include "hoNDArray_operators.h"
#include "hoNDArray_blas.h"
#include "hoNDArray_utils.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include <complex>
#include <ctime>
#include <omp.h>
#include <ismrmrd.h>
#include <cmath>
#include "hoMatrix_util.h"

#include "SomeFunctions.h"
#include "GlobalVar.h"
#include "Transform.h"
#include "hoNDKLT_CS.h"

#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron
{
	// abstract base class for the FOCUSS reconstruction
	class EXPORTCSLAB CS_FOCUSS : public Gadget2< ISMRMRD::ImageHeader, hoNDArray< std::complex< float > > >
	{

	public:

		// reconstruct the k-space data in the process(...) method
		virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader >* m1, GadgetContainerMessage< hoNDArray< std::complex< float > > >* m2) = 0;

		// read the flexible data header
		int process_config(ACE_Message_Block* mb);//MS 17/03/30 =0;

		// read user specific values from header and initialize variables
		void fInitVal(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);

		// calculating gradient of ESPReSSo
		virtual void fGradESPReSSo(hoNDArray< std::complex< float > > & hacfRho, hoNDArray< std::complex< float > > &hacfFullMask, hoNDArray< std::complex< float > > &hacfKSpace, hoNDArray <std::complex< float > > &hacfW, hoNDArray< std::complex< float > > &hacfQ) = 0;

		// init filter array and sampling masks for ESPReSSo constraint
		virtual void fInitESPReSSo(hoNDArray< bool >& habFullMask) = 0;

		// windowing incoming data for initial estimate
		virtual void fWindowing(hoNDArray<std::complex< float > > & hacfWWindowed) = 0;

		// FOCUSS reconstruction
		virtual int fRecon(hoNDArray< std::complex< float > >  &hacfInput, hoNDArray< std::complex< float > >  &hacfRecon) = 0;

		// method for setting up the transformation objects
		void fSetupTransformation();

	// properties
	#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			
		// ESPReSSo active?
		GADGET_PROPERTY(CSESPReSSo, int, "CSESPReSSo", 0);

		// header config or xml config control
		GADGET_PROPERTY(bXMLControl, int, "XMLControl", 0);

		// residual of CG method
		GADGET_PROPERTY(iCGResidual, int, "CG Beta", 0);		

		// k-t FOCUSS loops
		GADGET_PROPERTY(OuterIterations, int, "OuterIterations", 2);
	
		// CG loops
		GADGET_PROPERTY(InnerIterations, int, "InnerIterations", 20);

		// FFT_Sparse dimension
		GADGET_PROPERTY(fftSparseDim, int, "FFT_Sparse", 0);

		// DCT_Sparse dimension
		GADGET_PROPERTY(dctSparseDim, int, "DCT_Sparse", 0);

		// PCA Sparse dimension
		GADGET_PROPERTY(pcaSparseDim, int, "PCA_Sparse", 0);

		// Scrambling dimension
		GADGET_PROPERTY(scrambleDim, int, "Scrambe_Dim", 0);

		// Kernel_FFT dimension
		GADGET_PROPERTY(kernelFftDim, int, "Kernel_FFT_dim", 0);

		// Transform_fftBA dimension
		GADGET_PROPERTY(transformFftBaDim, int, "Transform_fftBA_dim", 0);

		// kSpaceOut dimension
		GADGET_PROPERTY(kSpaceOutDim, int, "kSpaceOut", 0);

		// energy normalization
		GADGET_PROPERTY(norm, int, "norm", 0);

		// FOCUSS
		GADGET_PROPERTY(lambda, double, "lambda", 0.01);
		GADGET_PROPERTY(lambdaESPReSSo, double, "lambdaESPReSSo", 0.0);

	#endif


	// bool:
		// data pointer
		bool *pbPtr_, *pbPtr2_, *pbPtr3_;

		// data is upper or lower Partial Fourier data
		bool bESPReSSoIsLower_;

		// Control Flag - indicates if class is used as standalone Gadget or called from Control class
		bool bControl_;

		// second control flag - indicates if class parameters are set by XML config or by accu gadget
		int bXMLControl_;

	// hoNDArray<bool>:
		// Masks for ESPReSSo
		hoNDArray<bool> habMaskConj_, habMaskConj2_, habMaskRight_, habMaskLeft_, habKSpaceCenter_;

	// int:
		// residual of CG method
		//int iCGResidual_;

		// number of dimensions
		int iDim_;

		// number of channels
		int iNChannels_;

		// channel/k-space normalization
		int iNorm_;

	// vector int
		std::vector<int> viCalibrationSize_;

	// vector size_t:
		// dimension vector
		std::vector<size_t> vtDim_;

	// float:
		//l_p norm (default: p = .5)
		float fP_;

		// noise level (default: 1e-6)
		float fEpsilon_;

		// CS acceleration factor
		float fCSAccel_;

	// complex float
		// data pointer
		std::complex<float> *pcfPtr_, *pcfPtr2_, *pcfPtr3_;

	// hoNDArray<cx_float>:
		// FilterArray for ESPReSSo filtering
		hoNDArray<std::complex<float> >  hacfFilter_;

	// objects
		// Transformation object for KernelTransform
		Transform *Transform_KernelTransform_;

		// Transformation object for fftBA
		Transform *Transform_fftBA_;

		// Transformation object after all - put k-space or image on stream
		Transform *Transform_fftAA_;

		// debug output on/off
		bool bDebug_;

		// debug output - MATLAB (true) or Gadgetron (false)
		bool bMatlab_;

};

// inherited class for a 2D acquisition
class EXPORTCSLAB CS_FOCUSS_2D : public CS_FOCUSS
{
	public:
		CS_FOCUSS_2D(){	bControl_ = false;};

		GADGET_DECLARE(CS_FOCUSS_2D)

		int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

		//int process_config(ACE_Message_Block* mb);

		// 2D FOCUSS CS reconstruction
		int fRecon(hoNDArray<std::complex<float> >  &hacfInput, hoNDArray<std::complex<float> >  &hacfRecon);

	protected:
		// calculating gradient of ESPReSSo - not used in 2Dt
		void fGradESPReSSo(hoNDArray<std::complex<float> > & hacfRho, hoNDArray<std::complex<float> > &hacfFullMask, hoNDArray<std::complex<float> > &hacfKSpace, hoNDArray<std::complex<float> > &hacfW, hoNDArray<std::complex<float> > &hacfQ){};

		// init filter array and sampling masks for ESPReSSo constraint - not used in 2Dt
		void fInitESPReSSo(hoNDArray<bool>& habFullMask){};

		// windowing incoming data for initial estimate
		void fWindowing(hoNDArray<std::complex<float> > & hacfWWindowed);

		// get calibration size
		void fGetCalibrationSize(const hoNDArray<bool> &habArray);
};

//// inherited class for a 2Dt acquisition
//class EXPORTCSLAB CS_FOCUSS_2Dt : public CS_FOCUSS
//{
//	public:
//		CS_FOCUSS_2Dt(){ bControl_ = false; };
//
//		GADGET_DECLARE(CS_FOCUSS_2Dt)
//
//		int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
//
//		int process_config(ACE_Message_Block* mb);
//
//		// 2Dt FOCUSS CS reconstruction
//		int fRecon(hoNDArray<std::complex<float> >  &hacfInput, hoNDArray<std::complex<float> >  &hacfRecon);
//
//	protected:
//		// calculating gradient of ESPReSSo - not used in 2Dt
//		void fGradESPReSSo(hoNDArray<std::complex<float> > & hacfRho, hoNDArray<std::complex<float> > &hacfFullMask, hoNDArray<std::complex<float> > &hacfKSpace, hoNDArray<std::complex<float> > &hacfW, hoNDArray<std::complex<float> > &hacfQ){};
//
//		// init filter array and sampling masks for ESPReSSo constraint - not used in 2Dt
//		void fInitESPReSSo(hoNDArray<bool>& habFullMask){};
//
//		// windowing incoming data for initial estimate
//		void fWindowing(hoNDArray<std::complex<float> > & hacfWWindowed);
//
//		// get calibration size
//		void fGetCalibrationSize(const hoNDArray<bool> &habArray);
//};

// inherited class for a 3D acquisition
class EXPORTCSLAB CS_FOCUSS_3D : public CS_FOCUSS
{
    public:
		CS_FOCUSS_3D(){	bControl_ = false; };

		GADGET_DECLARE(CS_FOCUSS_3D)

		int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
		// 3D FOCUSS CS reconstruction
		int fRecon(hoNDArray<std::complex<float> >  &hacfInput, hoNDArray<std::complex<float> >  &hacfRecon);

		//int process_config(ACE_Message_Block* mb);

	protected:
		// calculating gradient of ESPReSSo
		void fGradESPReSSo(hoNDArray<std::complex<float> > & hacfRho, hoNDArray<std::complex<float> > &hacfFullMask, hoNDArray<std::complex<float> > &hacfKSpace, hoNDArray<std::complex<float> > &hacfW, hoNDArray<std::complex<float> > &hacfQ);

		// init filter array and sampling masks for ESPReSSo constraint
		void fInitESPReSSo(hoNDArray<bool>& habFullMask);

		// windowing incoming data for initial estimate
		void fWindowing(hoNDArray<std::complex<float> > & hacfWWindowed);

		// get calibration size
		void fGetCalibrationSize(const hoNDArray<bool> &habArray);
};

class EXPORTCSLAB CS_FOCUSS_4D : public CS_FOCUSS
{
    public:
		CS_FOCUSS_4D(){	bControl_ = false; };

		GADGET_DECLARE(CS_FOCUSS_4D)
		
		int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);
		// 4D FOCUSS CS reconstruction
		int fRecon(hoNDArray<std::complex<float> >  &hacfInput, hoNDArray<std::complex<float> >  &hacfRecon);		

		//int process_config(ACE_Message_Block* mb);

	protected:
		// calculating gradient of ESPReSSo
		void fGradESPReSSo(hoNDArray<std::complex<float> > & hacfRho, hoNDArray<std::complex<float> > &hacfFullMask, hoNDArray<std::complex<float> > &hacfKSpace, hoNDArray<std::complex<float> > &hacfW, hoNDArray<std::complex<float> > &hacfQ);

		// init filter array and sampling masks for ESPReSSo constraint
		void fInitESPReSSo(hoNDArray<bool>& habFullMask);

		// windowing incoming data for initial estimate
		void fWindowing(hoNDArray<std::complex<float> > & hacfWWindowed);

		// get calibration size
		void fGetCalibrationSize(hoNDArray<bool> &habArray);
};

}
#endif //CS_FOCUSS_H
