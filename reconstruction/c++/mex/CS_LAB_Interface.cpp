//compile with: mex CS_FOCUSS.cpp -Iinclude\gadgetron -Iinclude\cuda -Isrc\FOCUSS -Iinclude\fftw -Iinclude -Llib\boost -Llib\ace -lACE -Llib\CS_LAB -lCS_LAB

/*==========================================================
 * CS_FOCUSS.cpp
 *
 * Author: Martin Schwartz
 *
 * Date: 18.12.2015
 *
 *========================================================*/
#define WIN32 1
#include "mex.h"
#include "matrix.h"

// include Gadgetron host based array type
#include <hoNDArray.h>

// include Gadgetron CS_LAB
#include "CS_LAB.h"

// declare Gadgetron namespace
using namespace Gadgetron;

// process Gadgetron array
void process(hoNDArray<std::complex<float>> &hoArray, const mxArray *matStruct)
{
    //hoNDArray<std::complex<float>> hoNewArray;
    double dFullySampled, dPartialFourierVal, dLambda, dLambdaESPReSSo, dNINNER, dNOUTER, dkSpaceOut, dDataset, dAlgorithm, dESPR_Active,dFFT_Sparse, dDCT_Sparse, dPCA_Sparse, dKernel_FFT_dim, dTransform_fftBA_dim;
    int iESPDirection;
    mxArray *mxPtr;
        
    /*********************************************************************/
    /************************* read Matlab struct ************************/
    /*********************************************************************/
    // get type of dataset
    mxPtr = mxGetField(matStruct, 0, "dataset");
    dDataset = *mxGetPr(mxPtr);
    mexPrintf("dDataset: %f\n", dDataset);mexEvalString("drawnow;");
    
    // get algorithm
    mxPtr = mxGetField(matStruct, 0, "algorithm");
    dAlgorithm = *mxGetPr(mxPtr);
    mexPrintf("dAlgorithm: %f\n", dAlgorithm);mexEvalString("drawnow;");
    
    // get fully sampled value
    mxPtr = mxGetField(matStruct, 0, "CSFullySampled");
    dFullySampled = *mxGetPr(mxPtr);
    mexPrintf("dFullySampled: %d\n", dFullySampled);mexEvalString("drawnow;");
    
    // get number of CG loops
    mxPtr = mxGetField(matStruct, 0, "iNINNER");
    dNINNER = *mxGetPr(mxPtr);
    mexPrintf("dNINNER: %f\n", dNINNER);mexEvalString("drawnow;");
    
    // get number of FOCUSS loops
    mxPtr = mxGetField(matStruct, 0, "iNOUTER");
    dNOUTER = *mxGetPr(mxPtr);
    mexPrintf("dNOUTER: %f\n", dNOUTER);mexEvalString("drawnow;");
    
    // get lambda - FOCUSS
    mxPtr = mxGetField(matStruct, 0, "lambda");
    dLambda = *mxGetPr(mxPtr);
    mexPrintf("dLambda: %f\n", dLambda);mexEvalString("drawnow;");
    
    // get lambda - ESPReSSo
    mxPtr = mxGetField(matStruct, 0, "lambdaESPReSSo");
    dLambdaESPReSSo = *mxGetPr(mxPtr);
    mexPrintf("dLambdaESPReSSo: %f\n", dLambdaESPReSSo);mexEvalString("drawnow;");
    
    // get partial fourier value
    mxPtr = mxGetField(matStruct, 0, "dPartialFourierVal");
    dPartialFourierVal = *mxGetPr(mxPtr);
    mexPrintf("dPartialFourierVal: %f\n", dPartialFourierVal);mexEvalString("drawnow;");
    
    // get ESPReSSo direction
    mxPtr = mxGetField(matStruct, 0, "ESPReSSoDirection");
    if (strcmp("y", (char*)mxGetPr(mxPtr)))
        iESPDirection = 1;
    else if (strcmp("z", (char*)mxGetPr(mxPtr)))
        iESPDirection = 2;
    else
        iESPDirection = 0;      
    mexPrintf("iESPDirection: %f\n", iESPDirection); mexEvalString("drawnow;"); 
    
    // ESPReSSo active
    mxPtr = mxGetField(matStruct, 0, "CS_ESPR");
    dESPR_Active = *mxGetPr(mxPtr);
    mexPrintf("dESPR_Active: %f\n", dESPR_Active);mexEvalString("drawnow;");
    
    // Transformations - FFT_Sparse
    mxPtr = mxGetField(matStruct, 0, "FFT_Sparse");
    dFFT_Sparse = *mxGetPr(mxPtr);
    mexPrintf("dFFT_Sparse: %f\n", dFFT_Sparse);mexEvalString("drawnow;");
    
    // Transformations - DCT_Sparse
    mxPtr = mxGetField(matStruct, 0, "DCT_Sparse");
    dDCT_Sparse = *mxGetPr(mxPtr);
    mexPrintf("dDCT_Sparse: %f\n", dDCT_Sparse);mexEvalString("drawnow;");
    
    // Transformations - PCA_Sparse
    mxPtr = mxGetField(matStruct, 0, "PCA_Sparse");
    dPCA_Sparse = *mxGetPr(mxPtr);
    mexPrintf("dPCA_Sparse: %f\n", dPCA_Sparse);mexEvalString("drawnow;");
    
    // Transformations - Kernel_FFT_dim
    mxPtr = mxGetField(matStruct, 0, "Kernel_FFT_dim");
    dKernel_FFT_dim = *mxGetPr(mxPtr);
    mexPrintf("dKernel_FFT_dim: %f\n", dKernel_FFT_dim);mexEvalString("drawnow;");
   
    // Transformations - Transform_fftBA_dim
    mxPtr = mxGetField(matStruct, 0, "Transform_fftBA_dim");
    dTransform_fftBA_dim = *mxGetPr(mxPtr);
    mexPrintf("dTransform_fftBA_dim: %f\n", dTransform_fftBA_dim);mexEvalString("drawnow;");
  
    // Transformations - kSpace out or image
    mxPtr = mxGetField(matStruct, 0, "kSpaceOut");
    dkSpaceOut = *mxGetPr(mxPtr);
    mexPrintf("dkSpaceOut: %f\n", dkSpaceOut);mexEvalString("drawnow;");
  
    // initialize CS_LAB algorithm parameters
    CS_LAB *pCS             = new CS_LAB();
    pCS->iNOuter_           = dNOUTER;
	pCS->iNInner_           = dNINNER;
	pCS->fP_ = .5;
    if (dESPR_Active != 0)
        pCS->bESPRActiveCS_ = true;
    else
        pCS->bESPRActiveCS_ = false;
    pCS->fEpsilon_          = (float)1e-6;	
    pCS->kSpaceOut_         = int(dkSpaceOut);
    pCS->iFFT_Sparse_       = int(dFFT_Sparse);
    pCS->iDCT_Sparse_       = int(dDCT_Sparse);
    pCS->iPCA_Sparse_       = int(dPCA_Sparse);
    pCS->iKernel_FFT_dim_   = int(dKernel_FFT_dim);
    pCS->iFFTBA_            = int(dTransform_fftBA_dim);
    pCS->bControl_          = false;
    pCS->iESPReSSoDirection_= 1;
    pCS->fFullySampled_     = (float)dFullySampled/100;
	pCS->fPartialFourierVal_= float(dPartialFourierVal);
    pCS->cfLambda_          = std::complex<float>(float(dLambda));
    pCS->cfLambdaESPReSSo_  = std::complex<float>(float(dLambdaESPReSSo));  
    pCS->bMatlab_           = true;                                         // enable debug output to MATLAB
    pCS->bDebug_            = true;                                         // enable debug output
    pCS->iDataset_          = int(dDataset);
    pCS->iAlgorithm_        = int(dAlgorithm);
    
    // initialize algorithm and start reconstruction
    pCS->fMatlabControl();
    int iRet = pCS->opCS_->fRecon(hoArray, hoArray);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    hoNDArray<std::complex<float>> hoInput;   // converted input array
    
    /* check for proper number of arguments */
    if(nrhs!=2) 
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    
    if(nlhs!=1) 
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    
    // first input argument has to be complex valued array
    if( !mxIsComplex(prhs[0])) 
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notComplex","Input matrix must be type complex.");
 
    if(!mxIsStruct(prhs[1]))
        mexErrMsgIdAndTxt( "MATLAB:phonebook:inputNotStruct","Input must be a structure.");
    
    /*********************************************************************/
    /********************* Matlab to Gadgetron array *********************/
    /*********************************************************************/
    mexPrintf("convert MATLAB to hoNDArray array..\n");mexEvalString("drawnow;");
    // get number of dimensions of Matlab array
    mwSize iNoDim = mxGetNumberOfDimensions(prhs[0]);
    mexPrintf("iNoDim: %i\n", iNoDim);mexEvalString("drawnow;");
    
    // get size of array
    const mwSize *iDim   = mxGetDimensions(prhs[0]);
    std::vector<size_t> vDim;
    for (int iI = 0; iI < iNoDim; iI++){
        vDim.push_back(iDim[iI]);
        mexPrintf("vDim[%i]: %i\n", iI, vDim[iI]);mexEvalString("drawnow;");
    }

    // get pointer to real and imaginary part of complex array
    float *pDataR = (float *)mxGetData(prhs[0]);
    float *pDataI = (float *)mxGetImagData(prhs[0]);

    // fill hoNDArray
    hoInput.create(vDim);                                               // create array
    std::complex<float> *pInput = hoInput.get_data_ptr();               // get pointer to array
    for (long iI = 0; iI < hoInput.get_number_of_elements(); iI++)      // loop over elements and copy element-wise
        pInput[iI] = std::complex<float>(pDataR[iI], pDataI[iI]);       // copy real and imaginary part in complex-valued array
    
    /*********************************************************************/
    /********************** process Gadgetron array **********************/
    /*********************************************************************/
    mexPrintf("process data..\n");mexEvalString("drawnow;");
    process(hoInput, prhs[1]);
    
    /*********************************************************************/
    /********************* Gadgetron to Matlab array *********************/
    /*********************************************************************/
    mexPrintf("convert hoNDArray to MATLAB array..");mexEvalString("drawnow;");
    // split array in real and imaginary part
    plhs[0] = mxCreateNumericArray(iNoDim, iDim, mxSINGLE_CLASS, mxCOMPLEX);
    float *pMATReal = (float *)mxGetPr(plhs[0]);
    float *pMATImag = (float *)mxGetPi(plhs[0]);    
    for (long iI = 0; iI < hoInput.get_number_of_elements(); iI++){
        pMATReal[iI] = pInput[iI].real();
        pMATImag[iI] = pInput[iI].imag();
    }
}