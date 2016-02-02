// ========================================================================
// ***
// *** fMeasCreateLUT.cpp
// ***
// *** Copyright 2014 Christian Wuerslin, University of Tuebingen and
// *** University of Stuttgart, Germany.
// *** Contact: christian.wuerslin@med.uni-tuebingen.de
// ***
// ========================================================================

#include "mex.h"
#include "matrix.h"
#include <io.h>
//#include <queue>
#include <cmath>
#include <fcntl.h>
#include <stdio.h>
//#include <sys/types.h>
//#include <sys/stat.h>

#define MAXSIZE 6000000
#define MDHLENGTH 16

typedef struct
{                                               //                  DrecksMDH 
    unsigned    int     ulFlagsAndDMALength;    // 4                NO
                int     lMeasUID;               // 8                NO
    unsigned    int     ulScanCounter;          // 12               NO
    unsigned    int     ulTimeStamp;            // 16               NO
    unsigned    int     ulPMUTimeStamp;         // 20               NO
    unsigned    int     aulEval[2];             // 28               NO
    unsigned    short   ushSamplesInScan;       // 30               NO 
    unsigned    short   ushUsedChannels;        // 32               NO 
    unsigned    short   aushLC[14];             // 60 = 32 + 14*2   YES
    unsigned    short   aushCutOff[2];          // 64               NO  
    unsigned    short   ushKSpaceCentre;        // 66               NO
    unsigned    short   ushCoilSelect;          // 68               NO
                float   fReadOutOffcenter;      // 72               YES
    unsigned    int     ulTimeSinceLastRF;      // 76               NO
    unsigned    short   ushKSpaceCentreLineNo;  // 78               YES
    unsigned    short   ushKSpaceCentrePartNo;  // 80               YES
    unsigned    short   aushIcePara[4];         // 88               YES
    unsigned    short   aushFreePara[4];        // 96               YES (Last is for float)
                float   afSliceData[7];         // 124 = 96 + 7*4   NO
    unsigned    short   ushChannelId;           // 126              NO
    unsigned    short   ushPTAPosNeg;           // 128              NO
}sMDH;

// ========================================================================
// ***
// *** MAIN MEX FUNCTION RegionGrowing_mex
// ***
// *** See m-file for description
// ***
// ========================================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // --------------------------------------------------------------------
    // Check the number of the input and output arguments.
    if(nrhs != 1) mexErrMsgTxt("Exactly one input argument required.");
    if(nlhs != 2) mexErrMsgTxt("Exactly two ouput arguments required.");
    // --------------------------------------------------------------------
    
    // --------------------------------------------------------------------
    // Get pointer/values to/of the input and outputs objects
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // 1st input: Filename
    if (!mxIsChar(prhs[0])) mexErrMsgTxt("First input argument must be of type char.");
    size_t buflen = mxGetN(prhs[0])*sizeof(mxChar)+1;
    char *pcPath = (char*) mxMalloc(buflen);
    mxGetString(prhs[0], pcPath, (mwSize)buflen);   
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Get pointer to output arguments and allocate memory for the corresponding objects
    const int pSize[2] = {38, MAXSIZE};
    plhs[0] = mxCreateNumericArray(2, pSize, mxDOUBLE_CLASS, mxREAL);	// create output array
    double *pdMDH = (double*) mxGetData(plhs[0]);						// get data pointer to mask
    
    const int pSize2[1] = {1};
    plhs[1] = mxCreateNumericArray(1, pSize2, mxDOUBLE_CLASS, mxREAL);	// create output array
    double *pdCount = (double*) mxGetData(plhs[1]);						// get data pointer to mask
    *pdCount = 0;
    
    FILE* pFile;
    pFile = fopen(pcPath, "rb");
    if(pFile == 0) mexErrMsgTxt("Invalid Filename");
    
    int iStatus, iSkip;
    sMDH mdh;
    
    iStatus = fread(&iSkip, sizeof(int), 1, pFile);
    if (iStatus == -1) {
        fclose(pFile);
        mexErrMsgTxt("Error!");
    }
    mexPrintf("Skipping %d bytes!\n", iSkip);
    fseek(pFile, iSkip, SEEK_SET);
    
    long long lMDHPos = static_cast<long long>(iSkip);
    
    float afDummy[8192];
    
    while (1) {
        
        //mexPrintf("Position is %u, ", lMDHPos);
    
        iStatus = fread(&mdh, sizeof(sMDH), 1, pFile);
        if (iStatus < 1) {
            mexPrintf("EOF reached after %f Iterations (fread)\n", *pdCount);
            fclose(pFile);
            return;
        }
        
        *pdMDH++ = (double) lMDHPos;
        *pdMDH++ = (double) mdh.aulEval[0];
        *pdMDH++ = (double) mdh.aulEval[1];
        *pdMDH++ = (double) mdh.ushSamplesInScan;
        *pdMDH++ = (double) mdh.ushUsedChannels;
        *pdMDH++ = (double) mdh.aushLC[0];
        *pdMDH++ = (double) mdh.aushLC[1];
        *pdMDH++ = (double) mdh.aushLC[2];
        *pdMDH++ = (double) mdh.aushLC[3];
        *pdMDH++ = (double) mdh.aushLC[4];
        *pdMDH++ = (double) mdh.aushLC[5];
        *pdMDH++ = (double) mdh.aushLC[6];
        *pdMDH++ = (double) mdh.aushLC[7];
        *pdMDH++ = (double) mdh.aushLC[8];
        *pdMDH++ = (double) mdh.aushLC[9];
        *pdMDH++ = (double) mdh.aushLC[10];
        *pdMDH++ = (double) mdh.aushLC[11];
        *pdMDH++ = (double) mdh.aushLC[12];
        *pdMDH++ = (double) mdh.aushLC[13];
        *pdMDH++ = (double) mdh.ushKSpaceCentreLineNo;
        *pdMDH++ = (double) mdh.ushKSpaceCentrePartNo;
        *pdMDH++ = (double) mdh.aushIcePara[0];
        *pdMDH++ = (double) mdh.aushIcePara[1];
        *pdMDH++ = (double) mdh.aushIcePara[2];
        *pdMDH++ = (double) mdh.aushIcePara[3];
        *pdMDH++ = (double) mdh.aushFreePara[0];
        *pdMDH++ = (double) mdh.aushFreePara[1];
        *pdMDH++ = (double) mdh.aushFreePara[2];
        *pdMDH++ = (double) mdh.aushFreePara[3];
        *pdMDH++ = (double) mdh.fReadOutOffcenter;
        *pdMDH++ = (double) mdh.afSliceData[0];
        *pdMDH++ = (double) mdh.afSliceData[1];
        *pdMDH++ = (double) mdh.afSliceData[2];
        *pdMDH++ = (double) mdh.afSliceData[3];
        *pdMDH++ = (double) mdh.afSliceData[4];
        *pdMDH++ = (double) mdh.afSliceData[5];
        *pdMDH++ = (double) mdh.afSliceData[6];
        *pdMDH++ = (double) mdh.ushPTAPosNeg;
        //mexPrintf("Samples: %u\n", mdh.ushSamplesInScan);
        //iStatus = fseek(pFile, long(mdh.ushSamplesInScan*8), SEEK_CUR);
        iStatus = fread(afDummy, sizeof(float), 2*mdh.ushSamplesInScan, pFile);
        if (iStatus < 2*mdh.ushSamplesInScan) {
            mexPrintf("EOF reached after %10.0f Iterations (fseek)\n", *pdCount);
            fclose(pFile);
            return;
        }
        lMDHPos += static_cast<long long>(2*sizeof(float)*mdh.ushSamplesInScan + sizeof(sMDH));
        (*pdCount)++;
    }
    fclose(pFile);
    mexPrintf("Terminating normally\n");
}
// ========================================================================
// *** END OF MAIN MEX FUNCTION fMeasCreateLUT
// ========================================================================