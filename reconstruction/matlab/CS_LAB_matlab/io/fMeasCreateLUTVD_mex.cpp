// ========================================================================
// ***
// *** fMeasCreateLUTVD.cpp
// ***
// *** Copyright 2018 Thomas Kuestner, University of Tuebingen and
// *** University of Stuttgart, Germany.
// *** Contact: thomas.kuestner@med.uni-tuebingen.de
// ***
// ========================================================================

#include "mex.h"
#include "matrix.h"
#include <io64.h>
//#include <queue>
#include <cmath>
#include <fcntl.h>
#include <stdio.h>
#include <cstdint>
//#include <sys/types.h>
#include <sys/stat.h>

#define MAXSIZE 6000000
#define MDHLENGTH 16
#define _LARGEFILE_SOURCE

typedef struct
{                                               //                    bytes     DrecksMDH 
    unsigned    int     ulFlagsAndDMALength;    // 1:4                  4           NO
                int     lMeasUID;               // 5:8                  4           NO
    unsigned    int     ulScanCounter;          // 9:12                 4           NO
    unsigned    int     ulTimeStamp;            // 13:16                4           NO
    unsigned    int     ulPMUTimeStamp;         // 17:20                4           NO
    // skip 20 bytes: system type, PTAB and reserved
    unsigned    short   ushSystemType;          // 21:22                2           NO
    unsigned    short   ushPTABPosDelay;        // 23:24                2           NO
                int     lPTABPosX;              // 24:27                4           NO
                int     lPTABPosY;              // 28:31                4           NO
                int     lPTABPosZ;              // 32:35                4           NO
    unsigned    int     ulReserved;             // 36:39                4           NO
    
    unsigned    int     aulEval[2];             // 40:47                8           NO
    unsigned    short   ushSamplesInScan;       // 48:49                2           NO 
    unsigned    short   ushUsedChannels;        // 50:51                2           NO 
    unsigned    short   aushLC[14];             // 52:79 = 52+14*2-1    28          YES
    unsigned    short   aushCutOff[2];          // 80:83                4           NO  
    unsigned    short   ushKSpaceCentre;        // 84:85                2           NO
    unsigned    short   ushCoilSelect;          // 86:87                2           NO
                float   fReadOutOffcenter;      // 88:91                4           YES
    unsigned    int     ulTimeSinceLastRF;      // 92:95                4           NO
    unsigned    short   ushKSpaceCentreLineNo;  // 96:97                2           YES
    unsigned    short   ushKSpaceCentrePartNo;  // 98:99                2           YES
				float   afSliceData[7];         // 100:127 = 100+7*4-1  28          NO
    unsigned    short   aushIcePara[24];        // 128:175 = 128+24*2-1 48          YES (only first 4 parameter = 8 byte)
    unsigned    short   aushFreePara[4];        // 176:183              8           YES (Last is for float)
    
    unsigned    short   ushApplicationCounter;  // 184:185              2           NO
    unsigned    short   ushApplicationMask;     // 186:187              2           NO
    unsigned    int     ulCRC;                  // 188:191              4           NO
                                                //                      = 192 bytes                    
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
    mexPrintf("Loading data: %s\n", pcPath);
    pFile = fopen(pcPath, "rb");
    if(pFile == 0) mexErrMsgTxt("Invalid Filename");
    
    int iStatus, iNScans, iSkip, iHeader;
    uint64_T uiSkip;
    int64_T iOffset = 0;
    int64_T iPosition = 0;
    sMDH mdh;
    //fseek(pFile, 0, SEEK_END);
    //int iFilesize = ftell(pFile);
    //struct __stat64 sFilestatus;
    //__stat64(pcPath, &sFilestatus);
    //long long lFilesize = static_cast<int>(sFilestatus.st_size);
    //mexPrintf("File size = %.0f bytes\n", double(lFilesize));
    //fseek(pFile, 0, SEEK_SET); // rewind to beginning of file
    structStat statbuf;
    int64_T iFileSize = 0;
    if(getFileFstat(fileno(pFile), &statbuf)==0)
    {
        iFileSize = statbuf.st_size;
        mexPrintf("File size = %" FMT64 "d bytes\n", iFileSize);
    }    
    
    iStatus = fread(&iSkip, sizeof(int), 1, pFile);
    if (iStatus == -1) {
        fclose(pFile);
        mexErrMsgTxt("Error!");
    }
    mexPrintf("First int = %d\n", iSkip);
    iStatus = fread(&iSkip, sizeof(int), 1, pFile);
    mexPrintf("Second int = %d\n", iSkip);
    iNScans = iSkip;
    uint64_T* measOffset = new uint64_T[iNScans];
    uint64_T* measLength = new uint64_T[iNScans];
    
    iStatus = fread(&iSkip, sizeof(int), 1, pFile);
    mexPrintf("Meas ID = %d\n", iSkip);
    iStatus = fread(&iSkip, sizeof(int), 1, pFile);
    mexPrintf("File ID = %d\n", iSkip);
    
    long long lMDHPos;
    //float afDummy[8192];
    int iDMASkip = 32;
    int i;
    unsigned int uiDMA = 0;
    
    for(i=0; i<iNScans; i++)
    {  
        iStatus = fread(&uiSkip, sizeof(uint64_T), 1, pFile);
        measOffset[i] = uiSkip;
        iStatus = fread(&uiSkip, sizeof(uint64_T), 1, pFile);
        measLength[i] = uiSkip;
        getFilePos(pFile, (fpos_T*) &iPosition);
        iOffset = iPosition + 152-16;
        setFilePos(pFile, (fpos_T*) &iOffset);
        //iStatus = fseek(pFile, 152 - 16, SEEK_CUR);
    }
    
    
    for(i=0; i<iNScans; i++)
    {    
        //mexPrintf("Only reading 2nd measurement!\n");
        mexPrintf("--- Scan %d ---\n", i);
        mexPrintf("Skipping %d bytes (meas offset %d)!\n", measOffset[i], i);
        //fseek(pFile, measOffset[i], SEEK_SET);
        setFilePos(pFile, (fpos_T*) &(measOffset[i]));
        iStatus = fread(&iHeader, sizeof(int), 1, pFile);
        mexPrintf("Header length = %d\n", iHeader);
        //fseek(pFile, measOffset[i]+iHeader, SEEK_SET);
        iOffset = measOffset[i] + iHeader;
        setFilePos(pFile, (fpos_T*) &iOffset);
        mexPrintf("Pointer starting at pos = %d\n", measOffset[i]+iHeader);
    
        lMDHPos = static_cast<long long>(measOffset[i]+iHeader);
        uiDMA = 0;

        while (1) {

            //mexPrintf("Position is MDHPos %u\n", lMDHPos);
            //mexPrintf("Position is file %u (before read)\n", ftell(pFile));

            iStatus = fread(&mdh, sizeof(sMDH), 1, pFile);
            if (iStatus < 1) {
                mexPrintf("EOF reached after %f Iterations (fread)\n", *pdCount);
                delete [] measOffset;
                delete [] measLength;
                fclose(pFile);
                return;
            }
            //mexPrintf("Position is file %u (after read)\n", ftell(pFile));

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
            *pdMDH++ = (double) mdh.ushPTABPosDelay;
            //mexPrintf("Samples: %u\n", mdh.ushSamplesInScan);
            //iStatus = fseek(pFile, long(mdh.ushSamplesInScan*8), SEEK_CUR);
            
            //mexPrintf("evalInfoMask[0] = %u\n", mdh.aulEval[0]);
            //mexPrintf("evalInfoMask[1] = %u\n", mdh.aulEval[1]);
            if((mdh.aulEval[0] & 1)) // (~mdh.ulFlagsAndDMALength & 0xFFFFFF00) // MDH_ACQEND | uiDMA == 0
            {
                // correct by 8 byte
                //iStatus = fseek(pFile, -8, SEEK_CUR);
                getFilePos(pFile, (fpos_T*) &iPosition);
                iOffset = iPosition - 8;
                setFilePos(pFile, (fpos_T*) &iOffset);
                lMDHPos -= 8;
                uiDMA = mdh.ulFlagsAndDMALength;
                //iStatus = fseek(pFile, uiDMA - 184, SEEK_CUR);
                getFilePos(pFile, (fpos_T*) &iPosition);
                iOffset = iPosition + uiDMA - 184;
                setFilePos(pFile, (fpos_T*) &iOffset);
                getFilePos(pFile, (fpos_T*) &iPosition);
                mexPrintf("MDH_ACQEND: Position is file %" FMT64 "d \n", iPosition);
                //mexPrintf("MDH_ACQEND: Position is file %u \n", ftell(pFile));
                lMDHPos += static_cast<long long>(uiDMA - 184);
                if(lMDHPos % 512) // jump to next full 512 bytes
                {
                    lMDHPos += static_cast<long long>(512 - (lMDHPos % 512) + sizeof(sMDH));
                }
                break;
            }
            if(mdh.aulEval[0] & 32) // MDH_SYNCDATA 
            {
                // correct by 8 byte
                //iStatus = fseek(pFile, -8, SEEK_CUR);
                getFilePos(pFile, (fpos_T*) &iPosition);
                iOffset = iPosition - 8;
                setFilePos(pFile, (fpos_T*) &iOffset);
                lMDHPos -= 8;
                //mexPrintf("MDH_SYNCDATA: Position is file %u (in)\n", ftell(pFile));
                uiDMA = mdh.ulFlagsAndDMALength; // & 0xFFFFFF01;
                //mexPrintf("mdh.ulFlagsAndDMALength = %u\n", mdh.ulFlagsAndDMALength);
                //mexPrintf("uiDMA = %u\n", uiDMA);
                //iStatus = fseek(pFile, uiDMA - 184, SEEK_CUR);
                getFilePos(pFile, (fpos_T*) &iPosition);
                iOffset = iPosition + uiDMA - 184;
                setFilePos(pFile, (fpos_T*) &iOffset);
                getFilePos(pFile, (fpos_T*) &iPosition);
                mexPrintf("MDH_SYNCDATA: Position is file %" FMT64 "d \n", iPosition);
                //mexPrintf("MDH_SYNCDATA: Position is file %u (out)\n", ftell(pFile));
                lMDHPos += static_cast<long long>(uiDMA - 184 + sizeof(sMDH));
                continue;
            }               
            
            uiDMA = (unsigned int) (8*mdh.ushSamplesInScan+iDMASkip)*mdh.ushUsedChannels;
            
            /*
            //iStatus = fseek(pFile, uiDMA, SEEK_CUR);
            //iStatus = fread(afDummy, 1, uiDMA, pFile);
            if (iStatus != 0) {
                if(lFilesize < static_cast<long long>(lMDHPos + uiDMA))
                {
                    mexPrintf("EOF reached after %10.0f Iterations at %u (fseek)\n", *pdCount, lMDHPos);
                    mexPrintf("uiDMA = %u\n", uiDMA);
                    delete [] measOffset;
                    delete [] measLength; 
                    fclose(pFile);
                    return;
                }
            }*/
            
            getFilePos(pFile, (fpos_T*) &iPosition);
            iOffset = iPosition + uiDMA;
            if(iOffset > iFileSize)
            {
                mexPrintf("EOF reached after %10.0f Iterations at %u (fseek)\n", *pdCount, lMDHPos);
                delete [] measOffset;
                delete [] measLength; 
                fclose(pFile);
                return;
            }
            setFilePos(pFile, (fpos_T*) &iOffset);       
            
            lMDHPos += static_cast<long long>(uiDMA + sizeof(sMDH));
            (*pdCount)++;
            //mexPrintf("%.2f %%", double(iPosition-measOffset[i])/double(measLength[i]));
        }
    }
    fclose(pFile);
    delete [] measOffset;
    delete [] measLength;
    mexPrintf("Terminating normally\n");
}
// ========================================================================
// *** END OF MAIN MEX FUNCTION fMeasCreateLUT
// ========================================================================