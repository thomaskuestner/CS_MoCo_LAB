
/***************************************************************************/
/* Name:          compute_weights_mex.c                                    */
/* Description:                                                            */
/* Date:          08-08-24                                                 */
/* Author:        Xavier Bresson (xbresson@math.ucla.edu)                  */
/***************************************************************************/


//  mex -v -g SBNLTV_mex.c 


#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <math.h>
#include <time.h>


 
#define YES 1
#define NO 0

#define X(ix,iy) (ix)*iNy+ (iy)
#define X2(ix,iy) (ix)*im+ (iy)
#define X3(ix,iy) (ix)*iw+ (iy)
#define X4(ix,iy,i) (i)*iNyNx+ (iy)*iNx+ (ix)
#define X5(iX,i) (i)*iNyNx+ iX
#define Xb(ix,iy) (iy)*iNx+ (ix)



#define ABS(x) ( x >= 0.0 ? x : -x )
#define SQR(x) (x)*(x)



float SQRT(float number) {
    long i;
    float x, y;
    const float f = 1.5F;
    
    x = number * 0.5F;
    y  = number;
    i  = * ( long * ) &y;
    i  = 0x5f3759df - ( i >> 1 );
    y  = * ( float * ) &i;
    y  = y * ( f - ( x * y * y ) );
    y  = y * ( f - ( x * y * y ) );
    return number * y;
}











/****************************************/
extern void mexFunction(int iNbOut, mxArray *pmxOut[],
int iNbIn, const mxArray *pmxIn[])
{
    
  /* iNbOut: number of outputs
     pmxOut: array of pointers to output arguments */
    
  /* iNbIn: number of inputs
     pmxIn: array of pointers to input arguments */
    
    
    float   *pfIm0, *pfu, *pfuNew, *pfv, *pfvNew, *pfVecGeneralParameters, *pfTemp, *pfW;
    float   *pfDivWP, *pfpy, *pfW2, f, fSum, fLambda, fMu, fuc, fw, fIm0c;
    float   *pfp, *pfpNew, fdij, fdji, *pfF, *pfGradWF, fF, *pfNormF, *pfInpaintRegion;
    float   fMaxDivWP, fMaxP, *pfd, *pfdNew, fctST, fSum1, fSum2, fDivWD, fDen, fNormU, fTemp;
    float   *pfuTemp, fDiff, fTol, *pfbd, *pfbdNew, fSqrtw, fGuij;
    float   *pfb, *pfbNew, fbij, fbji, *pfGu;
    int     iNy, iNx, iNdim, iDim[3], iDisplay;
    int     iNbNeigh, iN3, im, im2, iw, iw2, ic1;
    int     iy, ix, iy2, ix2, iy3, ix3, iy4, ix4, ixMax, iyMax, i, j;
    int     iNbIters, iIter, ixY, iyY, iCptX;
    int     *piY, iNyNx, iX, iXb, iXi, iIterU;
    time_t  start_time, end_time;
    
    
    start_time = clock();
    
    
    
  /* Get the  */
    pfu = mxGetData(pmxIn[0]);
    
  /* Get the  */
    pfd = mxGetData(pmxIn[1]);
    
  /* Get the  */
    pfb = mxGetData(pmxIn[2]);
    
  /* Get the  */
    pfIm0 = mxGetData(pmxIn[3]);
    
  /* Get the  */
    pfW = mxGetData(pmxIn[4]);
    
  /* Get the  */
    piY = mxGetData(pmxIn[5]);
    
  /* Get the vector of the general paremeters */
    pfVecGeneralParameters = mxGetData(pmxIn[6]);
    
    
  /* Get the displaying's indicator of different messages, values, etc */
    iDisplay = (int) pfVecGeneralParameters[0];
    iNy = (int) pfVecGeneralParameters[1];
    iNx = (int) pfVecGeneralParameters[2];
    im = (int) pfVecGeneralParameters[3];
    iw = (int) pfVecGeneralParameters[4];
    iNbNeigh = (int) pfVecGeneralParameters[5];
    fLambda = pfVecGeneralParameters[6];
    fMu = pfVecGeneralParameters[7];
    iNbIters = (int) pfVecGeneralParameters[8];
    fTol = pfVecGeneralParameters[9];
    //mexPrintf("iNy= %i, iNx= %i, im= %i, iw= %i, iNbNeigh= %i\n",iNy,iNx,im,iw,iNbNeigh);
    //mexPrintf("fLambda= %.6f, fMu= %.6f, iNbIters= %i, fTol= %.6f\n\n",fLambda,fMu,iNbIters,fTol);
    
    


    
    iNdim = 2;
    iDim[0] = iNy;
    iDim[1] = iNx;
    
    pmxOut[0] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfuNew = mxGetData(pmxOut[0]);
    
    
    iNdim = 3;
    iDim[0] = iNx;
    iDim[1] = iNy;
    iDim[2] = iNbNeigh;
    
    pmxOut[1] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfdNew = mxGetData(pmxOut[1]);
    
    pmxOut[2] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfbNew = mxGetData(pmxOut[2]);
    
    
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = 10;
    
    pmxOut[3] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfTemp = mxGetData(pmxOut[3]);
    
    
    
    
    im2 = (im-1)/2;
    iw2 = (iw-1)/2;
    ic1 = im2+iw2;
    //mexPrintf("im2= %i, iw2= %i, ic1= %i\n",im2,iw2,ic1);
    
    
   
    
    pfuTemp = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!pfuTemp)
        mexPrintf("Memory allocation failure\n");
    
    pfGu = (float *) calloc( (unsigned)(iNbNeigh), sizeof(float) );
    if (!pfGu)
        mexPrintf("Memory allocation failure\n");
    
    
    iNyNx = iNy* iNx;
    
    for (iy=0; iy< iNy; iy++)
        for(ix=0; ix< iNx; ix++)
    {
        pfuNew[X(ix,iy)] = pfu[X(ix,iy)];
        pfuTemp[X(ix,iy)] = pfu[X(ix,iy)];
        for (i=0; i<iNbNeigh; i++)
        {
            pfdNew[X4(ix,iy,i)] = pfd[X4(ix,iy,i)];
            pfbNew[X4(ix,iy,i)] = pfb[X4(ix,iy,i)];
        }
        }


    

    fctST = 1./ fLambda;
    iIter=0;
    
    for (iIter=0; iIter< iNbIters; iIter++)
//     do
    {
        
        // u
        for (iIterU=0; iIterU< 2; iIterU++)
            for (iy=0; iy< iNy; iy++)
                for(ix=0; ix< iNx; ix++)
        {
            
            iX = X(ix,iy);
            iXb = Xb(ix,iy);
            fSum1 = 0.0;
            fSum2 = 0.0;
            fDivWD = 0.0;
            for (i=0; i<iNbNeigh; i++)
            {
                ixY = piY[X5(iXb,1+2*i)];
                iyY = piY[X5(iXb,1+2*i+1)];
                fw = pfW[X5(iXb,2*i)];
                fSqrtw = pfW[X5(iXb,2*i+1)];
                iXi = X5(iXb,i);
                fdij = pfdNew[iXi];
                fbij = pfbNew[iXi];
                fdji = pfdNew[X4(ixY,iyY,i)]; 
                fbji = pfbNew[X4(ixY,iyY,i)];
                fSum1 += fw;
                fSum2 += fw* pfuNew[X(ixY,iyY)];
                fDivWD += fSqrtw*( fdij-fdji );
                fDivWD -= fSqrtw*( fbij-fbji );
            }
            fDen = fMu + fLambda*fSum1;
            pfuNew[iX] = (fLambda*fSum2 + fMu*pfIm0[iX] - fLambda*fDivWD)/ fDen;
                }


        // d
        for (iy=0; iy< iNy; iy++)
            for(ix=0; ix< iNx; ix++)
        {
            iXb = Xb(ix,iy);
            fNormU = 0.0;
            for (i=0; i<iNbNeigh; i++)
            {
                ixY = piY[X5(iXb,1+2*i)];
                iyY = piY[X5(iXb,1+2*i+1)];
                fw = pfW[X5(iXb,2*i)];
                fSqrtw = pfW[X5(iXb,2*i+1)];
                iXi = X5(iXb,i);
                fbij = pfbNew[iXi];
                fbji = pfbNew[X4(ixY,iyY,i)];
                fGuij = fSqrtw* (pfuNew[X(ixY,iyY)]-pfuNew[X(ix,iy)]);
                fNormU += SQR(fGuij+ fbij);
                pfdNew[iXi] = fGuij+ fbij;
                pfGu[i] = fGuij;
            }
            fNormU = SQRT(fNormU);
            //mexPrintf("fNormU= %.6f\n",fNormU);
            if ( fNormU<fctST )
                for (i=0; i<iNbNeigh; i++)
            {
                iXi = X5(iXb,i);
                pfdNew[iXi] = 0.0; // d
                pfbNew[iXi] += pfGu[i]; // b
                }
            else
            {
                fTemp = fNormU-fctST; fTemp /= fNormU;
                for (i=0; i<iNbNeigh; i++)
                {
                    iXi = X5(iXb,i);
                    pfdNew[iXi] *= fTemp; // d
                    pfbNew[iXi] += pfGu[i] - pfdNew[iXi]; // b
                }
            }
            }

        
        
        fDiff = 0.0;
        for (iy=0; iy< iNy; iy++)
            for(ix=0; ix< iNx; ix++)
        {
            f = SQR(pfuNew[X(ix,iy)] - pfuTemp[X(ix,iy)]);
            fDiff += f;
            pfuTemp[X(ix,iy)] = pfuNew[X(ix,iy)];
            }
        fDiff /= (float)(iNx*iNy);
        fDiff = sqrt(fDiff);
      if(iDisplay==YES)  mexPrintf("|u-uold|= %.6f\n",fDiff);


        
    } // END do
//     while(fDiff>fTol && iIter++<100);

    
    
    
    free( (float *) pfuTemp );
    free( (float *) pfGu );
    
    
    
    
    end_time = clock();
    if (iDisplay == YES)
    {
        mexPrintf("\nComputing Time for NL-TV with Split-Bregman= %.3f sec\n \n",difftime(end_time,start_time)/1000);
    }
    

}
/****************************************/






/**************************************** End of file ****************************************/
