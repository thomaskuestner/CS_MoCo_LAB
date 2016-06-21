/*************************************************************************
 * MATLAB MEX ROUTINE inplaceprod.c
 *
 * > inplaceprod(A,B) 
 * performs inplace (i.e., without allocate array) complex array product
 * > A(:) = A(:).*B(:);
 *
 * User must make sure A is not shared by other array.
 * B must have the same size as A, no check will be carried out
 *
 * A, B must be of class double or single (mixing is allowed)
 *
 * Author: Bruno Luong <brunoluong@yahoo.com>
 * History
 *  Original: 16/Sept/2009
 ************************************************************************/

#include "mex.h"
#include "matrix.h"

#define A (mxArray*)(prhs[0])
#define B prhs[1]

/* Product engine */
#define ENGINE(Ar, Ai, Br, Bi, Ari, Atype, Btype) { \
    if (Ai != NULL && Bi != NULL) { \
        for (i=0; i<n; i++) { \
            Ari = Ar[i]; \
            Ar[i] = (Atype)(Ar[i]*Br[i] - Ai[i]*Bi[i]); \
            Ai[i] = (Atype)(Ari*Bi[i] + Ai[i]*Br[i]); \
        } \
    } \
    else if (Ai != NULL) { \
        for (i=0; i<n; i++) { \
            Ar[i] *= (Atype)Br[i]; \
            Ai[i] *= (Atype)Br[i]; \
        } \
    } \
    else if (Bi != NULL) { \
        Ai = (Atype*)mxCalloc(n, sizeof(Atype)); \
        mxSetImagData(A, (void*)Ai); \
        if (Ai==NULL) \
            mexErrMsgTxt("INPLACEPROD: cannot allocate memory."); \
        for (i=0; i<n; i++) { \
            Ai[i] = (Atype)(Ar[i]*Bi[i]); \
            Ar[i] *= (Atype)Br[i]; \
        } \
    } \
    else { \
        for (i=0; i<n; i++) { \
            Ar[i] *= (Atype)Br[i]; \
        } \
    } \
} \

/* Gateway of inplaceprod */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    size_t i, n;
    double *Ard, *Aid, *Brd, *Bid;
    float *Ars, *Ais, *Brs, *Bis;
    double Ardi;
    float Arsi;
    mxClassID ClassA, ClassB;
    
    /* Get the classes of A and B */
    ClassA = mxGetClassID(A);
    ClassB = mxGetClassID(B);
    
    /* n = numel(A) */
    n = mxGetNumberOfElements(A);
    if (ClassA == mxDOUBLE_CLASS) {
        Ard = (double*)mxGetData(A);
        Aid = (double*)mxGetImagData(A);
    } else {
        Ars = (float*)mxGetData(A);
        Ais = (float*)mxGetImagData(A);
    }
    if (ClassB == mxDOUBLE_CLASS) {
        Brd = (double*)mxGetData(B);
        Bid = (double*)mxGetImagData(B);
    } else {
        Brs = (float*)mxGetData(B);
        Bis = (float*)mxGetImagData(B);
    }
    
    /* Call the macros depending of the classes of A and B */
    if ((ClassA == mxDOUBLE_CLASS) && (ClassB == mxDOUBLE_CLASS))
    {
        ENGINE(Ard, Aid, Brd, Bid, Ardi, double, double);
    } else if (ClassA == mxDOUBLE_CLASS)
    {
        ENGINE(Ard, Aid, Brs, Bis, Ardi, double, float);
    } else if (ClassB == mxDOUBLE_CLASS)
    {
        ENGINE(Ars, Ais, Brd, Bid, Arsi, float, double);
    } else
    {
        ENGINE(Ars, Ais, Brs, Bis, Arsi, float, float);
    }
    
    return;
} /* of inplaceprod */