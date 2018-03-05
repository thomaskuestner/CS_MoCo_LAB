#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"

#include "orderTree.h"



/*
 * In this file, we propose an O(n^2) algorithm for solving the problem:
 *
 * min   1/2 \|x - u\|^2
 * s.t.  x_i \ge x_j, (i,j) \in I,
 *
 * where I is the edge set of the tree
 *
 *
*/

/* 
 * We write the wrapper for calling from Matlab
 *
*/

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* u=                mxGetPr(prhs[1]);
    int     rootNum= (int )   mxGetScalar(prhs[2]);	
    int     n=       (int )   mxGetScalar(prhs[3]);	
    int yLength;
    char *FileName;    
    double *x;    
    
    yLength = mxGetN(prhs[0])+1;
    FileName = mxCalloc(yLength, sizeof(char));
    mxGetString(prhs[0], FileName, yLength);
        
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix( n, 1, mxREAL); 	
    x=  mxGetPr(plhs[0]);

	orderTree_without_nonnegative(x, FileName, u, rootNum, n);
}

