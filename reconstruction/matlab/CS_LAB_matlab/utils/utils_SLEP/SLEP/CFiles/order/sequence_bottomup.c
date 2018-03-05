#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"

#include "sequence.h"



/*
 * In this file, we propose an algorithm for solving the problem:
 *
 * min   1/2 \|x - u\|^2
 * s.t.  x1 \ge x2 \ge x3 \ge ... \ge xn \ge 0
*/

/* 
 * We write the wrapper for calling from Matlab
 *
 * x=sequence_bottomup(u,n);
*/

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* u=            mxGetPr(prhs[0]);
    int     n=   (int )   mxGetScalar(prhs[1]);	
    
    double *x;
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix( n, 1, mxREAL); 	
    x=  mxGetPr(plhs[0]);

	sequence_bottomup(x,u,n);
}

