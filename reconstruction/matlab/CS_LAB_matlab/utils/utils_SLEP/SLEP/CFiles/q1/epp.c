#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"
#include "epph.h" /* This is the head file that contains the implementation of the used functions*/


/*
 Lp Norm Regularized Euclidean Projection
 
        min  1/2 ||x- v||_2^2 + rho * ||x||_p
 
 Usage (in Matlab):
 [x, c, iter_step]=epp(v, n, rho, p, c0);

 Usage in C:
 epp(x, c, iter_step, v, n, rho, p, c0);

 The function epp implements the following three functions
 epp1(x, v, n, rho) for p=1
 epp2(x, v, n, rho) for p=2
 eppInf(x, c, iter_step, v,  n, rho, c0) for p=inf
 eppO(x, c, iter_step, v,   n, rho, p) for other p

------------------------------------------------------------

  Here, the input and output are of matrix form. Each row corresponds a group


 Written by Jun Liu, May 18th, 2009
 For any problem, please contact: j.liu@asu.edu
 
 */


void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* v=            mxGetPr(prhs[0]);
    int     n=     (int)  mxGetScalar(prhs[1]);
    double  rho=          mxGetScalar(prhs[2]);
	double  p  =          mxGetScalar(prhs[3]);
	double  c0  =         mxGetScalar(prhs[4]);
    
    double *x;
    double *c;
    double *iter_step;
    int steps;
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix( n,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix( 1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix( 1,1, mxREAL);
    
    x=mxGetPr(plhs[0]);
    c=mxGetPr(plhs[1]);
    iter_step=mxGetPr(plhs[2]);    

	epp(x, c, &steps, v, n, rho, p, c0);
    *iter_step=steps;
}

