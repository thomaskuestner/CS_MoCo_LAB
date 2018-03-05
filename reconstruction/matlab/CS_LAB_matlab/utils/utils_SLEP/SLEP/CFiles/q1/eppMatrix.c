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

void eppMatrix(double *X, double * V, int k, int n, double rho, double p){
    int i, j, *iter_step;
    double *v, *x;
    double c0, c;
    
    v=(double *)malloc(sizeof(double)*n);
    x=(double *)malloc(sizeof(double)*n);
    iter_step=(int *)malloc(sizeof(int)*2);
    
    /*
     *X and V are k x n matrices in matlab, stored in column priority manner
     *x corresponds a row of X
     */
    
            
    c0=0;
    for(i=0; i<k; i++){
       
        for(j=0; j<n; j++)
            v[j]=V[i + j*k];
        
        epp(x, &c, iter_step, v, n, rho, p, c0);
        c0=c;
    
        for(j=0; j<n; j++)
            X[i + j*k]=x[j];
    }
    
    free(v);
    free(x);
    free(iter_step);    
}

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* V=            mxGetPr(prhs[0]);
    int     k=     (int)  mxGetScalar(prhs[1]);
    int     n=     (int)  mxGetScalar(prhs[2]);
    double  rho=          mxGetScalar(prhs[3]);
	double  p  =          mxGetScalar(prhs[4]);
    
    double *X;
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix( k, n, mxREAL);
    
    X=mxGetPr(plhs[0]);
    

	eppMatrix(X, V, k, n, rho, p);
}

