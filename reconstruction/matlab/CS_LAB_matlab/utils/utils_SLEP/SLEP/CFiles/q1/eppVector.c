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

  Here, the input and output are of Vector form.


 Written by Jun Liu, May 18th, 2009
 For any problem, please contact: j.liu@asu.edu
 
 */

void eppVector(double *x, double * v, double * ind, int k, int n, double * rho, double p){
    int i, *iter_step;
    double c0, c;
	double *px, *pv;
    
    iter_step=(int *)malloc(sizeof(int)*2);
            
    c0=0;
    for(i=0; i<k; i++){
       
        px=x+(int)ind[i];
		pv=v+(int)ind[i];
        
        epp(px, &c, iter_step, pv, (int)(ind[i+1]-ind[i]), rho[i], p, c0);
    
    }
    
    free(iter_step);    
}

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* v=            mxGetPr(prhs[0]);
	double* ind=          mxGetPr(prhs[1]);
    int     k=   (int )   mxGetScalar(prhs[2]);
    int     n=   (int )   mxGetScalar(prhs[3]);
    double* rho=          mxGetPr(prhs[4]);
	double  p  =          mxGetScalar(prhs[5]);
    
    double *x;
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix( n, 1, mxREAL);
    
    x=mxGetPr(plhs[0]);
    

	eppVector(x, v, ind, k, n, rho, p);
}

