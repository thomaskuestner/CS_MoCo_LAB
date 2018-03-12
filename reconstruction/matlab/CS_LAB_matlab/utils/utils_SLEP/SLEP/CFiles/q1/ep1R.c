#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"


/*
 Euclidean Projection onto l_{2,1} Ball
 
        min  1/2 ||x- u||_2^2 + 1/2 ||t- v||_2^2
        s.t. |x|<=t

 
 Usage:
 [x, t]=ep1R(u, v, n);
 
 */


void ep1R(double * x, double *t, double * u, double * v, int n)
{
    int j;


	for(j=0;j<n;j++){

        if(fabs(u[j]) > fabs(v[j])){
			t[j]=(fabs(u[j]) + v[j])/2;

			if (u[j] >0)
				x[j]=t[j];
			else
				x[j]=-t[j];
        }
        else
           if(fabs(u[j]) <= v[j]){
               t[j]=v[j];
			   x[j]=u[j];
            }
            else{
                t[j]=x[j]=0;
            }
              
	}
}


void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* u=            mxGetPr(prhs[0]);
    double* v=            mxGetPr(prhs[1]);
    int     n=            (int) mxGetScalar(prhs[2]);
    
    double *x, *t;
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    
    x=mxGetPr(plhs[0]);
    t=mxGetPr(plhs[1]);
    

    ep1R(x, t, u, v, n);
}

