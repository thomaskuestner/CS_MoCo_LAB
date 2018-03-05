#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"
#include "epph.h" /* This is the head file that contains the implementation of the used functions*/

/*
 Euclidean Projection onto l_{2,1} Ball
 
        min  1/2 ||X- V||_2^2
        s.t. ||X||_{2,1} <= z
 
which is converted to the following zero finding problem
 
        f(lambda)= \sum_i ( max( |v^i|-lambda,0) )-z=0

		v^i denotes the i-th row of V
 
 Usage:
 [x, lambda, iter_step]=ep21d(y, n, k, z, lambda0);
 
 */


void ep21d(double * x, double *root, int * steps, double * v, int n, int k, double z, double lambda0)
{
    int i, j, tn=n*k;
    double *vnorm=(double *)malloc(sizeof(double)*n);
	double *vproj=(double *)malloc(sizeof(double)*n);
	double t;

	/* compute the 2 norm of each group
	*/

	for(j=0;j<n;j++){
		t=0;
		for(i=j; i< tn; i+=n)
			t+= v[i]* v[i];
		vnorm[j]=sqrt(t);
	}



    eplb(vproj, root, steps, vnorm, n, z, lambda0);

	/* compute x
	*/

	if (*root==0){
		for(i=0;i<tn;i++)
			x[i]=v[i];
	}
	else{
		for (j=0;j<n;j++){
			if ( vnorm[j] <= *root){
				for(i=j; i< tn; i+=n)
					x[i]=0;
			}
			else{
				t=1- *root/ vnorm[j];
				for(i=j; i< tn; i+=n)
					x[i]=t* v[i];
			}
		}
	}

	free(vnorm);
	free(vproj);

}


void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* v=            mxGetPr(prhs[0]);
    int     n=            (int) mxGetScalar(prhs[1]);
	int     k=            (int) mxGetScalar(prhs[2]);
    double  z=            mxGetScalar(prhs[3]);
    double  lambda0=      mxGetScalar(prhs[4]);
    
    double *x, *lambda;
    int *iter_step;
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix(n,k,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[2] = mxCreateNumericMatrix(1,1, mxINT32_CLASS, 0);
    
    x=mxGetPr(plhs[0]);
    lambda=mxGetPr(plhs[1]);
    iter_step=(int *) mxGetPr(plhs[2]);
    

    ep21d(x, lambda, iter_step, v, n, k, z, lambda0);
}

