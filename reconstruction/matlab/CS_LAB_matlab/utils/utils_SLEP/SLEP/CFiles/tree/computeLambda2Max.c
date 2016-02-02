#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "altra.h"

/*
 * compute
 *  lambda2_max=computeLambda2Max(x,n,ind,nodes);
 */


void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
    double*		x		=	mxGetPr(prhs[0]);
    int			n		=   (int) mxGetScalar(prhs[1]);
	double*		ind   	=	mxGetPr(prhs[2]);
	int			nodes	=   (int) mxGetScalar(prhs[3]);
    
	int i;
	double *lambda2_max;    
    
	/* set up output arguments */
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    
	lambda2_max = mxGetPr(plhs[0]);
    computeLambda2Max(lambda2_max,x,n,ind,nodes);
}