#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "altra.h"


/*
 * -------------------------------------------------------------------
 *                       Function and parameter
 * -------------------------------------------------------------------
 *
 * findLambdaMax compute
 * 
 * the lambda_{max} that achieves a zero solution for
 *
 *     min  1/2 \|x-v\|^2 +  \lambda_{\max} * \sum  w_i \|x_{G_i}\|,
 *
 * where x is of dimension n,
 *       w_i >=0, and G_i's follow the tree structure
 *
 * The file is implemented in the following in Matlab:
 *
 * lambdaMax=findLambdaMax(v, n, ind,nodes);
 */

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
    double*		v		=	mxGetPr(prhs[0]);
    int			n		=   (int) mxGetScalar(prhs[1]);
	double*		ind   	=	mxGetPr(prhs[2]);
	int			nodes	=   (int) mxGetScalar(prhs[3]);
    
	double *lambdaMax;    
    
	/* set up output arguments */
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    
	lambdaMax = mxGetPr(plhs[0]);
    findLambdaMax(lambdaMax,v,n,ind,nodes);
}