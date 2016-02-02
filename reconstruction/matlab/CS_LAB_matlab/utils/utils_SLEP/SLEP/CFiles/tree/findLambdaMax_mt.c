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
 * findLambdaMax_mt compute
 * 
 * the lambda_{max} that achieves a zero solution for
 *
 *     min  1/2 \|X-V\|^2 +  \sum_{row of X} \lambda_{\max} * \sum  w_i \|x_{G_i}\|,
 *
 * where x (a row of X) is of dimension n,
 *       w_i >=0, and G_i's follow the tree structure
 *
 * The file is implemented in the following in Matlab:
 *
 * lambdaMax=findLambdaMax_mt(v, n, ind,nodes);
 */

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
    double*		V		=	mxGetPr(prhs[0]);
    int			n		=   (int) mxGetScalar(prhs[1]);
    int         k       =   (int) mxGetScalar(prhs[2]);
    double*		ind   	=	mxGetPr(prhs[3]);
    int			nodes	=   (int) mxGetScalar(prhs[4]);
    
	double *lambdaMax;    
    
	/* set up output arguments */
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    
	lambdaMax = mxGetPr(plhs[0]);
    findLambdaMax_mt(lambdaMax,V,n,k, ind,nodes);
}
