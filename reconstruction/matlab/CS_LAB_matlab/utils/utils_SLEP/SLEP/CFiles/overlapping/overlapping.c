#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "overlapping.h"

/*
 * -------------------------------------------------------------------
 *                       Function and parameter
 * -------------------------------------------------------------------
 *
 * In this file, we focus solving the following problem
 *
 * 1/2 \|x-v\|^2 + \lambda_1 \|x\|_1 + \lambda_2 \sum w_i \|x_{G_i}\|,
 *
 * where x and v are of dimension p,
 *       w >0, and G_i contains the indices for the i-th group
 *
 * The file is implemented in the following in Matlab:
 *
 * x=overlapping(v, p, g, lambda1, lambda2, w, G, Y, maxIter, flag);
 *
 * x and v are vectors of dimension p
 *
 * g denotes the number of groups
 *
 * lambda1 and labmda2 are non-negative regularization paramter
 *
 * G is a vector containing the indices for the groups
 * G_1, G_2, ..., G_g
 *
 * w is a 3xg matrix
 * w(1,i) contains the starting index of the i-th group in G
 * w(2,i) contains the ending index   of the i-th group in G
 * w(3,i) contains the weight for the i-th group
 *
 * Y is the dual of \|x_{G_i}\|, it is of the same size as G
 *
 * maxIter is the maximal number of iteration
 *
 * flag=0, we apply the projected gradient descent
 *
 * flag=1, we apply the accelerated gradient descent with adaptive line search
 *           (see our KDD'09 paper)
 *
 * Note: 
 * 
 *  1. One should ensure w(2,i)-w(1,i)+1=|G_i|. 
 *      !! The program does not check w(2,i)-w(1,i)+1=|G_i|.!!
 *
 *  2. The index in G and w starts from 0
 *
 * -------------------------------------------------------------------
 *                       History:
 * -------------------------------------------------------------------
 *
 * Composed by Jun Liu on May 17, 2010
 *
 * For any question or suggestion, please email j.liu@asu.edu or
 *                                              jun.liu.80@gmail.com
 *
 */

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
    double*		v		=	mxGetPr(prhs[0]);
    int			p		=   (int) mxGetScalar(prhs[1]);
    int			g		=   (int) mxGetScalar(prhs[2]);
    double		lambda1 =         mxGetScalar(prhs[3]);
    double		lambda2	=         mxGetScalar(prhs[4]);
	double*		w   	=	      mxGetPr(prhs[5]);
	double*		G   	=	      mxGetPr(prhs[6]);
	double*		Y   	=	      mxGetPr(prhs[7]);
	int			maxIter	=   (int) mxGetScalar(prhs[8]);
	int			flag	=   (int) mxGetScalar(prhs[9]);
	double		tol	=             mxGetScalar(prhs[10]);
    
	double *x;
    double *gap;
    double *penalty2;
    
	/* set up output arguments */
	plhs[0] = mxCreateDoubleMatrix(p,1,mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1,5,mxREAL);
    
	x           = mxGetPr(plhs[0]);
    gap         = mxGetPr(plhs[1]);
    penalty2    = mxGetPr(plhs[2]);
    /*
     * On May 27, 2010, we add three other elements to penalty2: pp, gg, and iterNum
     * 
     * where pp is the number of non-zero elements and 
     *       gg is the number of non-zero groups after identifying groups
     *
     * and   iterNum is the number of iteration consumed for solving the dual problem
     */
    
    overlapping(x, gap, penalty2,
            v,  p, g, lambda1, lambda2,
            w, G, Y, maxIter, flag,tol);
}

