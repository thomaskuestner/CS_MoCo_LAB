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
 * treeNorm compute
 *
 *        \sum \lambda_i \|x_{G_i}\|,
 *
 * where x is of dimension n,
 *       \lambda_i >=0, and G_i's follow the tree structure
 *
 * The file is implemented in the following in Matlab:
 *
 * tree_norm=treeNorm(x, n, ind,nodes);
 *
 * ind is a 3 x nodes matrix.
 *       Each column corresponds to a node.
 *
 *       The first element of each column is the starting index,
 *       the second element of each column is the ending index
 *       the third element of each column corrreponds to \lambbda_i.
 *
 * -------------------------------------------------------------------
 *                       Notices:
 * -------------------------------------------------------------------
 *
 * 1. The nodes in the parameter "ind" should be given in the 
 *    either
 *           the postordering of depth-first traversal
 *    or 
 *           the reverse breadth-first traversal.
 *
 * 2. When each elements of x are penalized via the same L1 
 *    (equivalent to the L2 norm) parameter, one can simplify the input
 *    by specifying 
 *           the "first" column of ind as (-1, -1, lambda)
 *
 *    In this case, we treat it as a single "super" node. Thus in the value
 *    nodes, we only count it once.
 *
 * 3. The values in "ind" are in [1,n].
 *
 * 4. The third element of each column should be positive. The program does
 *    not check the validity of the parameter. 
 *
 *    It is still valid to use the zero regularization parameter.
 *    In this case, the program does not change the values of 
 *    correponding indices.
 *    
 *
 * -------------------------------------------------------------------
 *                       History:
 * -------------------------------------------------------------------
 *
 * Composed by Jun Liu on April 21, 2010
 *
 * For any question or suggestion, please email j.liu@asu.edu.
 *
 */






void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
    double*		x		=	mxGetPr(prhs[0]);
    int			n		=   (int) mxGetScalar(prhs[1]);
	double*		ind   	=	mxGetPr(prhs[2]);
	int			nodes	=   (int) mxGetScalar(prhs[3]);
    
	double *tree_norm;    
    
	/* set up output arguments */
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    
	tree_norm = mxGetPr(plhs[0]);
	treeNorm(tree_norm, x, n, ind, nodes);
}

