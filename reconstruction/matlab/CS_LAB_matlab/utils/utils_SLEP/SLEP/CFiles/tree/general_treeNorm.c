#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "general_altra.h"

/*
 * Important Notice: September 20, 2010
 *
 * In this head file, we deal with the case that the features might not be well ordered.
 * 
 * If the features in the tree strucutre are well ordered, i.e., the indices of the left nodes is always less
 * than the right nodes, please refer to "altra.h".
 *
 * The advantage of "altra.h" is that, we donot need to use an explicit
 * variable for recording the indices.
 *
 *
 */

/*
 * -------------------------------------------------------------------
 *                       Functions and parameter
 * -------------------------------------------------------------------
 *
 * general_altra solves the following problem
 *
 * 1/2 \|x-v\|^2 + \sum \lambda_i \|x_{G_i}\|,
 *
 * where x and v are of dimension n,
 *       \lambda_i >=0, and G_i's follow the tree structure
 *
 * It is implemented in Matlab as follows:
 *
 * x=general_altra(v, n, G, ind, nodes);
 *
 * G contains the indices of the groups.
 *   It is a row vector. Its length equals to \sum_i \|G_i\|.
 *   If all the entries are penalized with L1 norm,
 *      its length is \sum_i \|G_i\| - n.
 *
 * ind is a 3 x nodes matrix.
 *       Each column corresponds to a node.
 *
 *       The first element of each column is the starting index,
 *       the second element of each column is the ending index
 *       the third element of each column corrreponds to \lambbda_i.
 *
 *
 *
 * The following example shows how G and ind works:
 *
 * G={ {1, 2}, {4, 5}, {3, 6}, {7, 8},
 *     {1, 2, 3, 6}, {4, 5, 7, 8}, 
 *     {1, 2, 3, 4, 5, 6, 7, 8} }.
 *
 * ind={ [1, 2, 100]', [3, 4, 100]', [5, 6, 100]', [7, 8, 100]',
 *       [9, 12, 100]', [13, 16, 100]', [17, 24, 100]' }
 *
 * Before calling the functions in this folder,
 * we have set G=G-1, and ind(1:2,:)=ind(1:2,:)-1
 *
 *
 * -------------------------------------------------------------------
 *                       Notices:
 * -------------------------------------------------------------------
 *
 * 1. The features in the tree might not be well ordered. Otherwise, you are
 *    suggested to use "altra.h".
 *
 * 2. When each elements of x are penalized via the same L1 
 *    (equivalent to the L2 norm) parameter, one can simplify the input
 *    by specifying 
 *           the "first" column of ind as (-1, -1, lambda)
 *
 *    In this case, we treat it as a single "super" node. Thus in the value
 *    nodes, we only count it once.
 *
 * 3. The values in "ind" are in [1,length(G)].
 *
 * 4. The third element of each column should be positive. The program does
 *    not check the validity of the parameter. 
 *
 * 5. The values in G should be within [0, n-1].
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
 * Composed by Jun Liu on April 20, 2010
 *
 * For any question or suggestion, please email j.liu@asu.edu.
 *
 */



void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
    double*		x		=	mxGetPr(prhs[0]);
    int			n		=   (int) mxGetScalar(prhs[1]);
	double*		G   	=	mxGetPr(prhs[2]);
	double*		ind   	=	mxGetPr(prhs[3]);
	int			nodes	=   (int) mxGetScalar(prhs[4]);
    
	double *tree_norm;    
    
	/* set up output arguments */
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    
	tree_norm = mxGetPr(plhs[0]);
	general_treeNorm(tree_norm, x, n, G, ind, nodes);
}

