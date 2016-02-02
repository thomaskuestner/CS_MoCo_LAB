/* Resampling along the column (Type 1 and 2)
 *
 * Created by: Minh N. Do, March 2000
 *
 * Modified by: Yue M. Lu, Decembe 2008
 */

#include "mex.h"

/*
  function y = resampc(x, type, shift, extmod)
  % RESAMPC	Resampling along the column
  %
  %	y = resampc(x, type, shift, extmod)
  %
  % Input:
  %	x:	image that is extendable along the column direction
  %	type:	either 1 or 2 (1 for shuffering down and 2 for up)
  %	shift:	amount of shifts (typically 1)
  %     extmod: extension mode:
  %		'per' 	periodic
  %		'ref1'	reflect about the edge pixels
  %		'ref2'	reflect, doubling the edge pixels 
  %
  % Output:
  %	y:	resampled image with:
  %		R1 = [1, shift; 0, 1] or R2 = [1, -shift; 0, 1]
*/
void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *x, *px;		/* input matrix and pointer */
    double *y, *py;		/* result matrix and pointer */
    int type;			/* type of resampling */
    int s;			/* amount of shifts */
    char extmod[10];		/* extension mode */

    int i, j, k, m, n;

    /* Parse input */
    if (nrhs < 4)
	mexErrMsgTxt("Not enough input for RESAMPC!");

    x = mxGetPr(prhs[0]);
    m = (int) mxGetM(prhs[0]);
    n = (int) mxGetN(prhs[0]);

    type = (int) mxGetPr(prhs[1])[0];
    if ((type != 1) && (type != 2))
	mexErrMsgTxt("The second input (type) must be either 1 or 2");

    s = (int) mxGetPr(prhs[2])[0];

    if (!mxIsChar(prhs[3]))
	mexErrMsgTxt("EXTMOD arg must be a string");

    mxGetString(prhs[3], extmod, 10);
    
    /* Create output */
    plhs[0] = mxCreateDoubleMatrix((mwSize) m, (mwSize) n, mxREAL);
    y = mxGetPr(plhs[0]);

    px = x;
    py = y;

    if (strcmp(extmod, "per") == 0)
    {
	/* Resampling column-wise:
	 * 		y[i, j] = x[<i+sj>, j] 	if type == 1
	 * 		y[i, j] = x[<i-sj>, j] 	if type == 2
	 */
	
	for (j = 0; j < n; j++)
	{
	    /* Circular shift in each column */
	    if (type == 1)
		k = (s * j) % m;
	    else
		k = (-s * j) % m;
	    
	    /* Convert to non-negative mod if needed */
	    if (k < 0)
		k += m;
	    
	    for (i = 0; i < m; i++)
	    {
		if (k >= m)
		    k -= m;
		
		py[i] = px[k];
		
		k++;
	    }
	    
	    px += m;
	    py += m;
	}
    }

    else
	mexErrMsgTxt("Invalid EXTMOD");
}
