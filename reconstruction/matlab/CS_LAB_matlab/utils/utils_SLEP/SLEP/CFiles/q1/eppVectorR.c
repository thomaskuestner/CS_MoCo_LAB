#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"

/*
   min  1/2 ( ||x- u||_2^2 + ||t-v||_2^2 )
   s.t.  ||x_j||_2 <= t_j
 
 */

void eppVectorR(double *x, double * t, double * u, double * v, double * ind, int n, int k){
    int i, j;
    double temp;

	/* compute the 2 norm of each group
	*/

	for(j=0;j<k;j++){
		temp=0;
		for(i=(int) (ind[j]); i< (int) (ind[j+1]); i++)
			temp+= u[i]* u[i];
        temp=sqrt(temp);
        /*temp contains the 2-norm of of each row of u*/

        if(temp > fabs(v[j])){
           t[j]=(temp + v[j])/2;
           
           for(i=(int) (ind[j]); i< (int) (ind[j+1]); i++)
               x[i]= t[j] / temp * u[i];
        }
        else
           if(temp <= v[j]){
               t[j]=v[j];
                
               for(i=(int) (ind[j]); i< (int) (ind[j+1]); i++)
                   x[i]= u[i];
            }
            else{
                t[j]=0;
                
               for(i=(int) (ind[j]); i< (int) (ind[j+1]); i++)
                   x[i]=0;
            }
              
	}    
}

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* u=            mxGetPr(prhs[0]);
    double* v=            mxGetPr(prhs[1]);
	double* ind=          mxGetPr(prhs[2]);
    int     n=   (int )   mxGetScalar(prhs[3]);
    int     k=   (int )   mxGetScalar(prhs[4]);
    
    double *x, *t;
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix( n, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix( k, 1, mxREAL);
    
    x=mxGetPr(plhs[0]);
    t=mxGetPr(plhs[1]);    

	eppVectorR(x, t, u, v, ind, n, k);
}

