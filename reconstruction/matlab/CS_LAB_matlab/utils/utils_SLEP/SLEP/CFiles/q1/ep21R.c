#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"


/*
 Euclidean Projection onto l_{2,1} Ball
 
        min  1/2 ||x- u||_2^2 + 1/2 ||t- v||_2^2
        s.t. ||x^j||_{2,1} <= t^j

 
 Usage:
 [x, t]=ep21R(u, v, n, k);
 
 */


void ep21R(double * x, double *t, double * u, double * v, int n, int k)
{
    int i, j, tn=n*k;
    double temp;

	/* compute the 2 norm of each group
	*/

	for(j=0;j<n;j++){
		temp=0;
		for(i=j; i< tn; i+=n)
			temp+= u[i]* u[i];
                temp=sqrt(temp);
                /*temp contains the 2-norm of of each row of u*/

        if(temp > fabs(v[j])){
           t[j]=(temp + v[j])/2;
           for (i=j; i<tn; i+=n)
               x[i]= t[j] / temp * u[i];
        }
        else
           if(temp <= v[j]){
               t[j]=v[j];
                
               for (i=j; i<tn; i+=n)
                   x[i]= u[i];
            }
            else{
                t[j]=0;
                
               for (i=j; i<tn; i+=n)
                   x[i]=0;
            }
              
	}
}


void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* u=            mxGetPr(prhs[0]);
    double* v=            mxGetPr(prhs[1]);
    int     n=            (int) mxGetScalar(prhs[2]);
	int     k=            (int) mxGetScalar(prhs[3]);
    
    double *x, *t;
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix(n,k,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    
    x=mxGetPr(plhs[0]);
    t=mxGetPr(plhs[1]);
    

    ep21R(x, t, u, v, n, k);
}

