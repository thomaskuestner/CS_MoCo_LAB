#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"

/*
-------------------------- Function proximalRegC -----------------------------

 min_x 0.5*\|x - d\|^2 + \sum_i r_i(x)

 Usage (in matlab):
 [x]=proximalRegC(d, n, lambda, theta,type); n is the dimension of the vector d

 type = 1: Capped L1 regularizer (CapL1) 
           r_i(x) = lambda*\min(|x_i|,theta), (theta > 0, lambda >= 0)

 type = 2: Log Sum Penalty (LSP)
           r_i(x) = lambda*\sum_i log(1 + |x_i|/theta), (theta > 0, lambda >= 0)

 type = 3: Smoothly Clipped Absolute Deviation (SCAD)
           r_i(x) = lambda*|x_i|, if |x_i|<=lambda
           r_i(x) = (-x_i^2 + 2*theta*lambda*|x_i| - lambda^2)/(2(theta - 1)), if lambda<=|x_i|<=theta*lambda
           r_i(x) = 0.5*(theta + 1)*lambda^2, if |x_i| > theta*lambda, (theta > 2, lambda >= 0)

 type = 4: Minimax Concave Penalty (MCP)
           r_i(x) = lambda*|x_i| - 0.5*x_i^2/theta, if |x_i|<=theta*lambda
           r_i(x) = 0.5*theta*lambda^2, if |x_i| > theta*lambda, (theta > 0, lambda >= 0)

 default: type = 1

-------------------------- Function proximalRegC -----------------------------

-------------------------- Reference -----------------------------------------

[1] Pinghua Gong, Changshui Zhang, Zhaosong Lu, Jianhua Huang, Jieping Ye,
    A General Iterative Shrinkage and Thresholding Algorithm for Non-convex
    Regularized Optimization Problems. ICML 2013.

-----------------------------------------------------------------------------

   Copyright (C) 2012-2013 Pinghua Gong

   For any problem, please contact Pinghua Gong via pinghuag@gmail.com
*/

long mymin(double *x, long n)
{
    long i,ind;
	double temp;
	temp = x[0];
	ind = 0;
	for (i=1;i<n;i++)
	{
		if (x[i] < temp)
		{
			ind = i;
			temp = x[i];
		}

	}
	return ind;
}


void proximalCapL1(double *x, double *d, long n, double lambda, double theta)
{
    long i;
	double u,x1,x2;
	for(i=0;i<n;i++)
	{   
		u = fabs(d[i]);
		x1 = max(u,theta); 
		x2 = min(theta, max(0, u - lambda));
		if (0.5*(x1 + x2 - 2*u)*(x1 - x2) + lambda*(theta - x2) < 0)
			x[i] = x1;
		else
			x[i] = x2;
		x[i] = d[i]>= 0 ? x[i] : -x[i];

	}	
	return;
}

void proximalLSP(double *x, double *d, long n, double lambda, double theta)
{
    long i;
	double v,u,z,sqrtv;
	double xtemp[3],ytemp[3];
	xtemp[0] = 0.0;

	for(i=0;i<n;i++)
	{ 
		u = fabs(d[i]);
		z = u - theta;
		v = z*z - 4.0*(lambda - u*theta);

		if (v < 0)
			x[i] = 0;
		else
		{
			sqrtv = sqrt(v);
			xtemp[1] = max(0,0.5*(z + sqrtv));
			xtemp[2] = max(0,0.5*(z - sqrtv));

			ytemp[0] = 0.5*u*u;
			ytemp[1] = 0.5*(xtemp[1] - u)*(xtemp[1] - u) + lambda*log(1.0 + xtemp[1]/theta);
			ytemp[2] = 0.5*(xtemp[2] - u)*(xtemp[2] - u) + lambda*log(1.0 + xtemp[2]/theta);

			x[i] = xtemp[mymin(ytemp,3)];
			x[i] = d[i]>= 0 ? x[i] : -x[i];
		}
	}	
	return;
}


void proximalSCAD(double *x, double *d, long n, double lambda, double theta)
{
    long i;
	double u,z,w;
	double xtemp[3],ytemp[3];
	z = theta*lambda;
	w = lambda*lambda;
	for(i=0;i<n;i++)
	{ 
		u = fabs(d[i]);
		xtemp[0] = min(lambda,max(0,u - lambda));
	    xtemp[1] = min(z,max(lambda,(u*(theta-1.0)-z)/(theta-2.0)));
	    xtemp[2] = max(z,u);

		ytemp[0] = 0.5*(xtemp[0] - u)*(xtemp[0] - u) + lambda*xtemp[0];
		ytemp[1] = 0.5*(xtemp[1] - u)*(xtemp[1] - u) + 0.5*(xtemp[1]*(-xtemp[1] + 2*z) - w)/(theta-1.0);
		ytemp[2] = 0.5*(xtemp[2] - u)*(xtemp[2] - u) + 0.5*(theta+1.0)*w;

		x[i] = xtemp[mymin(ytemp,3)];
		x[i] = d[i]>= 0 ? x[i] : -x[i];

	}	
	return;
}

void proximalMCP(double *x, double *d, long n, double lambda, double theta)
{
    long i;
	double x1,x2,v,u,z;
	z = theta*lambda;
	if (theta  >  1)
	{
		for(i=0;i<n;i++)
		{ 
			u = fabs(d[i]);
			x1 = min(z,max(0,theta*(u - lambda)/(theta - 1.0)));
			x2 = max(z,u);
			if (0.5*(x1 + x2 - 2*u)*(x1 - x2) + x1*(lambda - 0.5*x1/theta) - 0.5*z*lambda < 0)
				x[i] = x1;
			else
				x[i] = x2;
			x[i] = d[i]>= 0 ? x[i] : -x[i];
		}
    }
	else if (theta < 1)
	{
		for(i=0;i<n;i++)
		{ 
			u = fabs(d[i]);
			v = theta*(u - lambda)/(theta - 1);
			x1 = fabs(v) > fabs(v - z) ? 0.0 : z;
			x2 = max(z,u);
			if (0.5*(x1 + x2 - 2*u)*(x1 - x2) + x1*(lambda - 0.5*x1/theta) - 0.5*z*lambda < 0)
				x[i] = x1;
			else
				x[i] = x2;
			x[i] = d[i]>= 0 ? x[i] : -x[i];
		}
	}
	else
	{
		for (i=0;i<n;i++)
		{
			u = fabs(d[i]);
			x1 = lambda > u ? 0.0 : z;
			x2 = max(z,u);
			if (0.5*(x1 + x2 - 2*u)*(x1 - x2) + x1*(lambda - 0.5*x1/theta) - 0.5*z*lambda < 0)
				x[i] = x1;
			else
				x[i] = x2;
			x[i] = d[i]>= 0 ? x[i] : -x[i];	
		}
	}
	return;
}




void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* d       =            mxGetPr(prhs[0]);
    long     n      =      (long)mxGetScalar(prhs[1]);
    double  lambda  =            mxGetScalar(prhs[2]);
	double  theta   =            mxGetScalar(prhs[3]);
	int     type    =       (int)mxGetScalar(prhs[4]);
    
    double *x;


    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    
    x=mxGetPr(plhs[0]);

	switch (type)
	{
	case 1:
		proximalCapL1(x, d, n, lambda, theta);
		break;
	case 2:
		proximalLSP(x, d, n, lambda, theta);
		break;
	case 3:
		proximalSCAD(x, d, n, lambda, theta);
		break;
	case 4:
		proximalMCP(x, d, n, lambda, theta);
		break;
	default:
		proximalCapL1(x, d, n, lambda, theta);
	}
}


