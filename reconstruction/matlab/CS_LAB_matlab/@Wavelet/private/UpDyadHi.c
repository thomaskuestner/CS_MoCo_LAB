/*

 UpDyadHi.C	Sample .MEX file corresponding to UpDyadHi.m
		HIPass Dyadic Downsampling Transform

  The calling syntax is:

			wc = UpDyadHi(signal,lpf)


  David Donoho
  Copyright (c) 1993  David Donoho
  All Rights Reserved

*/

#include <math.h>
#include "mex.h"
#include "wavelab.h"

#define DOUBLE		double
#define INT			int

/* Input Arguments */

#define	 Sig_IN	prhs[0]
#define  LPF_IN prhs[1]

/* Output Arguments */

#define	HP_OUT	plhs[0]

#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))


INT nlhs, nrhs;
mxArray *plhs[], *prhs[];

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	DOUBLE	*hpp, *lpf, *hpf, *sig;
	unsigned int	m,n,nn,mm;
	int nr,nc,J,kk,lenfil;
	mxArray *HPF_MAT;

	/* Check for proper number of arguments */

	if (nrhs != 2) {
		mexErrMsgTxt("UpDyadHi requires two input arguments.");
	} else if (nlhs != 1) {
		mexErrMsgTxt("UpDyadHi requires one output argument.");
	}


	/* Check the dimensions of signal.  signal can be n X 1 or 1 X n. */

	m = mxGetM(Sig_IN);
	n = mxGetN(Sig_IN);
	if(m == 1){
		nr = (int) n;
		nc = 1;
		nn = 2*n;
		mm = 1;
	} else {
		nr = (int) m;
		nc = (int) n;
		nn = 1;
		mm = 2*m;
	}

	J = 0;
	for( kk = 1; kk < nr;  kk *= 2 )  
		 J ++;
	if(  kk  !=  nr){
		mexErrMsgTxt("UpDyadHi requires dyadic length");
	}


	/* Create a matrix for the return argument */
	HP_OUT = mxCreateDoubleMatrix(mm,nn, mxREAL);


	/* Assign pointers to the various parameters */

	hpp = mxGetPr(HP_OUT);

	sig = mxGetPr(Sig_IN);
	lpf = mxGetPr(LPF_IN);
    lenfil = (int) (mxGetM(LPF_IN) * mxGetN(LPF_IN));   /* should check this */
	HPF_MAT = mxCreateDoubleMatrix(lenfil,1,mxREAL);
    hpf = mxGetPr(HPF_MAT);

	/* Do the actual computations in a subroutine */
    mirrorfilt(lpf,hpf,lenfil);
	uphi(sig,nr,hpf,lenfil,hpp);
}


#include "uphi.c"
#include "mirrorfilt.c"
