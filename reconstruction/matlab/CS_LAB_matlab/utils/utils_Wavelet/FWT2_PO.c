/*

 FWT2_PO.C	.MEX file corresponding to fwt2_po.m
		Periodized Wavelet Transform

  The calling syntax is:

			wc = fwt2_po(image,L,qmf)


  David Donoho
  Copyright (c) 1993  David Donoho
  All Rights Reserved

*/

#include <math.h>
#include "mex.h"
#include "wavelab.h"

void dpwt2(double *sig,int nr,int nc,int ell,int J,
   double *hpf,double *lpf,int lenfil,double *wc,double *temp);


#define DOUBLE		double
#define INT			int

/* Input Arguments */

#define	Sig_IN	prhs[0]
#define	LLL_IN	prhs[1]
#define  LPF_IN prhs[2]


/* Output Arguments */

#define	WC_OUT	plhs[0]

INT nlhs, nrhs;
mxArray *plhs[], *prhs[];

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	DOUBLE	*hpf,*lpf;
	DOUBLE	*sig,*wcp,*tmp;
	unsigned int	m,n;
	int nr,nc,nn,J,lenfil,ell;
	mxArray *temp, *hpfmat;



	/* Check for proper number of arguments */

	if (nrhs != 3) {
		mexErrMsgTxt("FWT2_PO requires 3 input arguments.");
	} else if (nlhs != 1) {
		mexErrMsgTxt("FWT2_PO requires one output argument.");
	}


	/* Check the dimensions of signal.  signal can be n X 1 or 1 X n. */

	m  = mxGetM(Sig_IN);
	n = mxGetN(Sig_IN);
  if(m != n){
      mexErrMsgTxt("FWT2_PO requires a square array");
  }
  nr = (int) m; nc = nr;
  J = 0;
  for( nn = 1; nn < nr;  nn *= 2 )  
         J ++;
  if(  nn  !=  nr){
		mexErrMsgTxt("FWT2_PO requires dyadic length sides");
	}


	/* Create a matrix for the return argument */

	WC_OUT = mxCreateDoubleMatrix(nr, nr, mxREAL);
	temp   = mxCreateDoubleMatrix(nr, 3, mxREAL);

	/* Assign pointers to the various parameters */

	wcp = mxGetPr(WC_OUT);
	tmp = mxGetPr(temp);

	sig = mxGetPr(Sig_IN);
    ell = floor ((mxGetPr(LLL_IN))[0] + .5);   /* should check whether this is in range */
    lpf = mxGetPr(LPF_IN);
    lenfil = (int) (mxGetM(LPF_IN) * mxGetN(LPF_IN));   /* should check this */
	hpfmat = mxCreateDoubleMatrix((unsigned int) lenfil,  1, mxREAL);
	hpf    = mxGetPr(hpfmat);
	mirrorfilt(lpf,hpf,lenfil);


	/* Do the actual computations in a subroutine */

	dpwt2(sig,nr,nc,ell,J,hpf,lpf,lenfil,wcp,tmp);
	mxDestroyArray(temp);
	mxDestroyArray(hpfmat);
}


void dpwt2(sig,nr,nc,ell,J,hpf,lpf,lenfil,wc,temp)
DOUBLE sig[],hpf[],lpf[],wc[],temp[];
int  nr,nc,ell,lenfil,J;
{
        DOUBLE *wcplo,*wcphi,*templo,*temphi;
        int k,j,nj;
        copydouble(sig,wc,nr*nc);
		templo = &temp[nr];
		temphi = &temp[2*nr];

               nj = nr;  
               for( j=(J-1); j >= ell; --j){
					   for( k=0; k < nj; k++){
					       wcplo = &wc[k*nr];
                       	   wcphi = &wc[k*nr + nj/2];
					       copydouble(wcplo,temp,nj);
                           downlo(temp, nj, lpf,lenfil,wcplo);
					       downhi(temp, nj, hpf,lenfil,wcphi);
					   }
					   for( k=0; k < nj; k++){
					       unpackdouble(wc,nj,nc,k,temp);
                           downlo(temp, nj, lpf,lenfil,templo);
					       downhi(temp, nj, hpf,lenfil,temphi);
						   packdouble(templo,nj/2,nc,k,wc);
						   packdouble(temphi,nj/2,nc,k,&wc[nj/2*nr]);
					   }
					   nj = nj/2;
               }
}

void unpackdouble(x,n,nc,k,y)
DOUBLE x[],*y;
int n,nc,k;
{  int i;
   for( i=0; i < n; i++)
   		*y++ = x[k+nc*i];
 }
 
void packdouble(x,n,nc,k,y)
DOUBLE *x,y[];
int n,nc,k;
{  int i;
   for( i=0; i < n; i++)
		 y[k+nc*i] = *x++;
 }


void copydouble(x,y,n)
DOUBLE *x,*y;
int n;
{
   while(n--) *y++ = *x++;
 }
 
#include "downhi.c"
#include "downlo.c"
#include "mirrorfilt.c"
         


			
          
