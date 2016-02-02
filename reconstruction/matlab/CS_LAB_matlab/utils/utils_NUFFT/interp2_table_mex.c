/*
* interp2_table_mex.c
* Mex file for 2D periodic interpolation using tabulated interpolator.
* Copyright 03-30-2004 Yingying Zhang and Jeff Fessler, The University of Michigan
*/
#include "mex.h"
#include "math.h"

#define mxIsInt32n(mx, n) \
	( (n == mxGetM(mx) * mxGetN(mx)) && mxIsInt32(mx) )

static void interp2_table(
const double *r_ck,	/* [K1,K2] in */
const double *i_ck,
const int K1,
const int K2,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *i_h1,
const double *r_h2,	/* [J2*L2+1,1] in */
const double *i_h2,
const int J1,
const int J2,
const int L1,
const int L2,
const double *p_tm,	/* [M,2] in */
const int M,
double *r_fm,		/* [M,1] out */
double *i_fm)
{
	int mm;
	const int J_shift1 = (J1 % 2) ? (J1+1)/2 : J1/2; /* nufft_offset */
	const int J_shift2 = (J2 % 2) ? (J2+1)/2 : J2/2; /* from 0 in C */

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2);
	r_h1 += ncenter1;
	i_h1 += ncenter1;
	}
	{
	const int ncenter2 = floor(J2 * L2/2);
	r_h2 += ncenter2;
	i_h2 += ncenter2;
	}

	/* interp */
    for (mm=0; mm<M; mm++) {
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;

	/* put t_m in range [0,K-1] */
	const double tm1 = t1 - K1 * floor(t1 / K1);
	const double tm2 = t2 - K2 * floor(t2 / K2);

	int k1 = 1 + ((J1%2==1)
		? ( round(tm1) - J_shift1 )
		: ( floor(tm1) - J_shift1 ));
	const int koff2 = 1 + ((J2%2==1)
		? ( round(tm2) - J_shift2 )
		: ( floor(tm2) - J_shift2 ));

	register double sum1r = 0.;
	register double sum1i = 0.;
	int j1, j2;

	for (j1=0; j1<J1; j1++, k1++) {
		const int n1 = /* ncenter1 + */ round((tm1 - k1) * L1);
		register double coef1r = r_h1[n1];
		register double coef1i = i_h1[n1];
		const int k1mod = (k1 + K1) % K1;

		register double sum2r=0.;
		register double sum2i=0.;
		int k2 = koff2;

		for (j2=0; j2<J2; j2++, k2++) {
			const int n2 = /* ncenter2 + */ round((tm2 - k2) * L2);
			register double coef2r = r_h2[n2];
			register double coef2i = i_h2[n2];
			const int k2mod = (k2 + K2) % K2;
			const int kk = k2mod*K1+k1mod;

			/* sum2 += coef2 * ck */
			sum2r += coef2r * r_ck[kk] - coef2i * i_ck[kk];
			sum2i += coef2r * i_ck[kk] + coef2i * r_ck[kk];
		}
		/* sum1 += coef1 * sum2 */
		sum1r += coef1r * sum2r - coef1i * sum2i;
		sum1i += coef1r * sum2i + coef1i * sum2r;
	}

	*r_fm++ = sum1r;
	*i_fm++ = sum1i;
    }
}


int interp2_table_mex(
mxArray *plhs[],
const mxArray *mx_ck,
const mxArray *mx_h1,
const mxArray *mx_h2,
const mxArray *mx_J,
const mxArray *mx_L,
const mxArray *mx_tm)
{
	int nn;
	const int K1 = mxGetM(mx_ck);	/* [K1,K2] DFT coefficients */
	const int K2 = mxGetN(mx_ck);
	const int ndim = mxGetNumberOfDimensions(mx_ck);
	const int N = (ndim > 2) ?  (mxGetDimensions(mx_ck))[2] : 1;
	const int M = mxGetM(mx_tm);	/* # of time samples */

	const int *Jd = (int *) mxGetData(mx_J);
	const int *Ld = (int *) mxGetData(mx_L);

	const double *p_tm = mxGetPr(mx_tm);
	const double *r_h1 = mxGetPr(mx_h1);
	const double *i_h1 = mxGetPi(mx_h1);
	const double *r_h2 = mxGetPr(mx_h2);
	const double *i_h2 = mxGetPi(mx_h2);
	const double *r_ck = mxGetPr(mx_ck);
	const double *i_ck = mxGetPi(mx_ck);

	double *r_fm, *i_fm;

	if (N != 1)
		fprintf(stderr, "Caution: multiple realizations?");

	if (!mxIsComplex(mx_ck))
		mexErrMsgTxt("ck must be complex");

	/* J,L must be [1,2] */
	if (!mxIsInt32n(mx_J, 2))
		mexErrMsgTxt("J must be [1,2]");
	if (!mxIsInt32n(mx_L, 2))
		mexErrMsgTxt("L must be [1,2]");

	/* check size & type of tables */
	if (mxGetM(mx_h1) != Jd[0]*Ld[0]+1 || mxGetN(mx_h1) != 1) {
		fprintf(stderr, "J=%d L=%d tablelength=%d\n",
			Jd[0], Ld[0], mxGetM(mx_h1));
		mexErrMsgTxt("h1 size problem");
	}
	if (!mxIsComplex(mx_h1))
		mexErrMsgTxt("h1 must be complex");

	if (mxGetM(mx_h2) != Jd[1]*Ld[1]+1 || mxGetN(mx_h2) != 1) {
		fprintf(stderr, "J=%d L=%d tablelength=%d\n",
			Jd[1], Ld[1], mxGetM(mx_h2));
		mexErrMsgTxt("h2 size problem");
	}
	if (!mxIsComplex(mx_h2))
		mexErrMsgTxt("h2 must be complex");

	if (mxGetN(mx_tm) != 2)
		mexErrMsgTxt("tm must have two columns");

	/* create a new array and set the output pointer to it */
	if (N != 1)
		mexErrMsgTxt("N=1 done only");
	else
		plhs[0] = mxCreateDoubleMatrix(M, 1, mxCOMPLEX);
	r_fm = mxGetPr(plhs[0]);
	i_fm = mxGetPi(plhs[0]);

	for (nn=0; nn < N; ++nn) {
		interp2_table(r_ck, i_ck, K1, K2, r_h1, i_h1, r_h2, i_h2,
			Jd[0], Jd[1], Ld[0], Ld[1], p_tm, M, r_fm, i_fm);
		r_ck += K1*K2; i_ck += K1*K2;
		r_fm += M; i_fm += M;
	}

	return 1;
}


/* The gateway routine. */
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
	/* check for the proper number of arguments */
	if (nrhs != 6)
		mexErrMsgTxt("6 inputs required: (ck, h1, h2, J, L, tm)");
	if (nlhs > 1)
		mexErrMsgTxt("Less than one output arguments.");

	if (!interp2_table_mex(plhs,
		prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], prhs[5]))
		mexErrMsgTxt("interp2_table_mex failed");

	return;
}
