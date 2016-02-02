/*
* interp1_table_mex.c
* Mex file for 1D periodic interpolation using table lookup.
*
* forward direction: (for m = 1,...,M)
* f(t_m) = \sum_{k=0}^{K-1} c_k h( (t_m - k) mod K )
*
* The interpolator h is nonzero (and tabulated) for -J/2 <= t <= J/2.
*
* Copyright 2004-3 Yingying Zhang and Jeff Fessler, The University of Michigan
*/
#include "mex.h"
#include "math.h"

#define mxIsScalarInt32(mx) \
	( (1 == mxGetM(mx)) && (1 == mxGetN(mx)) && mxIsInt32(mx) )

static void interp1_table_per(
const double *r_ck,	/* [K,1] in */
const double *i_ck,
const int K1,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *i_h1,
const int J1,
const int L1,
const double *p_tm,	/* [M,1] in */
const int M,
double *r_fm,		/* [M,1] out */
double *i_fm)
{
	int mm;
	const int J_shift1 = (J1 % 2) ? (J1+1)/2 : J1/2; /* nufft_offset */

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = floor(J1 * L1/2);
	r_h1 += ncenter1;
	i_h1 += ncenter1;
	}

	/* interp */
    for (mm=0; mm<M; mm++) {
	const double t1 = *p_tm++;
	/* put t_m in range [0,K-1] */
	const double tm1 = t1 - K1 * floor(t1 / K1);
	int k1 = 1 + ((J1%2==1)
		? ( round(tm1) - J_shift1 )
		: ( floor(tm1) - J_shift1 ));
	register double sum1r = 0.;
	register double sum1i = 0.;
	int j1;

	for (j1=0; j1<J1; j1++, k1++) {
		const int n1 = /* ncenter1 + */ round((tm1 - k1) * L1);
		register double coef1r = r_h1[n1];
		register double coef1i = i_h1[n1];
		const int k1mod = (k1 + K1) % K1;

		/* complex: sum += coef * ck */
		sum1r += coef1r * r_ck[k1mod] - coef1i * i_ck[k1mod];
		sum1i += coef1r * i_ck[k1mod] + coef1i * r_ck[k1mod];
	}
	*r_fm++ = sum1r;
	*i_fm++ = sum1i;
    }
}


/*
* interp1_table_per_mex()
*/
static int interp1_table_per_mex(
mxArray *plhs[],
const mxArray *mx_ck,
const mxArray *mx_h1,
const mxArray *mx_J,
const mxArray *mx_L,
const mxArray *mx_tm)
{
	const int K = mxGetM(mx_ck);    /* # of DFT coefficients */
	const int N = mxGetN(mx_ck);    /* # of realizations */
	const int M = mxGetM(mx_tm);    /* # of time samples */

	const int J1 = *((int *) mxGetData(mx_J));
	const int L1 = *((int *) mxGetData(mx_L));

	const double *r_h1 = mxGetPr(mx_h1);
	const double *i_h1 = mxGetPi(mx_h1);
	const double *p_tm = mxGetPr(mx_tm);
	const double *r_ck = mxGetPr(mx_ck);
	const double *i_ck = mxGetPi(mx_ck);
	double *r_fm, *i_fm;
	int nn;

	if (N != 1)
		fprintf(stderr, "Caution: multiple columns?");

	/* J, L, must be scalar */
	if (!mxIsScalarInt32(mx_J))
		mexErrMsgTxt("J must be scalar int32");
	if (!mxIsScalarInt32(mx_L))
		mexErrMsgTxt("L must be scalar int32");

	/* check h table size */
	if (mxGetM(mx_h1) != J1*L1+1 || (mxGetN(mx_h1) != 1)) {
		fprintf(stderr, "J=%d L=%d tablelength=%d\n",
			J1, L1, mxGetM(mx_h1));
		mexErrMsgTxt("h size problem");
	}
	if (!mxIsComplex(mx_h1))
		mexErrMsgTxt("h must be complex.");

	if (!mxIsComplex(mx_ck))
		mexErrMsgTxt("ck must be complex.");

	if (mxGetN(mx_tm) != 1)
		mexErrMsgTxt("t_m must be col vector.");
	if (mxIsComplex(mx_tm))
		mexErrMsgTxt("tm must be real");

	/* create a new array and set the output pointer to it */
	plhs[0] = mxCreateDoubleMatrix(M, N, mxCOMPLEX);
	r_fm = mxGetPr(plhs[0]);
	i_fm = mxGetPi(plhs[0]);

	/* call the C subroutine N times; once for each realization */
	for (nn=0; nn < N; ++nn) {
		interp1_table_per(r_ck, i_ck, K, r_h1, i_h1,
			J1, L1, p_tm, M, r_fm, i_fm);
		r_ck += K; i_ck += K;
		r_fm += M; i_fm += M;
	}

	return 1;
}


/*
* The gateway routine.
* Usage: fm = function(ck, h_table, J, L, tm)
*/
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
	/* check for the proper number of arguments */
	if (nrhs != 5 )
		mexErrMsgTxt("5 inputs needed: (ck, h, J, L, tm)");
	if (nlhs > 1)
		mexErrMsgTxt("Less than one output arguments.");

	if (!interp1_table_per_mex(plhs,
		prhs[0], prhs[1], prhs[2], prhs[3], prhs[4]))
			mexErrMsgTxt("interp1_table_mex() failed");

	return;
}
