/*
* interp1_table_adj_mex.c
* Mex file for *adjoint* of 1D periodic interpolation using table lookup.
*
* forward direction: (for m = 1,...,M)
* f(t_m) = \sum_{k=0}^{K-1} c_k h( (t_m - k) mod K )
*
* adjoint direction: (for k=0,...,K-1) (note complex conjugate!)
* c_k = \sum_{m=1}^M f(t_m) h^*( (t_m - k) mod K )
*
* The interpolator h is nonzero (and tabulated) for -J/2 <= t <= J/2.
*
* Copyright 2004-4-1 Jeff Fessler and Yingying Zhang, The University of Michigan
*/
#include "mex.h"
#include "math.h"

#define mxIsScalarInt32(mx) \
	( (1 == mxGetM(mx)) && (1 == mxGetN(mx)) && mxIsInt32(mx) )

static void interp1_table_adj(
double *r_ck,		/* [K,1] out */
double *i_ck,
const int K,
const double *hr,	/* [J*L+1,1] in */
const double *hi,
const int J,
const int L,
const double *pt,	/* [M,1] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm, j1;
	/* Precompute some params: Note index begins from 0 in C */
	const int ncenter = floor(J * L/2);	/* ? */
	const int J_shift = (J%2) ? (J+1)/2 : J/2;	/* nufft_offset */

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm<M; mm++) {
	const double tval = *pt++;
	const double fmr = *r_fm++;
	const double fmi = *i_fm++;

	/* put t in range [0,K-1] */
	const double tm = tval - K * floor(tval / K);
	const int koff = (J%2==1)
		? ( round(tm) - J_shift )
		: ( floor(tm) - J_shift );

	for (j1=0; j1<J; j1++) {
		const int k1 = koff + j1 + 1;
		const int int_n = ncenter + round((tm - k1) * L);
		register double coefr = *(hr+int_n);
		register double coefi = *(hi+int_n);
		const int kk = (k1 + K) % K;

		/* instead of f = h c, we have c += h^* f */
		r_ck[kk] += coefr * fmr + coefi * fmi;
		i_ck[kk] += coefr * fmi - coefi * fmr;
	}
    }
}


/*
* Usage: ck = function(fm, h_table, J, L, tm, K)
*/
static int interp1_table_adj_mex(
mxArray *plhs[],
const mxArray *mx_fm,
const mxArray *mx_ht,
const mxArray *mx_J,
const mxArray *mx_L,
const mxArray *mx_tm,
const mxArray *mx_K)
{
	const int M = mxGetM(mx_fm);	/* # of time samples */
	const int N = mxGetN(mx_fm);	/* # of realizations */
	const int J = *((int *) mxGetData(mx_J));
	const int K = *((int *) mxGetData(mx_K));
	const int L = *((int *) mxGetData(mx_L));

	const double *r_fm = mxGetPr(mx_fm);
	const double *i_fm = mxGetPi(mx_fm);
	const double *p_tm = mxGetPr(mx_tm);
	const double *r_ht = mxGetPr(mx_ht);
	const double *i_ht = mxGetPi(mx_ht);
	double *r_ck, *i_ck;
	int nn;

	if (N != 1)
		fprintf(stderr, "Caution: multiple columns?");

	/* J, L, K must be scalar */
	if (!mxIsScalarInt32(mx_J))
		mexErrMsgTxt("J must be scalar int32");
	if (!mxIsScalarInt32(mx_K))
		mexErrMsgTxt("K must be scalar int32");
	if (!mxIsScalarInt32(mx_L))
		mexErrMsgTxt("L must be scalar int32");

	/* check h table size */
	if (mxGetM(mx_ht) != J*L+1) {
		fprintf(stderr, "J=%d L=%d tablelength=%d\n", J, L, mxGetM(mx_ht));
		mexErrMsgTxt("h size problem");
	}
	if (mxGetN(mx_ht) != 1)
		mexErrMsgTxt("h must be col vector");
	if (!mxIsComplex(mx_ht))
		mexErrMsgTxt("h must be complex");

	if (!mxIsComplex(mx_fm))
		mexErrMsgTxt("fm must be complex");

	if (M != mxGetM(mx_tm) || 1 != mxGetN(mx_tm))
		mexErrMsgTxt("t_m must be Mx1 col vector");

	if (mxIsComplex(mx_tm))
		mexErrMsgTxt("tm must be real");

	/* create a new array and set the output pointer to it */
	plhs[0] = mxCreateDoubleMatrix(K, N, mxCOMPLEX);
	r_ck = mxGetPr(plhs[0]);
	i_ck = mxGetPi(plhs[0]);

	/* call the C subroutine N times; once for each realization */
	for (nn=0; nn < N; ++nn) {
		interp1_table_adj(r_ck, i_ck, K,
			r_ht, i_ht, J, L, p_tm, M, r_fm, i_fm);
		r_ck += K; i_ck += K;
		r_fm += M; i_fm += M;
	}

	return 1;
}


/*
* The gateway routine.
* Usage: ck = function(fm, h_table, J, L, tm, K)
*/
void mexFunction(
const int nlhs, mxArray *plhs[],
const int nrhs, const mxArray *prhs[])
{
	/* check for the proper number of arguments */
	if (nrhs != 6)
		mexErrMsgTxt("6 inputs needed: (f, h, J, L, t, K)");
	if (nlhs > 1)
		mexErrMsgTxt("Less than one output arguments.");

	if (!interp1_table_adj_mex(plhs,
		prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], prhs[5]))
		mexErrMsgTxt("interp1_table_adj_mex() failed");

	return;
}
