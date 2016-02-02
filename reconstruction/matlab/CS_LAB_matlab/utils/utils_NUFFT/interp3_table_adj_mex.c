/*
* interp3_table_adj_mex.c
* ck = interp3_table_adj_mex(fm, h1_table, h2_table J, L, tm, K)
* Mex file for *adjoint* of 2D periodic interpolation using table lookup.
*
* forward direction: (for m = 1,...,M)
* f(t_m) = \double_sum_{k1,2=0}^{K1,2-1} c_k1_k2
*		h_1( (t_m - k1) mod K1 ) h_2( (t_m - k2) mod K2 )
*
* adjoint direction: (for k1,2=0,...,K1,2-1) (note complex conjugate!)
* c_k1_k2 = \sum_{m=1}^M f(t_m) h1^*( (t_m - k1) mod K1 )
*		h2^*( (t_m - k2) mod K2 )
*
* Interpolators h1,2 is nonzero (and tabulated) for -J1,2/2 <= t <= J1,2/2.
*
* Copyright 3004-4-2 Jeff Fessler and Yingying Zhang, The University of Michigan
*/
#include "mex.h"
#include "math.h"

#define mxIsInt32n(mx, n) \
	( (n == mxGetM(mx) * mxGetN(mx)) && mxIsInt32(mx) )

static void interp3_table_adj(
double *r_ck,		/* [K1*K2 1] out */
double *i_ck,
const int K1,
const int K2,
const int K3,
const double *r_h1,	/* [J1*L1+1,1] in */
const double *i_h1,
const double *r_h2,	/* [J2*L2+1,1] in */
const double *i_h2,
const double *r_h3,	/* [J3*L3+1,1] in */
const double *i_h3,
const int J1,
const int J2,
const int J3,
const int L1,
const int L2,
const int L3,
const double *p_tm,	/* [M,3] in */
const int M,
const double *r_fm,	/* [M,1] in */
const double *i_fm)
{
	int mm;
	const int J_shift1 = (J1 % 2) ? (J1+1)/2 : J1/2; /* nufft_offset */
	const int J_shift2 = (J2 % 2) ? (J2+1)/2 : J2/2; /* from 0 in C */
	const int J_shift3 = (J3 % 2) ? (J3+1)/2 : J3/2;
	const int K12 = K1*K2;

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
	{
	const int ncenter3 = floor(J3 * L3/2);
	r_h3 += ncenter3;
	i_h3 += ncenter3;
	}

	/* initialize output to zero */
	(void) memset((void *) r_ck, 0, K1*K2*K3*sizeof(*r_ck));
	(void) memset((void *) i_ck, 0, K1*K2*K3*sizeof(*i_ck));

	/* interp */
    for (mm=0; mm<M; mm++) {
	const double t3 = p_tm[2*M];
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;

	const double fmr = *r_fm++;
	const double fmi = *i_fm++;

	/* put t_m in range [0,K-1] */
	const double tm1 = t1 - K1 * floor(t1 / K1);
	const double tm2 = t2 - K2 * floor(t2 / K2);
	const double tm3 = t3 - K3 * floor(t3 / K3);

	int k1 = 1 + ((J1%2==1)
		? ( round(tm1) - J_shift1 )
		: ( floor(tm1) - J_shift1 ));
	const int koff2 = 1 + ((J2%2==1)
		? ( round(tm2) - J_shift2 )
		: ( floor(tm2) - J_shift2 ));
	const int koff3 = 1 + ((J3%2==1)
		? ( round(tm3) - J_shift3 )
		: ( floor(tm3) - J_shift3 ));

	int j1, j2, j3;

	for (j1=0; j1<J1; j1++, k1++) {
		const int n1 = /* ncenter1 */ + round((tm1 - k1) * L1);
		register double coef1r = r_h1[n1];
		register double coef1i = i_h1[n1];
		const int k1mod = (k1 + K1) % K1;

		int k2 = koff2;

	for (j2=0; j2<J2; j2++, k2++) {
		const int n2 = /* ncenter2 + */ round((tm2 - k2) * L2);
		register double coef2r = r_h2[n2];
		register double coef2i = i_h2[n2];
		const int k2mod = (k2 + K2) % K2;
		const int k12mod = k2mod*K1 + k1mod;

		const double v2r = coef2r * fmr + coef2i * fmi;
		const double v2i = coef2r * fmi - coef2i * fmr;

		int k3 = koff3;

	for (j3=0; j3<J3; j3++, k3++) {
		const int n3 = /* ncenter3 + */ round((tm3 - k3) * L3);
		register double coef3r = r_h3[n3];
		register double coef3i = i_h3[n3];
		const int k3mod = (k3 + K3) % K3;

		/* const int kk = k3mod*K1*K2 + k2mod*K1 + k1mod; */
		const int kk = k3mod*K12 + k12mod; /* 3D array index */

		const double v3r = coef3r * v2r + coef3i * v2i;
		const double v3i = coef3r * v2i - coef3i * v2r;

		r_ck[kk] += coef1r * v3r + coef1i * v3i;
		i_ck[kk] += coef1r * v3i - coef1i * v3r;

	} /* j3 */

	} /* j2 */

	} /* j1 */
    }
}


/*
* Usage: ck = function(fm, h1, h2, h3, J, L, tm, K)
*/
static int interp3_table_adj_mex(
mxArray *plhs[],
const mxArray *mx_fm,
const mxArray *mx_h1,
const mxArray *mx_h2,
const mxArray *mx_h3,
const mxArray *mx_J,
const mxArray *mx_L,
const mxArray *mx_tm,
const mxArray *mx_K)
{
	int nn;
	const int ndim = mxGetNumberOfDimensions(mx_fm);
	const int N = (ndim > 2) ?  (mxGetDimensions(mx_fm))[2] : 1;
	const int M = mxGetM(mx_fm);	/* # of time samples */

	const int *Jd = (int *) mxGetData(mx_J);
	const int *Ld = (int *) mxGetData(mx_L);
	const int *Kd = (int *) mxGetData(mx_K);

	const double *p_tm = mxGetPr(mx_tm);
	const double *r_fm = mxGetPr(mx_fm);
	const double *i_fm = mxGetPi(mx_fm);
	const double *r_h1 = mxGetPr(mx_h1);
	const double *i_h1 = mxGetPi(mx_h1);
	const double *r_h2 = mxGetPr(mx_h2);
	const double *i_h2 = mxGetPi(mx_h2);
	const double *r_h3 = mxGetPr(mx_h3);
	const double *i_h3 = mxGetPi(mx_h3);

	double *r_ck, *i_ck;

	if (N != 1)
		fprintf(stderr, "Caution: multiple realizations?");

	if (!mxIsInt32n(mx_J, 3))
		mexErrMsgTxt("J must be int32 [1,3]");
	if (!mxIsInt32n(mx_L, 3))
		mexErrMsgTxt("L must be int32 [1,3]");
	if (!mxIsInt32n(mx_K, 3))
		mexErrMsgTxt("K must be int32 [1,3]");

	/* check h1,h2,h3 tables' sizes */
	if (mxGetM(mx_h1) != Jd[0]*Ld[0]+1 || (mxGetN(mx_h1) != 1)) {
		fprintf(stderr, "J1=%d L1=%d tablelength=%d\n",
			Jd[0], Ld[0], mxGetM(mx_h1));
		mexErrMsgTxt("h1 size problem");
	}
	if (!mxIsComplex(mx_h1))
		mexErrMsgTxt("h1 must be complex");

	if (mxGetM(mx_h2) != Jd[1]*Ld[1]+1 || (mxGetN(mx_h2) != 1)) {
		fprintf(stderr, "J2=%d L2=%d tablelength=%d\n",
			Jd[1], Ld[1], mxGetM(mx_h2));
		mexErrMsgTxt("h2 size problem");
	}
	if (!mxIsComplex(mx_h2))
		mexErrMsgTxt("h2 must be complex");

	if (mxGetM(mx_h3) != Jd[2]*Ld[2]+1 || (mxGetN(mx_h3) != 1)) {
		fprintf(stderr, "J3=%d L3=%d tablelength=%d\n",
			Jd[2], Ld[2], mxGetM(mx_h3));
		mexErrMsgTxt("h3 size problem");
	}
	if (!mxIsComplex(mx_h3))
		mexErrMsgTxt("h3 must be complex");

	if (!mxIsComplex(mx_fm))
		mexErrMsgTxt("fm must be complex");

	if (M != mxGetM(mx_tm) || 3 != mxGetN(mx_tm))
		mexErrMsgTxt("t_m must be Mx3 matrix");

	if (mxIsComplex(mx_tm))
		mexErrMsgTxt("tm must be real");

	/* create a new array and set the output pointer to it */
	if (N != 1)
		mexErrMsgTxt("N=1 done only");
	else
		plhs[0] = mxCreateDoubleMatrix(Kd[0]*Kd[1]*Kd[2], 1, mxCOMPLEX);
	r_ck = mxGetPr(plhs[0]);
	i_ck = mxGetPi(plhs[0]);

	/* call the C subroutine N times; once for each realization */
	for (nn=0; nn < N; ++nn) {
		interp3_table_adj(r_ck, i_ck, Kd[0], Kd[1], Kd[2],
			r_h1, i_h1, r_h2, i_h2, r_h3, i_h3,
			Jd[0], Jd[1], Jd[2], Ld[0], Ld[1], Ld[2],
			p_tm, M, r_fm, i_fm);
		r_ck += Kd[0]*Kd[1]*Kd[2]; i_ck += Kd[0]*Kd[1]*Kd[2];
		r_fm += M; i_fm += M;
	}

	return 1;
}


/*
* The gateway routine.
* Usage: ck = function(fm, h1_table, h2_table, h3_table, J, L, tm, K)
*/
void mexFunction(
const int nlhs, mxArray *plhs[],
const int nrhs, const mxArray *prhs[])
{

	/* check for the proper number of arguments */
	if (nrhs != 8)
		mexErrMsgTxt("8 inputs needed: (f, h1, h2, h3, J, L, t, K)");
	if (nlhs > 1)
		mexErrMsgTxt("Less than one output arguments.");

	if (!interp3_table_adj_mex(plhs,
		prhs[0], prhs[1], prhs[2], prhs[3], prhs[4],
		prhs[5], prhs[6], prhs[7]))
		mexErrMsgTxt("interp3_table_adj_mex() failed");

	return;
}
