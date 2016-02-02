/*
* interp3_table_mex.c
* Mex file for 2D periodic interpolation using tabulated interpolator.
* Copyright 03-30-2004 Yingying Zhang and Jeff Fessler, University of Michigan
*/
#include "mex.h"
#include "math.h"

#if 0
#define ifloor(x)	( (int) (x) ) /* wrong for negatives */
#define iround(x)	ifloor( (x) + 0.5 ) /* wrong for negatives! */
#define iround(x)	round(x)
#endif
#define iround(x)	( (int) ((x) >= 0 ? ((x) + 0.5) : ((x) - 0.5)) )
#define ifloor(x)	floor(x)

#define mxIsInt32n(mx, n) \
	( (n == mxGetM(mx) * mxGetN(mx)) && mxIsInt32(mx) )

static void interp3_table(
const double *r_ck,	/* [K1,K2,K3] in */
const double *i_ck,
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
double *r_fm,		/* [M,1] out */
double *i_fm)
{
	int mm;
	const int J_shift1 = (J1 % 2) ? (J1+1)/2 : J1/2; /* nufft_offset */
	const int J_shift2 = (J2 % 2) ? (J2+1)/2 : J2/2; /* from 0 in C */
	const int J_shift3 = (J3 % 2) ? (J3+1)/2 : J3/2;

	/* trick: shift table pointer to center */
	{
	const int ncenter1 = ifloor(J1 * L1/2);
	r_h1 += ncenter1;
	i_h1 += ncenter1;
	}
	{
	const int ncenter2 = ifloor(J2 * L2/2);
	r_h2 += ncenter2;
	i_h2 += ncenter2;
	}
	{
	const int ncenter3 = ifloor(J3 * L3/2);
	r_h3 += ncenter3;
	i_h3 += ncenter3;
	}

	/* interp */
    for (mm=0; mm<M; mm++) {
	const double t3 = p_tm[2*M];
	const double t2 = p_tm[M];
	const double t1 = *p_tm++;

	/* put t_m in range [0,K) */
	const double tm1 = t1 - K1 * ifloor(t1 / K1);
	const double tm2 = t2 - K2 * ifloor(t2 / K2);
	const double tm3 = t3 - K3 * ifloor(t3 / K3);

	const int koff1 = 1 + ((J1%2==1)
		? ( iround(tm1) - J_shift1 )
		: ( ifloor(tm1) - J_shift1 ));
	const int koff2 = 1 + ((J2%2==1)
		? ( iround(tm2) - J_shift2 )
		: ( ifloor(tm2) - J_shift2 ));
	int k3 = 1 + ((J3%2==1)
		? ( iround(tm3) - J_shift3 )
		: ( ifloor(tm3) - J_shift3 ));

	register double sum3r = 0.;
	register double sum3i = 0.;
	int j1, j2, j3;

	for (j3=0; j3<J3; j3++, k3++) {
		const int n3 = /* ncenter3 + */ iround((tm3 - k3) * L3);
		register double coef3r = r_h3[n3];
		register double coef3i = i_h3[n3];
		const int k3mod = (k3 + K3) % K3;

		register double sum2r = 0.;
		register double sum2i = 0.;
		int k2 = koff2;

	for (j2=0; j2<J2; j2++, k2++) {
		const int n2 = /* ncenter2 + */ iround((tm2 - k2) * L2);
		register double coef2r = r_h2[n2];
		register double coef2i = i_h2[n2];
		const int k2mod = (k2 + K2) % K2;
		const int k23mod = (k3mod*K2 + k2mod)*K1;

		register double sum1r = 0.;
		register double sum1i = 0.;
		int k1 = koff1;

	for (j1=0; j1<J1; j1++, k1++) {
		const int n1 = /* ncenter1 + */ iround((tm1 - k1) * L1);
		register double coef1r = r_h1[n1];
		register double coef1i = i_h1[n1];
		const int k1mod = (k1 + K1) % K1;

		/* const int kk = k3mod*K1*K2 + k2mod*K1 + k1mod; */
		const int kk = k23mod + k1mod;

		/* sum1 += coef1 * ck */
		sum1r += coef1r * r_ck[kk] - coef1i * i_ck[kk];
		sum1i += coef1r * i_ck[kk] + coef1i * r_ck[kk];

	} /* j1 */

		/* sum2 += coef2 * sum1 */
		sum2r += coef2r * sum1r - coef2i * sum1i;
		sum2i += coef2r * sum1i + coef2i * sum1r;

	} /* j2 */

		/* sum3 += coef3 * sum2 */
		sum3r += coef3r * sum2r - coef3i * sum2i;
		sum3i += coef3r * sum2i + coef3i * sum2r;

	} /* j3 */

	*r_fm++ = sum3r;
	*i_fm++ = sum3i;
    }
}


int interp3_table_mex(
mxArray *plhs[],
const mxArray *mx_ck,	/* [K1,K2,K3] DFT coefficients */
const mxArray *mx_h1,
const mxArray *mx_h2,
const mxArray *mx_h3,
const mxArray *mx_J,
const mxArray *mx_L,
const mxArray *mx_tm)
{
	int nn;
	const int ndim = mxGetNumberOfDimensions(mx_ck);
	const int K1 = (ndim > 0) ? (mxGetDimensions(mx_ck))[0] : 1;
	const int K2 = (ndim > 1) ? (mxGetDimensions(mx_ck))[1] : 1;
	const int K3 = (ndim > 2) ? (mxGetDimensions(mx_ck))[2] : 1;
	const int N = (ndim > 3) ?  (mxGetDimensions(mx_ck))[3] : 1;
	const int M = mxGetM(mx_tm);	/* # of time samples */

	const int *Jd = (int *) mxGetData(mx_J);
	const int *Ld = (int *) mxGetData(mx_L);

	const double *p_tm = mxGetPr(mx_tm);
	const double *r_h1 = mxGetPr(mx_h1);
	const double *i_h1 = mxGetPi(mx_h1);
	const double *r_h2 = mxGetPr(mx_h2);
	const double *i_h2 = mxGetPi(mx_h2);
	const double *r_h3 = mxGetPr(mx_h3);
	const double *i_h3 = mxGetPi(mx_h3);
	const double *r_ck = mxGetPr(mx_ck);
	const double *i_ck = mxGetPi(mx_ck);

	double *r_fm, *i_fm;

	if (N != 1)
		fprintf(stderr, "Caution: multiple realizations?");

	if (!mxIsComplex(mx_ck))
		mexErrMsgTxt("ck must be complex");

	/* J,L must be [1,3] */
	if (!mxIsInt32n(mx_J, 3))
		mexErrMsgTxt("J must be [1,3]");
	if (!mxIsInt32n(mx_L, 3))
		mexErrMsgTxt("L must be [1,3]");

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

	if (mxGetM(mx_h3) != Jd[2]*Ld[2]+1 || mxGetN(mx_h3) != 1) {
		fprintf(stderr, "J=%d L=%d tablelength=%d\n",
			Jd[2], Ld[2], mxGetM(mx_h3));
		mexErrMsgTxt("h3 size problem");
	}
	if (!mxIsComplex(mx_h3))
		mexErrMsgTxt("h3 must be complex");

	if (mxGetN(mx_tm) != 3)
		mexErrMsgTxt("tm must have 3 columns");

	/* create a new array and set the output pointer to it */
	if (N != 1)
		mexErrMsgTxt("N=1 done only");
	else
		plhs[0] = mxCreateDoubleMatrix(M, 1, mxCOMPLEX);
	r_fm = mxGetPr(plhs[0]);
	i_fm = mxGetPi(plhs[0]);

	for (nn=0; nn < N; ++nn) {
		interp3_table(r_ck, i_ck, K1, K2, K3,
			r_h1, i_h1, r_h2, i_h2, r_h3, i_h3,
			Jd[0], Jd[1], Jd[2], Ld[0], Ld[1], Ld[2],
			p_tm, M, r_fm, i_fm);
		r_ck += K1*K2*K3; i_ck += K1*K2*K3;
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
	if (nrhs != 7)
		mexErrMsgTxt("7 inputs required: (ck, h1, h2, h3, J, L, tm)");
	if (nlhs > 1)
		mexErrMsgTxt("Less than one output arguments.");

	if (!interp3_table_mex(plhs,
		prhs[0], prhs[1], prhs[2], prhs[3], prhs[4], prhs[5], prhs[6]))
		mexErrMsgTxt("interp3_table_mex failed");

	return;
}
