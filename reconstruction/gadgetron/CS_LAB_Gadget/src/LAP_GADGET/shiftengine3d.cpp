#include "shiftengine3d.h"

// define unix flag for gplot.h to work (bug in sigpack)
#if defined (__unix__) || defined (__unix)
	#define unix
#endif

#include "sigpack/sigpack.h"
Gadgetron::ShiftEngine3D::ShiftEngine3D()
{
}

Gadgetron::ShiftEngine3D::ShiftEngine3D(CubeType i_, CubeType ux_, CubeType uy_, CubeType uz_) {
	if (i_.n_rows!=ux_.n_rows ||
			i_.n_rows != uy_.n_rows ||
			i_.n_rows != uz_.n_rows) {
		throw std::runtime_error("ShiftEngine3D needs image and flow fields of same size. Rows not equal");
	}

	if (i_.n_cols!=ux_.n_cols ||
			i_.n_cols != uy_.n_cols ||
			i_.n_cols != uz_.n_cols) {
		throw std::runtime_error("ShiftEngine3D needs image and flow fields of same size. Cols not equal");
	}

	if (i_.n_slices!=ux_.n_slices ||
			i_.n_slices != uy_.n_slices ||
			i_.n_slices != uz_.n_slices) {
		throw std::runtime_error("ShiftEngine3D needs image and flow fields of same size. Slices not equal");
	}

	image = i_;
	ux = ux_;
	uy = uy_;
	uz = uz_;
	a = image.n_rows;
	b = image.n_cols;
	c = image.n_slices;

	//Generate x, y, z
	x.set_size(a,b,c);

	for (uword k = 0; k < x.n_slices; k++) {
		for (uword j = 0; j < x.n_cols; j++) {
			x.slice(k).col(j) = cs_lab_lap_arma::regspaceSimple(1, a);
		}
	}

	y.set_size(b,a,c);

	for (uword k = 0; k < y.n_slices; k++) {
		for (uword j = 0; j < y.n_cols; j++) {
			y.slice(k).col(j) = cs_lab_lap_arma::regspaceSimple(1,b);
		}
	}
	y = cs_lab_lap_arma::permuteSimple(y, 213);

	z.set_size(c,a,b);

	for (uword k = 0; k < z.n_slices; k++) {
		for (uword j = 0; j < z.n_cols; j++) {
			z.slice(k).col(j) = cs_lab_lap_arma::regspaceSimple(1,c);
		}
	}
	z = cs_lab_lap_arma::permuteSimple(z, 321);
}

CubeType Gadgetron::ShiftEngine3D::execCubicShift()
{
	return cubicInterp(x-ux, y-uy, z-uz, image);
}

CubeType Gadgetron::ShiftEngine3D::execLinShift()
{
	//Parameters for shifted-linear interpolation
	float tau = 0.2113;
	int L1 = floor(-1+tau);
	int L2 = ceil(1+tau);

	//Minimum and maximum row index needed in the interpolation formula
	CubeType x1 = x-ux;
	CubeType y1 = y-uy;
	CubeType z1 = z-uz;

	//Minimum and maximum row index needed in the interpolation formula
	//std::cout << x1.min() << std::endl;
	//std::cout << y1.min() << std::endl;
	//std::cout << z1.min() << std::endl;
	int k0 = floor(x1.min()-L2+1);
	int k1 = floor(x1.max()-L1);
	int l0 = floor(y1.min()-L2+1);
	int l1 = floor(y1.max()-L1);
	int m0 = floor(z1.min()-L2+1);
	int m1 = floor(z1.max()-L1);

	//Smallest box enclosing the image and the (x,y) positions
	int kk0 = std::min(k0,1);
	int kk1 = std::max(k1,a);
	int ll0 = std::min(l0,1);
	int ll1 = std::max(l1,b);
	int mm0 = std::min(m0,1);
	int mm1 = std::max(m1,c);

	// Convert entries to float to suppress narrowing conversion warning
	ColType extsize = {
		static_cast<float>(1-kk0),
		static_cast<float>(kk1-a),
		static_cast<float>(1-ll0),
		static_cast<float>(ll1-b),
		static_cast<float>(1-mm0),
		static_cast<float>(mm1-c)
	};

	//extsize.print("extsize");
	CubeType I0 = ext(image, extsize);
	//I0.subcube(span(0,3),span(0,3),span(0,3)).print("I0");
	I0 = I0(span(1-kk0+k0-1, 1-kk0+k1-1), span(1-ll0+l0-1, 1-ll0+l1-1), span(1-mm0+m0-1, 1-mm0+m1-1));

	int a0 = I0.n_rows;
	int b0 = I0.n_cols;
	int c0 = I0.n_slices;

	//Generate new x,y,z
	x.resize(a0,b0,c0);
	for (uword k = 0; k < x.n_slices; k++) {
		for (uword j = 0; j < x.n_cols; j++) {
			x.slice(k).col(j) = cs_lab_lap_arma::regspaceSimple(k0,k1);
		}
	}

	y.resize(b0,a0,c0);
	for (uword k = 0; k < y.n_slices; k++) {
		for (uword j = 0; j < y.n_cols; j++) {
			y.slice(k).col(j) = cs_lab_lap_arma::regspaceSimple(l0,l1);
		}
	}
	y = cs_lab_lap_arma::permuteSimple(y, 213);

	z.resize(c0,b0,a0);
	for (uword k = 0; k < z.n_slices; k++) {
		for (uword j = 0; j < z.n_cols; j++) {
			z.slice(k).col(j) = cs_lab_lap_arma::regspaceSimple(m0,m1);
		}
	}
	z = cs_lab_lap_arma::permuteSimple(z, 321);

	//Prefiltering image I0 so as to have shifted-linear interpolation
	float z0 = tau/(1-tau);

	//along 1st dimension
	//I0.subcube(span(0,3),span(0,3),span(0,3)).print("I0_before_filtering");
	I0 = 1/(1-tau)*filtering(ColType{1}, ColType({1, z0}), I0);
	//I0.subcube(span(0,3),span(0,3),span(0,3)).print("I0_after_filtering");

	//then 2nd
	I0 = cs_lab_lap_arma::permuteSimple(I0, 213);
	I0 = 1/(1-tau)*filtering(ColType{1}, ColType({1, z0}), I0);

	//finally 3rd
	I0 = cs_lab_lap_arma::permuteSimple(I0, 321);
	I0 = 1/(1-tau)*filtering(ColType{1}, ColType({1, z0}), I0);
	I0 = cs_lab_lap_arma::permuteSimple(I0, 231);

	//Andriy Myronenko's interpolation code:
	y1 = y1 - (min(vectorise(y))-1);
	x1 = x1 - (min(vectorise(x))-1);
	z1 = z1 - (min(vectorise(z))-1);

	y1 = y1-tau;
	x1 = x1-tau;
	z1 = z1-tau;

	PixelType *Z, *S, *T, *W, *F;
	int i, MN, nrows, ncols, npages, vol, ndim, newXndim, Xndim, *newdims;
	int dims[3], Xdims[3];

	Xndim = 3;
	Xdims[0] = x1.n_rows;
	Xdims[1] = x1.n_cols;
	Xdims[2] = x1.n_slices;

	ndim = 3;
	dims[0] = I0.n_rows;
	dims[1] = I0.n_cols;
	dims[2] = I0.n_slices;
	newdims = (int*) calloc(ndim-1, sizeof(int));

	MN = 1;
	for (i = 0; i < Xndim; i++) {
		MN =MN*Xdims[i];
	};

	vol = 1;
	newXndim = Xndim;
	if (ndim > 3) {
		// Check if interpolate along column vectors
		if ((Xndim == 2) && (Xdims[1] == 1)) {
			newXndim = newXndim-1;
		}

		// Allocate space for the new number of dimensions for output
		newdims = (int*) calloc(newXndim + 1, sizeof(int));

		// Copy original dimenstions
		for (i = 0; i < newXndim; i++) {
			newdims[i] = Xdims[i];
		};

		// Add the number of images as a last dimenstion
		newdims[newXndim] = dims[3];

		// Set the new number of dimenstions
		newXndim = newXndim + 1;

		vol = dims[3];
	} else {
		newdims = (int*) calloc(newXndim, sizeof(int));

		for (i = 0; i < newXndim; i++) {
			newdims[i] = Xdims[i];
		};
	}

	CubeType out(newdims[0], newdims[1], newdims[2], fill::zeros);

	nrows = dims[0];
	ncols = dims[1];
	npages = dims[2];

	Z = I0.memptr();
	S = y1.memptr();
	T = x1.memptr();
	W = z1.memptr();

	F = out.memptr();

	mirt3D_mexinterp(Z, S, T, W, F, MN, nrows, ncols, npages, vol);

	return out;
}

CubeType Gadgetron::ShiftEngine3D::ext(CubeType i_, ColType extsize_)
{
	int newa = a+extsize_(0)+extsize_(1);
	int newb = b+extsize_(2)+extsize_(3);
	int newc = c+extsize_(4)+extsize_(5);

	// Extend in x
	CubeType J;
	if (extsize_(0)>extsize_(1)) {
		J = cs_lab_lap_arma::padarraySymmetricSimple(i_, extsize_(0), 0, 0);

		J = J.tube(span(0, newa-1), span::all);
	} else {
		J = cs_lab_lap_arma::padarraySymmetricSimple(i_, extsize_(1), 0, 0);
		J = J.tube(span(J.n_rows-newa, J.n_rows-1), span::all);
	}

	// Extend in y:
	if (extsize_(2)>extsize_(3)) {
		J = cs_lab_lap_arma::padarraySymmetricSimple(J, 0, extsize_(2), 0);
		J = J.tube(span::all, span(0, newb-1));
	} else {
		J = cs_lab_lap_arma::padarraySymmetricSimple(J, 0, extsize_(3), 0);
		J = J.tube(span::all, span(J.n_cols-newb, J.n_cols-1));
	}

	// Extend in z:
	if (extsize_(4) > extsize_(5)) {
		J = cs_lab_lap_arma::padarraySymmetricSimple(J, 0, 0, extsize_(4));
		J = J.slices(0,newc-1);
	} else {
		J = cs_lab_lap_arma::padarraySymmetricSimple(J, 0, 0, extsize_(5));
		J = J.slices(J.n_slices-newc, J.n_slices-1);
	}

	return J;
}

CubeType Gadgetron::ShiftEngine3D::filtering(ColType numerator, ColType denumerator, CubeType i_)
{
	CubeType out(i_.n_rows, i_.n_cols, i_.n_slices, fill::zeros);

	// bsxfun minus
	for (uword k = 0; k < out.n_slices; k++) {
		out.slice(k) = i_.slice(k).each_row()-i_.slice(k).row(0);
	}

	for (uword k = 0; k < out.n_slices; k++) {
		for (uword j = 0; j < out.n_cols; j++) {
			sp::IIR_filt<PixelType, PixelType, PixelType> iir_filt;
			iir_filt.set_coeffs(numerator,denumerator);
			//I_.slice(k).col(j) = iir_filt.filter(I_.slice(k).col(j));
			ColType x = out.slice(k).col(j);
			ColType y = iir_filt.filter(x);
			out.slice(k).col(j) = y;
		}
	}

	// bsxfun plus
	for (uword k = 0; k < out.n_slices; k++) {
		out.slice(k) = out.slice(k).each_row()+i_.slice(k).row(0)*1/accu(denumerator);
	}

	return out;
}

CubeType Gadgetron::ShiftEngine3D::symfilter(PixelType a, PixelType b, CubeType x)
{
	int N = x.n_rows;
	int K1 = x.n_cols;
	int K2 = x.n_slices;

	ColType z0 = cs_lab_lap_arma::roots(ColType({b, a, b}));

	if (abs(z0(0)) < 1) {
		z0 = z0(0);
	} else {
		z0 = z0(1);
	}

	PixelType z0f = z0(0);

	if (abs(z0f) >= 1) {
		throw std::runtime_error("The filter has poles on the unit circle!");
	}

	PixelType A = 1/b/(z0f-1/z0f);

	// One-pole IIR filtering of a symmetrized version of x
	ColType numerator = ColType({1});
	ColType denumerator = ColType({1, -z0f});
	sp::IIR_filt<PixelType, PixelType, PixelType> iir_filt;
	iir_filt.set_coeffs(numerator, denumerator);

	ColType one(1, fill::ones);
	ColType zeros(2*N, fill::zeros);
	ColType onezeros = join_vert(one, zeros);
	ColType z0n = iir_filt.filter(onezeros);
	//ColType z0n = cs_lab_lap_filter::iir(onezeros, numerator, denumerator);
	CubeType t(2*N+1, K1, K2);

	for (uword k = 0; k < x.n_slices; k++) {
		MatType temp = join_vert(x.slice(k), flipud(x.slice(k)));
		t.slice(k) = join_vert(temp, x.slice(k).row(0));
	}

	CubeType out(2*N+1, K1, K2);

	for (uword k = 0; k < out.n_slices; k++) {
		for (uword j = 0; j < out.n_cols; j++) {
			sp::IIR_filt<PixelType, PixelType, PixelType> iir_filt2;
			iir_filt2.set_coeffs(numerator, denumerator);
			ColType x = t.slice(k).col(j);
			ColType y = iir_filt2.filter(x);
			out.slice(k).col(j) = y;
		}
	}

	CubeType a1 = (out.tube(span(out.n_rows-1), span::all) - out.tube(span(0), span::all)) / (1-z0n(z0n.n_elem-1));

	//%**************************************************************************
	//%*****  A naive approach to N-Dim multiplication using For loops	 ******
	//%***				Fine provided K2 is not too large				   ***
	//%**************************************************************************
	CubeType holder(z0n.n_elem, K1, K2, fill::zeros);

	for (int l1 = 0; l1 < K2; l1++) {
		holder.slice(l1) = z0n*a1.slice(l1);
	}

	out = out+holder;
	CubeType out1 = out.tube(span(0,N-1), span::all);
	CubeType out2 = out.tube(span(N, 2*N-1), span::all);

	for (uword k = 0; k < out2.n_slices; k++) {
		out2.slice(k) = flipud(out2.slice(k));
	}

	out = A * (out1+ out2-x);

	return out;
}

bool Gadgetron::ShiftEngine3D::mirt3D_mexinterp(PixelType *Z, PixelType *S, PixelType *T, PixelType *W, PixelType *F, int MN, int nrows, int ncols, int npages, int ndim)
{
	int n, in1, in2, in3, in4, in5, in6, in7, in8;
	double t, s, s1, w, w1, tmp, nan;
	double m1, m2, m3, m4, m5, m6, m7, m8;
	int ndx, nw, Zshift, i, nrowsncols, ft, fs, fw;

	nw = nrows*ncols;
	nrowsncols=nw*npages;
	nan=datum::nan;

	for (n=0; n < MN; n++) {
		t=T[n];
		s=S[n];
		w=W[n];

		ft=(int) floor(t);
		fs=(int) floor(s);
		fw=(int) floor(w);

		if (fs < 1 || s > ncols || ft < 1 || t > nrows || fw < 1 || w > npages) {
			 /* Put nans if outside*/
			for (i = 0; i < ndim; i++) {
				F[n+i*MN]=nan;
			}
		} else {
			ndx = ft + (fs - 1) * nrows + (fw - 1) * nw;

			if (s == ncols) {
				s = s + 1;
				ndx = ndx - nrows;
			}
			s = s - fs;

			if (t == nrows) {
				t = t + 1;
				ndx = ndx - 1;
			}
			t = t - ft;

			if (w == npages) {
				w = w + 1;
				ndx = ndx - nw;
			}
			w = w - fw;

			in1 = ndx - 1;
			in2 = ndx;
			// in3=ndx+nrows-1;
			in3 = in1 + nrows;
			// in4=ndx+nrows;
			in4 = in3 + 1;
			// in5=ndx+nw-1;
			in5 = in1 + nw;
			// in6=ndx+nw;
			in6 = in5 + 1;
			// in7=ndx+nrows+nw-1;
			in7 = in5 + nrows;
			// in8=ndx+nrows+nw;
			in8 = in7 + 1;

			////////////
			s1 = 1 - s;
			w1 = 1 - w;

			tmp = s1 * w1;
			m2 = t * tmp;
			m1 = tmp - m2;

			tmp = s * w1;
			m4 = t * tmp;
			m3 = tmp - m4;

			tmp = s1 * w;
			m6 = t * tmp;
			m5 = tmp - m6;

			tmp = s * w;
			m8 = t * tmp;
			m7 = tmp - m8;

			 for (i = 0; i < ndim; i++) {
				Zshift = i * nrowsncols;
				F[n+i*MN] = Z[in1+Zshift] * m1
							+ Z[in2+Zshift] * m2
							+ Z[in3+Zshift] * m3
							+ Z[in4+Zshift] * m4
							+ Z[in5+Zshift] * m5
							+ Z[in6+Zshift] * m6
							+ Z[in7+Zshift] * m7
							+ Z[in8+Zshift] * m8;
			}
		}
	} // cycle end

	return true;
}

CubeType Gadgetron::ShiftEngine3D::cubicInterp(CubeType x, CubeType y, CubeType z, CubeType I)
{
	int L1 = -2;
	int L2 = 2;

	//Minimum and maximum row index needed in the interpolation formula
	int k0 = floor(x.min()-L2+1);
	int k1 = floor(x.max()-L1);
	int l0 = floor(y.min()-L2+1);
	int l1 = floor(y.max()-L1);
	int m0 = floor(z.min()-L2+1);
	int m1 = floor(z.max()-L1);

	// Convert entries to float to suppress narrowing conversion warning
	ColType offset({
		static_cast<float>(1-k0),
		static_cast<float>(1-l0),
		static_cast<float>(1-m0)
	});

	// Smallest box enclosing the image and the (x,y) positions
	int kk0 = std::min(k0,1);
	int kk1 = std::max(k1,a);
	int ll0 = std::min(l0,1);
	int ll1 = std::max(l1,b);
	int mm0 = std::min(m0,1);
	int mm1 = std::max(m1,c);

	// Indices used in the interpolation formula
	CubeType k=arma::floor(x-L2+1);
	CubeType l=arma::floor(y-L2+1);
	CubeType m=arma::floor(z-L2+1);

	// Convert entries to float to suppress narrowing conversion warning
	ColType extsize = {
		static_cast<float>(1-kk0),
		static_cast<float>(kk1-a),
		static_cast<float>(1-ll0),
		static_cast<float>(ll1-b),
		static_cast<float>(1-mm0),
		static_cast<float>(mm1-c)
	};
	CubeType I0 = ext(image, extsize);
	I0 = I0(span(1-kk0+k0-1, 1-kk0+k1-1), span(1-ll0+l0-1, 1-ll0+l1-1), span(1-mm0+m0-1, 1-mm0+m1-1));

	int a0 = I0.n_rows;
	int b0 = I0.n_cols;

	// now prefiltering
	// along 1st dimension
	CubeType J = symfilter((float)13/(float)21,(float)4/(float)21,image);

	// then 2nd
	J = cs_lab_lap_arma::permuteSimple(J, 213);
	J = symfilter((float)13/(float)21,(float)4/(float)21,J);

	// finnaly 3rd
	J = cs_lab_lap_arma::permuteSimple(J, 321);
	J = symfilter((float)13/(float)21,(float)4/(float)21,J);
	J = cs_lab_lap_arma::permuteSimple(J, 231);

	// Convert entries to float to suppress narrowing conversion warning
	extsize = {
		static_cast<float>(1-kk0),
		static_cast<float>(kk1-a),
		static_cast<float>(1-ll0),
		static_cast<float>(ll1-b),
		static_cast<float>(1-mm0),
		static_cast<float>(mm1-c)
	};

	I0 = ext(J, extsize);
	I0 = I0(span(1-kk0+k0-1, 1-kk0+k1-1), span(1-ll0+l0-1, 1-ll0+l1-1), span(1-mm0+m0-1, 1-mm0+m1-1));

	//% Kernel-based interpolation formula

	//%**** This is probably not the best approach *****
	//%*  The majority of the computation time is here *

	CubeType JJ(x.n_rows, x.n_cols, x.n_slices, fill::zeros);
	ColType h = cs_lab_lap_arma::regspaceSimple(0, L2-L1-1);
	CubeType DK(h.n_elem, h.n_elem, h.n_elem, fill::zeros);

	for (uword k = 0; k < DK.n_slices; k++) {
		DK.slice(k) = repmat(h, 1, DK.n_cols);
	}

	CubeType DL = cs_lab_lap_arma::permuteSimple(DK, 213);
	CubeType DM = cs_lab_lap_arma::permuteSimple(DK, 321);

	for (uword i = 0; i < DK.n_elem; i++) {
		int dk = DK(i);
		int dl = DL(i);
		int dm = DM(i);

		CubeType ind = k+dk+offset(0)+a0*(l+dl+offset(1)-1) + a0*b0*(m+dm+offset(2)-1);
		ind = ind - 1;

		CubeType x1 = 1-abs(x-k-dk);
		CubeType x13 = x1%x1%x1;
		CubeType x2=2-abs(x-k-dk);
		CubeType x23=x2%x2%x2;
		CubeType x2_12 = x2;
		x2_12.transform([](double val) {return (val > 1 && val <=2) ? double(1): double(0);});
		CubeType x2_01 = x2;
		x2_01.transform([](double val) {return (val > 0 && val <=1) ? double(1): double(0);});
		CubeType phi1=((float)1/6*x23-(float)2/3*x13+(float)1/42*x2-(float)2/21*x1)%x2_12+((float)1/6*x23+(float)1/42*x2)%x2_01;

		x1 = 1-abs(y-l-dl);
		x13 = x1%x1%x1;
		x2=2-abs(y-l-dl);
		x23=x2%x2%x2;
		x2_12 = x2;
		x2_12.transform([](double val) {return (val > 1 && val <=2) ? double(1): double(0);});
		x2_01 = x2;
		x2_01.transform([](double val) {return (val > 0 && val <=1) ? double(1): double(0);});
		CubeType phi2=((float)1/6*x23-(float)2/3*x13+(float)1/42*x2-(float)2/21*x1)%x2_12+((float)1/6*x23+(float)1/42*x2)%x2_01;

		x1 = 1-abs(z-m-dm);
		x13 = x1%x1%x1;
		x2=2-abs(z-m-dm);
		x23=x2%x2%x2;
		x2_12 = x2;
		x2_12.transform([](double val) {return (val > 1 && val <=2) ? double(1): double(0);});
		x2_01 = x2;
		x2_01.transform([](double val) {return (val > 0 && val <=1) ? double(1): double(0);});
		CubeType phi3=((float)1/6*x23-(float)2/3*x13+(float)1/42*x2-(float)2/21*x1)%x2_12+((float)1/6*x23+(float)1/42*x2)%x2_01;

		uvec ind_vec = vectorise(conv_to<ucube>::from(ind));

		ColType subI0_vec = I0.elem(ind_vec);

		CubeType subI0(subI0_vec.memptr(), ind.n_rows, ind.n_cols, ind.n_slices);
		JJ = JJ+phi1 % phi2 % phi3 % subI0;
	}

	return JJ;
}

