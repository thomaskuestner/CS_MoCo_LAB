#include "mynamespaces.h"
//#include "sigpack/sigpack.h"

ColType myarma::regspaceSimple(int start, int end){
    return linspace<ColType>(start, end, fabs(start-end)+1);
}

//Own roots method (copied from matlab)
ColType myarma::roots(ColType c){

    //c = c.t();
    int n = c.n_elem;

    uvec inz = find(c);
    if (inz.is_empty()){
        // All elements are zero
        return false;
    }

    // Strip leading zeros and throw away.
    // Strip trailing zeros, but remember them as roots at zero.
    int nnz = inz.n_elem;
    c = c(span(inz(0),inz(nnz-1)));
    int t = inz(nnz-1);
    ColType r;

    // Prevent relatively small leading coefficients from introducing Inf
    // by removing them.
    ColType c_ = c(span(1,c.n_elem-1));
    ColType d = c_/c(0);
    while (d.has_inf()){
        c = c(span(1,c.n_elem-1));
        d = c(span(1, c.n_elem-1))/c(0);
    }

    // Polynomial roots via a companion matrix
    n = c.n_elem;
    if (n > 1){
        MatType a = diagmat(ColType(n-2, fill::ones),-1);
        //mat a = diag(ones(1,n-2,class(c)),-1);
        a.row(0) = -d.t();
        CxColType e = eig_gen(a);
        ColType e_ = conv_to<ColType>::from(e);
        r = join_vert(r, e_);
    }

    return r;
}

//Permute the Matrix dimension, very slow
CubeType myarma::permuteSimple(CubeType &m, int order){
    CubeType temp;
    switch (order) {
    case 123:{
        CubeType temp = m;
        return temp;
        break;}
    case 231:{
        temp.set_size(m.n_cols, m.n_slices, m.n_rows);
        for(int c = 0; c < m.n_cols; ++c){
            for(int r = 0; r < m.n_rows; ++r){
                for(int s = 0; s < m.n_slices; ++s){
                    temp(c, s, r) = m(r, c, s);
                }
            }
        }
        //return temp;
        break;}
    case 312:{
        temp.set_size(m.n_slices, m.n_rows, m.n_cols);
        for(int c = 0; c < m.n_cols; ++c){
            for(int r = 0; r < m.n_rows; ++r){
                for(int s = 0; s < m.n_slices; ++s){
                    temp(s, r, c) = m(r, c, s);
                }
            }
        }
        //return temp;
        break;}
    case 321:{
        temp.set_size(m.n_slices, m.n_cols, m.n_rows);
        for(int c = 0; c < m.n_cols; ++c){
            for(int r = 0; r < m.n_rows; ++r){
                for(int s = 0; s < m.n_slices; ++s){
                    temp(s, c, r) = m(r, c, s);
                }
            }
        }
        //return temp;
        break;}
    case 213:{
        temp.set_size(m.n_cols, m.n_rows, m.n_slices);
        for(int c = 0; c < m.n_cols; ++c){
            for(int r = 0; r < m.n_rows; ++r){
                for(int s = 0; s < m.n_slices; ++s){
                    temp(c, r, s) = m(r, c, s);
                }
            }
        }
        //return temp;
        break;}
    case 132:{
        temp.set_size(m.n_rows, m.n_slices, m.n_cols);
        for(int c = 0; c < m.n_cols; ++c){
            for(int r = 0; r < m.n_rows; ++r){
                for(int s = 0; s < m.n_slices; ++s){
                    temp(r, s, c) = m(r, c, s);
                }
            }
        }
        //return temp;
        break;}
    default:
        break;
    }
    return temp;
}

CubeType myarma::reshapeSimpleVecToCube(ColType v, int r, int c, int s){
    CubeType ret(v.n_elem, 1, 1);
    ret.slice(0).col(0) = v;
    ret.reshape(r, c, s);
    return ret;
}

ColType myarma::padarraySymmetricSimple(const ColType &in_, int rows){
    CubeType temp(in_.n_rows, 1, 1);
    temp.slice(0).col(0) = in_;
    temp = padarraySymmetricSimple(temp, rows, 0, 0);
    return temp.slice(0).col(0);
}

ColType myarma::padarrayConstantPost(const ColType &in_, int elems, float padval){
    ColType ret(in_.n_elem+elems);
    ret(span(0, in_.n_elem-1)) = in_;
    ret(span(in_.n_elem, ret.n_elem-1)).fill(padval);
    return ret;
}

//Symmetric Array Padding
CubeType myarma::padarraySymmetricSimple(const CubeType &in_, int rows, int cols, int slices){
    int oldrows = in_.n_rows;
    int oldcols = in_.n_cols;
    int oldslices = in_.n_slices;

    int newrows = in_.n_rows+2*rows;
    int newcols = in_.n_cols+2*cols;
    int newslices = in_.n_slices+2*slices;

    CubeType ret(newrows, newcols, newslices, fill::zeros);
    ret.subcube(span(rows, rows+oldrows-1), span(cols, cols+oldcols-1), span(slices, slices+oldslices-1)) = in_;

    //Pad above and below
    for(int p = 0; p < rows; p++){
        ret.tube(span(rows-p-1), span::all) = ret.tube(span(rows+(p % oldrows)), span::all);

        ret.tube(span((oldrows+rows)+p), span::all) = ret.tube(span((oldrows+rows)-1-(p % oldrows)), span::all);
    }

    //Pad left and right
    for(int p = 0; p < cols; p++){
        ret.tube(span::all, span(cols-p-1)) = ret.tube(span::all, span(cols+(p % oldcols)));

        ret.tube(span::all, span((oldcols+cols)+p)) = ret.tube(span::all, span((oldcols+cols)-1-(p % oldcols)));
    }

    //Pad front and back
    for(int p = 0; p < slices; p++){
        ret.slice(slices-p-1) = ret.slice(slices+(p % oldslices));

        ret.slice((oldslices+slices)+p) = ret.slice((oldslices+slices)-1-(p % oldslices));
    }
    return ret;
}

//Slow filter. 1D Kernel x Dimension (col)
bool myfilter::xConv(CubeType &m, ColType kernel){
    int orig_rows = m.n_rows;
    int u = floor(kernel.n_elem/2);
    CubeType m_ = myarma::padarraySymmetricSimple(m, u, 0, 0);

//    m_.each_slice([kernel](MatType &slice){slice.each_col([kernel](ColType &c){c = conv(c, kernel, "same");});});

//    m = m_.subcube(span(u, u+orig_rows-1), span::all, span::all);


    #pragma omp parallel for
    for(unsigned int k = 0; k < m_.n_slices; k++){
        for(unsigned int j = 0; j < m_.n_cols; j++){
            ColType v = m_.slice(k).col(j);
            v = arma::conv(v, kernel, "same");
            m.slice(k).col(j) = v.subvec(span(u, u+orig_rows-1));
        }
    }
    return true;
}

//Slow filter. 1D Kernel y Dimension (row)
bool myfilter::yConv(CubeType &m, ColType kernel){
    int orig_cols = m.n_cols;
    int u = floor(kernel.n_elem/2);
    CubeType m_ = myarma::padarraySymmetricSimple(m, 0, u, 0);

//   m_.each_slice([kernel](MatType &slice){slice.each_row([kernel](RowType &c){c = conv(c, kernel, "same");});});

//   m = m_.subcube(span::all, span(u, u+orig_cols-1), span::all);
    #pragma omp parallel for
    for(unsigned int k = 0; k < m_.n_slices; k++){
        for(unsigned int i = 0; i < m_.n_rows; i++){
            RowType rv = m_.slice(k).row(i);
            ColType v(rv.memptr(), rv.n_elem);
            v = arma::conv(v, kernel, "same");
            rv = reshape(v, 1, v.n_elem);
            m.slice(k).row(i) = rv.subvec(span(u, u+orig_cols-1));
        }
    }
    return true;
}

//Slow filter. 1D Kernel z Dimension (tube)
bool myfilter::zConv(CubeType &m, ColType kernel){
    int orig_slices = m.n_slices;

    int u = floor(kernel.n_elem/2);
    CubeType m_ = myarma::padarraySymmetricSimple(m, 0, 0, u);
    #pragma omp parallel for
    for(unsigned int i = 0; i < m_.n_rows; i++){
        for(unsigned int j = 0; j < m_.n_cols; j++){
            CubeType tv = m_.tube(i,j);
            ColType v = reshape(tv, tv.n_elem, 1, 1);
            v = arma::conv(v, kernel, "same");
            CubeType tv_(v.memptr(), 1,1,v.n_elem);
            m.tube(i,j) = tv_.subcube(span::all, span::all, span(u, u+orig_slices-1));
        }
    }
    return true;
}

//Slow filter combines all 3 Dims
CubeType myfilter::conv(const CubeType &m, ColType kernel){
    return myfilter::conv(m, kernel, kernel, kernel);
}

//Slow filter combines all 3 Dims
CubeType myfilter::conv(const CubeType &m, ColType kernelX, ColType kernelY, ColType kernelZ){
    CubeType ret = m;
    myfilter::xConv(ret, kernelX);
    myfilter::yConv(ret, kernelY);
    myfilter::zConv(ret, kernelZ);
    return ret;
}

//Slow filter all Dimension with linear Gaussian Kernel
CubeType myfilter::gaussian(const CubeType &m, float sigma, ColType kernel){
    ColType F = exp(-kernel % kernel / 2/ pow(sigma,2));
    F = F / accu(F);
    CubeType ret = myfilter::conv(m, F);
    return ret;
}


CubeType myfilter::gaussian(const CubeType &m, float sigma, ColType kernelX, ColType kernelY, ColType kernelZ){
    ColType Fx = exp(-kernelX % kernelX / 2/ pow(sigma,2));
    Fx = Fx / accu(Fx);
    ColType Fy = exp(-kernelY % kernelY / 2/ pow(sigma,2));
    Fy = Fy / accu(Fy);
    ColType Fz = exp(-kernelZ % kernelZ / 2/ pow(sigma,2));
    Fz = Fz / accu(Fz);
    return myfilter::conv(m, Fx, Fy, Fz);
}

//ColType myfilter::iir(const ColType &c, ColType numerator, ColType denumerator){
//    sp::IIR_filt<PixelType, PixelType, PixelType> iir_filt;
//    iir_filt.set_coeffs(numerator, denumerator);
//    return iir_filt.filter(c);
//}

//bool myfilter::xIir(CubeType &m, ColType numerator, ColType denumerator){


//    for(arma::uword k = 0; k < m.n_slices; k++){
//        for(arma::uword j = 0; j < m.n_cols; j++){
//            sp::IIR_filt<PixelType, PixelType, PixelType> iir_filt;
//            iir_filt.set_coeffs(numerator, denumerator);
//            ColType v = m.slice(k).col(j);
//            v = iir_filt.filter(v);
//            m.slice(k).col(j) = v;
//        }
//    }
//}

//bool myfilter::yIir(CubeType &m, ColType numerator, ColType denumerator){


//    for(arma::uword k = 0; k < m.n_slices; k++){
//        for(arma::uword i = 0; i < m.n_rows; i++){
//            sp::IIR_filt<PixelType, PixelType, PixelType> iir_filt;
//            iir_filt.set_coeffs(numerator, denumerator);
//            RowType rv = m.slice(k).row(i);
//            ColType v(rv.memptr(), rv.n_elem);
//            v = iir_filt.filter(v);
//            rv = reshape(v, 1, v.n_elem);
//            m.slice(k).row(i) = rv;
//        }
//    }
//}

//bool myfilter::zIir(CubeType &m, ColType numerator, ColType denumerator){


//    for(arma::uword i = 0; i < m.n_rows; i++){
//        for(arma::uword j = 0; j < m.n_cols; j++){
//            sp::IIR_filt<PixelType, PixelType, PixelType> iir_filt;
//            iir_filt.set_coeffs(numerator, denumerator);
//            CubeType tv = m.tube(i,j);
//            ColType v = reshape(tv, tv.n_elem, 1, 1);
//            v = iir_filt.filter(v);
//            CubeType tv_(v.memptr(), 1, 1, v.n_elem);
//            m.tube(i,j) = tv_;
//        }
//    }

//}

//CubeType myfilter::iir(const CubeType &m, ColType numerator, ColType denumerator){
//    CubeType ret = m;
//    myfilter::xIir(ret, numerator, denumerator);
//    myfilter::yIir(ret, numerator, denumerator);
//    myfilter::zIir(ret, numerator, denumerator);
//    return ret;
//}




//g Functions for Filter Basis
ColType myfunctions::g(ColType &v, float K){
    float s = (K+2)/4;
    ColType e = -v % v / 2 / pow(s, 2);
    return arma::exp(e);
}

MatType myfunctions::g(MatType &m, float K){
    float s = (K+2)/4;
    MatType e = -m % m / 2 / pow(s, 2);
    return arma::exp(e);
}

CubeType myfunctions::g(CubeType &c, float K){
    float s = (K+2)/4;
    CubeType e = -c % c / 2 / pow(s, 2);
    return arma::exp(e);
}




