#include "lap3d.h"
#include "shiftengine3d.h"
#include "gaussianfilterbasis.h"

LAP3D::LAP3D(){

}

LAP3D::LAP3D(const CubeType &I1_, const CubeType &I2_, int levelMin_, int levelMax_){

    //Check if images have same size
    if(I1_.n_rows != I2_.n_rows || I1_.n_rows != I2_.n_rows || I1_.n_rows != I2_.n_rows){
        throw std::runtime_error("Images must have same size for LAP3D Algorithm");
    }


    I1 = I1_;
    I2 = I2_;

    //Image dimensions
    m = I1.n_rows;
    n = I1.n_cols;
    p = I1.n_slices;
    numel = I1.n_elem;

    //Gaussian Smoothing Kernel
    hGaussian = "0.0545 0.2442 0.4026 0.2442 0.0545";

    //Calculate Filtersizes
    FilterSizes = linspace<ColType>(levelMax_, levelMin_, fabs(levelMax_-levelMin_)+1);
    FilterSizes = exp2(FilterSizes);

    //Build holders for the estimated flow
    u_holder.set_size(3);
    for(int i = 0; i < 3; i++){
        u_holder(i).set_size(m,n,p);
        u_holder(i).fill(0);
    }
}

void LAP3D::setMovingImage(const CubeType &I2_){
	I2 = I2_;
}

field<CubeType> LAP3D::exec(){

    //Prefilter I1
    I1 = I1 - myfilter::conv(I1, hGaussian);



    for(int l = 0; l < FilterSizes.n_elem; l++){
        std::cout << "Level " << l+1 << "/" << FilterSizes.n_elem << endl;

        //Current Filtersize
        int K = FilterSizes(l);

        //Warp Image
        wall_clock warping_time;
        warping_time.tic();
        if(l==0){
            I2_shifted = I2-myfilter::conv(I2,hGaussian);
        }else{
            //Create Shiftengine
            ShiftEngine3D shifter(I2, -u_holder(0), -u_holder(1), -u_holder(2));
            I2_shifted = shifter.execLinShift();
            I2_shifted = I2_shifted - myfilter::conv(I2_shifted, hGaussian);
        }
        cout << "Filter size= " << 2*K+1 << ", Warping time = " << warping_time.toc() << endl;

        //Estimate the Optical Flow
        wall_clock lap_time;
        lap_time.tic();
        field<CubeType> u_est = estimateOpticalFlow3D(I1, I2_shifted, K);
        cout << "Filter size= " << 2*K+1 << ", LAP Algorithm time = " << lap_time.toc() << endl;

        //Postprocessing
        wall_clock post_processing_time;
        post_processing_time.tic();
        u_est = cleanOF3D(u_est); //Not completely implemented
        int R = round(2*K);
        ColType k1 = linspace<ColType>(-R, R, fabs(-R-R)+1);

        //Smoothing Filter for optical Flow
        for(int i = 0; i < 3; i++){
            u_est(i) = myfilter::gaussian(u_est(i), R, k1);
        }
        cout << "Filter size= " << 2*K+1 << ", Post-Processing time = " << post_processing_time.toc() << endl;

        //Add estimatet Flow to overall Flow
        for(int i = 0; i < 3; i++){
            u_holder(i) = u_holder(i) + u_est(i);
        }
        if(K<=2){
            //Cleaing procedure
        }
    }

//    //Warp the image with final estimated flow
//    ShiftEngine3D finalShifting(I_input, u_holder[0], u_holder[1], u_holder[2]);
//    std::cout << "doing final shifting" << std::endl;
//    wall_clock final_shifting_time;
//    final_shifting_time.tic();
//    CubeType I_hat_fast = finalShifting.execCubicShift();
//    cout << "Time taken fast-adding u on image I1: " << final_shifting_time.toc() << endl;
return u_holder;
}

field<CubeType> LAP3D::estimateOpticalFlow3DKSpace(CubeType &I1_k_, CubeType &I2_k_, int sampleStepSize, int K_){
    int N = 4;
    //Create the GaussianFilterBasis with current Filtersize K
    GaussianFilterBasis mBasis(K_);
    int K1 = 2*K+1; //Blocksize


}



field<CubeType> LAP3D::estimateOpticalFlow3D( CubeType &I1_,  CubeType &I2_, int K_){
    int N = 4;
    //Create the GaussianFilterBasis with current Filtersize K
    GaussianFilterBasis mBasis(K_);
    //Filter Images, very slow, to be improved
    CubeType A(m,n,p);
    CubeType B(m,n,p);
    MatType II(m*n*p, N);
    for(int i = 0; i < N; i++){
        switch (i) {
        case 0:
            A = myfilter::conv(I1_, mBasis.getG());
            B = myfilter::conv(I2_, mBasis.getGi());


            break;
        case 1:
            A = myfilter::conv(I1_, mBasis.getGd(), mBasis.getG(), mBasis.getG());
            B = myfilter::conv(I2_, mBasis.getGdi(), mBasis.getGi(), mBasis.getGi());
            break;
        case 2:
            A = myfilter::conv(I1_, mBasis.getG(), mBasis.getGd(), mBasis.getG());
            B = myfilter::conv(I2_, mBasis.getGi(), mBasis.getGdi(), mBasis.getGi());
            break;
        case 3:
            A = myfilter::conv(I1_, mBasis.getG(), mBasis.getG(), mBasis.getGd());
            B = myfilter::conv(I2_, mBasis.getGi(), mBasis.getGi(), mBasis.getGdi());
            break;
        default:
            break;
        }
        II.col(i) = vectorise(A-B);
    }
    ColType J = II.col(0);
    //matrices needed in linear system
    A.set_size(m*n*p, N-1, N-1);
    MatType b(m*n*p, N-1);
    for(int k = 0; k < N-1; k++){
        for(int l = k; l < N-1; l++){
            A.slice(l).col(k) = average(II.col(k+1)%II.col(l+1), K_);
            A.slice(k).col(l) = A.slice(l).col(k);
        }
        b.col(k) =average(II.col(k+1) % J, K_);
    }

    MatType coeffs(m*n*p, N-1);
    ColType c(m*n*p);
    for(int k = 0; k < N-1; k++){
        for(int l = k+1; l < N-1; l++){
            c = A.slice(k).col(l) / A.slice(k).col(k);
            for(int o = k+1; o < N-1; o++){
                A.slice(o).col(l) = A.slice(o).col(l) - c % A.slice(o).col(k);
            }
            A.slice(k).col(l).fill(0);
            b.col(l) = b.col(l)-c % b.col(k);
        }
    }
    for(int k = N-2; k >=0; k--){
        coeffs.col(k) = b.col(k);
        for(int m = k+1; m < N-1; m++){
            coeffs.col(k) = coeffs.col(k) - A.slice(m).col(k) % coeffs.col(m);
        }
        coeffs.col(k) = coeffs.col(k) / A.slice(k).col(k);
    }

    MatType coeffs_ = coeffs;
    coeffs.set_size(coeffs.n_rows, coeffs.n_cols+1);
    coeffs.col(0).fill(1);
    coeffs.cols(1, coeffs.n_cols-1) = coeffs_;

    //calculation
    ColType u1(I1.n_elem, fill::zeros);
    ColType u11(I1.n_elem, fill::zeros);
    ColType u2(I1.n_elem, fill::zeros);
    ColType u22(I1.n_elem, fill::zeros);
    ColType u3(I1.n_elem, fill::zeros);
    ColType u33(I1.n_elem, fill::zeros);
    for(int i = 0; i < N; i++){
        u1 = u1 - as_scalar(mBasis.getBasis().col(i).t() * vectorise(mBasis.getK())) * coeffs.col(i);
        u11= u11 + accu(mBasis.getBasis().col(i)) * coeffs.col(i);

        u2 = u2 - as_scalar(mBasis.getBasis().col(i).t() * vectorise(mBasis.getL())) * coeffs.col(i);
        u22 = u22 + accu(mBasis.getBasis().col(i)) * coeffs.col(i);

        u3 = u3 - as_scalar(mBasis.getBasis().col(i).t() * vectorise(mBasis.getP())) * coeffs.col(i);
        u33 = u33 + accu(mBasis.getBasis().col(i)) * coeffs.col(i);
    }
    MatType u_hold(u1.n_rows, 3);
    u_hold.col(0) = 2*u1/u11;
    u_hold.col(1) = 2*u2/u22;
    u_hold.col(2) = 2*u3/u33;
    u_hold.rows(find(arma::sqrt(arma::sum(arma::pow(arma::abs(u_hold),2),1))>K_)).fill(datum::nan);
    field<CubeType> u_est(3);
    u_est(0) = myarma::reshapeSimpleVecToCube(u_hold.col(0), m, n, p);
    u_est(1) = myarma::reshapeSimpleVecToCube(u_hold.col(1), m, n, p);
    u_est(2) = myarma::reshapeSimpleVecToCube(u_hold.col(2), m, n, p);

    return u_est;
}

field<CubeType> LAP3D::cleanOF3D(field<CubeType> &u_est){
    //find nans in flow and replace with zeros
    for(int i = 0; i < 3; i++){
        u_est(i).elem(find_nonfinite(u_est(i))).zeros();
    }
    ColType kernelX = "0.0000 0.0002 0.2576 0.0002 0.0000";
    ColType kernelY = "0.0000 0.0010 1.4516 0.0010 0.0000";
    ColType kernelZ = "0.0000 0.0008 1.2258 0.0008 0.0000";
    return u_est;
}

bool LAP3D::cleaningProcedure(){
    for(int i = 0; i < 3; i++){

    }
}


ColType LAP3D::average(const ColType &I_, int K_){
    K_=K_*2+1;
    CubeType II(I_.n_elem, 1, 1);
    II.slice(0).col(0) = I_;
    II.reshape(m, n, p);

    ColType kernel(K_, fill::ones);
    CubeType J_2 = myfilter::conv(II, kernel);
    return vectorise(J_2);
}


