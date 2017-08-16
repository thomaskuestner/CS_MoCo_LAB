#ifndef LAP3D_H
#define LAP3D_H
#include <armadillo>
#include "mynamespaces.h"
using namespace arma;
class LAP3D{
private:

    //Private member variables


    int m, n, p, numel;
    CubeType I1, I_input, I2, I2_shifted;
    ColType hGaussian, FilterSizes;
    field<CubeType> u_holder;
    CubeType mask[3];

    //Private member functions
    field<CubeType> estimateOpticalFlow3D( CubeType &I1_,  CubeType &I2_, int K_);
    field<CubeType> estimateOpticalFlow3DKSpace( CubeType &I1_k_, CubeType &I2_k_, int sampleStepSize, int K_);

    ColType average(const ColType &I_, int K_);
    field<CubeType> cleanOF3D(field<CubeType> &u_est_);
    bool cleaningProcedure();

public:
    LAP3D();
    LAP3D(const CubeType &I1_, const CubeType &I2_, int levelMin_, int levelMax_);

    field<CubeType> exec();
};

#endif // LAP3D_H
