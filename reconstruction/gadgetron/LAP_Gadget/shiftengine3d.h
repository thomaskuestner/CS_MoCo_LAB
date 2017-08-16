#ifndef SHIFTENGINE3D_H
#define SHIFTENGINE3D_H

#include "mynamespaces.h"

using namespace arma;

class ShiftEngine3D{
private:
    CubeType image;
    CubeType ux;
    CubeType uy;
    CubeType uz;

    CubeType x;
    CubeType y;
    CubeType z;

    int a;
    int b;
    int c;

    //Private member functions
    CubeType ext(CubeType i_, ColType extsize_);
    CubeType filtering(ColType numerator, ColType denumerator, CubeType i_);
    CubeType symfilter(PixelType a, PixelType b, CubeType x);
    bool mirt3D_mexinterp(PixelType* Z, PixelType* S, PixelType* T, PixelType* W, PixelType* F,
                            int	MN, int nrows, int ncols, int npages, int ndim);
    CubeType cubicInterp(CubeType x, CubeType y, CubeType z, CubeType I);

public:
    ShiftEngine3D();
    ShiftEngine3D(CubeType i_, CubeType ux_, CubeType uy_, CubeType uz_);
    CubeType execLinShift();
    CubeType execCubicShift();
};

#endif // SHIFTENGINE3D_H
