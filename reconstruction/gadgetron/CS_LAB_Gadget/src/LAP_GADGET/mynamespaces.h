#ifndef MYNAMESPACES_H
#define MYNAMESPACES_H
#include <armadillo>
typedef float PixelType;
typedef arma::cx_float Cx_PixelType;
typedef arma::Cube<PixelType> CubeType;
typedef arma::Mat<PixelType> MatType;
typedef arma::Col<PixelType> ColType;
typedef arma::Col<Cx_PixelType> CxColType;
typedef arma::Row<PixelType> RowType;
typedef arma::uword SizeType;
using namespace arma;
namespace myarma {
ColType roots(ColType c);
CubeType permuteSimple(CubeType &in, int order);
//cube repmatSimple(mat &m, int r, int c, int s);
//vec average(vec I, int K, int s1, int s2, int s3);
CubeType reshapeSimpleVecToCube(ColType v, int r, int c, int s);
ColType regspaceSimple(int start, int end);
ColType padarrayConstantPost(const ColType &in_, int rows, float padval);
ColType padarraySymmetricSimple(const ColType &in_, int rows);
CubeType padarraySymmetricSimple(const CubeType &in_, int rows, int cols, int slices);
}
namespace myfilterITK {
bool x(CubeType &m, ColType kernel);
bool y(CubeType &m, ColType kernel);
bool z(CubeType &m, ColType kernel);
CubeType filter(CubeType &m,  ColType kernel);
CubeType filter(CubeType &m, ColType kernelX, ColType kernelY, ColType kernelZ);
CubeType filterGaussian(CubeType &m, float sigma, ColType kernel);
CubeType filterGaussian(CubeType &m, float sigma, ColType kernelX, ColType kernelY, ColType kernelZ);
}
namespace myfilter {
bool xConv(CubeType &m, ColType kernel);
bool yConv(CubeType &m, ColType kernel);
bool zConv(CubeType &m, ColType kernel);
CubeType conv(const CubeType &m, ColType kernel);
CubeType conv(const CubeType &m, ColType kernelX, ColType kernelY, ColType kernelZ);
CubeType gaussian(const CubeType &m, float sigma, ColType kernel);
CubeType gaussian(const CubeType &m, float sigma, ColType kernelX, ColType kernelY, ColType kernelZ);

//ColType iir(const ColType &c, ColType numerator, ColType denumerator);
//bool xIir(CubeType &m, ColType numerator, ColType denumerator);
//bool yIir(CubeType &m, ColType numerator, ColType denumerator);
//bool zIir(CubeType &m, ColType numerator, ColType denumerator);
//CubeType iir(const CubeType &m, ColType numerator, ColType denumerator);
}
namespace myfunctions {
ColType g(ColType &v, float K);
MatType g(MatType &v, float K);
CubeType g(CubeType &c, float K);
}
#endif // MYNAMESPACES_H

