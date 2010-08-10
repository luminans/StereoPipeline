#ifndef __ASP_MULTIVIEWTK_MATRIXIMAGEMATH_H__
#define __ASP_MULTIVIEWTK_MATRIXIMAGEMATH_H__

namespace vw {

template <class MatrixT, class VectorT, int i, int j>
struct ProductTypeSpecialization<Matrix<MatrixT, i, j>, Vector<VectorT, j> > {
  typedef Vector<typename ProductType<MatrixT, VectorT>::type, i> type;
};

inline ImageView<Vector2> 
generate_gaussian_derivative_kernel(double x_sigma, double y_sigma, 
                                    double angle, int32 size)
{
  ImageView<Vector2> result(size, size);
  select_channel(result, 0) = generate_gaussian_derivative_kernel(x_sigma, 1, 
                                                                  y_sigma, 0, 
                                                                  angle, size);
  select_channel(result, 1) = generate_gaussian_derivative_kernel(x_sigma, 0, 
                                                                  y_sigma, 1, 
                                                                  angle, size);
  return result;
}

} // namespace vw

#endif 
