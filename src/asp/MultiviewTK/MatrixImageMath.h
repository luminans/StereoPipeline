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

struct PlaneJacobianFunctor {
  typedef Matrix<double, 4, 2> result_type;
  cartography::GeoReference m_georef;
  Vector4 m_plane;
  PlaneJacobianFunctor(cartography::GeoReference const& georef, Vector4 const& plane) :
    m_georef(georef), m_plane(plane) {}
  result_type operator()(double i, double j, int32 /*p*/) const {
    return result_type();
  }
};

} // namespace vw

#endif 
