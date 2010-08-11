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
/*
template <class SrcT>
SeparableConvolutionView<SrcT, Vector<typename DefaultKernelT<typename SrcT::pixel_type>::type, 2> >
derivative_gaussian_filter(ImageViewBase<SrcT> src, double x_sigma, double y_sigma, double angle) {
  std::vector<Vector<typename DefaultKernelT<typename SrcT::pixel_type>::type, 2> > x_kernel, y_kernel;


} */

template <class SrcT>
ImageView<Vector2>
derivative_gaussian_filter(ImageViewBase<SrcT> const& src, double x_sigma, double y_sigma, double angle, int32 size) {
  ImageView<Vector2> kern = generate_gaussian_derivative_kernel(x_sigma, y_sigma, angle, size);
  ImageView<Vector2> result(src.impl().cols(), src.impl().rows());
  select_channel(result, 0) = convolution_filter(src.impl(), select_channel(kern, 0));
  select_channel(result, 1) = convolution_filter(src.impl(), select_channel(kern, 1));
  return result;
}

struct PlaneJacobianFunctor {
  typedef Matrix<double, 4, 2> result_type;
  Vector2 m_center_lonlat;
  cartography::GeoReference m_georef;
  camera::PinholeModel const& m_camera;
  Vector4 m_plane;

  PlaneJacobianFunctor(Vector2 const& center_lonlat,
                       cartography::GeoReference const& georef,
                       camera::PinholeModel const& camera, 
                       Vector4 const& plane) :
    m_center_lonlat(center_lonlat), m_georef(georef), 
    m_camera(camera), m_plane(plane) {}

  result_type operator()(double i, double j, int32 /*p*/) const {
    Matrix2x2 S_inv = inverse(submatrix(m_georef.transform(), 0, 0, 2, 2));
    Vector3 n = subvector(m_plane, 0, 3);
    Vector3 e = cartography::lonlat_to_normal(m_center_lonlat);
    Vector3 x = m_plane[3] * e;
    Vector2 curr_lonlat = m_georef.pixel_to_lonlat(Vector2(i, j));
    Vector3 e_ = cartography::lonlat_to_normal(curr_lonlat);
    Vector3 x_ = m_plane[3] / dot_prod(n, e_) * e_;

    Matrix<double, 3, 4> P(m_camera.camera_matrix());
    
    Matrix3x3 Q[3];
    Vector3 p[3];
    Vector3 u, u_;
    for (unsigned i = 0; i < 3; i++) {
      p[i] = subvector(select_row(P, i), 0, 3);
      Q[i] = x_ * transpose(p[i]) + P(i, 3) * identity_matrix(3);
      u[i] = dot_prod(n, Q[i] * e);
      u_[i] = dot_prod(n, Q[i] * e_);
    }

    Matrix3x3 A[2];
    Matrix<double, 2, 4> B;
    Matrix<double, 2, 3> C_half;
    for (unsigned i = 0; i < 2; i++) {
      A[i] = u_[2] * Q[i] - u_[i] * Q[2];

      Matrix<double, 3, 4> tmp;
      submatrix(tmp, 0, 0, 3, 3) = A[i];
      select_col(tmp, 3) = u[2] * p[i] - u[i] * p[2];
      select_row(B, i) = transpose(tmp) * e_;

      select_row(C_half, i) = transpose(A[i]) * n;
    }

    Vector2 cos_curr_lonlat = curr_lonlat * M_PI / 180;
    Vector2 sin_curr_lonlat = curr_lonlat * M_PI / 180;
    for (unsigned i = 0; i < 2; i++) {
      cos_curr_lonlat[i] = cos(cos_curr_lonlat[i]);
      sin_curr_lonlat[i] = sin(sin_curr_lonlat[i]);
    }

    Matrix<double, 3, 2> de_d0;
    de_d0(0, 0) = -cos_curr_lonlat[1] * sin_curr_lonlat[0];
    de_d0(1, 0) = cos_curr_lonlat[1] * cos_curr_lonlat[0];
    de_d0(2, 0) = 0;

    de_d0(0, 1) = -sin_curr_lonlat[1] * cos_curr_lonlat[0];
    de_d0(1, 1) = -sin_curr_lonlat[1] * sin_curr_lonlat[0];
    de_d0(2, 1) = cos_curr_lonlat[1];

    return transpose(-S_inv * inverse(C_half * de_d0) * B);
  }
};

inline PerPixelIndexView<PlaneJacobianFunctor>
plane_jacobian_view(Vector2 const& center_lonlat,
                    cartography::GeoReference const& georef,
                    camera::PinholeModel const& camera,
                    Vector4 const& plane, int32 cols, int32 rows) {
  typedef PerPixelIndexView<PlaneJacobianFunctor> result_type;
  return result_type(PlaneJacobianFunctor(center_lonlat, georef, camera, plane),
                     cols, rows, 1);
}

} // namespace vw

#endif 
