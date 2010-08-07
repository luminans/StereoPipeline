#ifndef __ASP__MULTIVIEWTK__GEOMETRYOPTIMIZER_H__
#define __ASP__MULTIVIEWTK__GEOMETRYOPTIMIZER_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Manipulation.h>

#include <boost/foreach.hpp>

#include <vw/FileIO.h>

#include <vw/Math/ConjugateGradient.h>
#include <vw/Math/Vector.h>

#include <asp/MultiviewTK/multiview.h>

namespace vw {
namespace multiview {

template <class ImageT>  
struct GeometryOptimizer {
  static const int half_kern = 8;
  static const int whole_kern = half_kern * 2 + 1;

  typedef float32 result_type;

  // Domain and result type are planes. In this case, the last entry is
  // the radius in sphereical coords (unlike the stuff in gen_synth_scene)
  typedef Vector4 domain_type;
  typedef Vector4 gradient_type;

  const Vector2 m_lonlat;

  cartography::GeoReference m_georef;

  std::vector<ImageT> m_image_list;
  std::vector<camera::PinholeModel> m_camera_list;

  const unsigned m_num_patches;

  ImageView<float32> m_kernel, m_kernel_deriv_x, m_kernel_deriv_y;

  GeometryOptimizer(int32 x, int32 y, cartography::GeoReference const& georef,
                    std::vector<ImageT> const& image_list, 
                    std::vector<camera::PinholeModel> const& camera_list) :
    m_lonlat(georef.pixel_to_lonlat(Vector2(x, y))), 
    m_georef(get_crop_georef(georef, BBox2i(x - half_kern, y - half_kern, whole_kern, whole_kern))), 
    m_image_list(image_list), m_camera_list(camera_list),
    m_num_patches(image_list.size()),
    m_kernel(generate_gaussian_derivative_kernel(whole_kern / 6.0, 0,
                                                 whole_kern / 6.0, 0,
                                                 M_PI / 2, whole_kern)),
    m_kernel_deriv_x(generate_gaussian_derivative_kernel(whole_kern / 6.0, 1, 
                                                         whole_kern / 6.0, 0, 
                                                         M_PI / 2, whole_kern)),
    m_kernel_deriv_y(generate_gaussian_derivative_kernel(whole_kern / 6.0, 0, 
                                                         whole_kern / 6.0, 1, 
                                                         M_PI / 2, whole_kern))
  {}
     
  result_type operator()(domain_type const& x) const {
    float32 sum = 0;

    std::vector<ImageView<float32> > cost_list = get_cost_patches(get_ortho_patches(x));

    BOOST_FOREACH(ImageView<float32> cost, cost_list) {
      sum += sum_of_pixel_values(m_kernel * cost);
    }

    return sum;
  } 

  gradient_type gradient(domain_type const& x) const {
    std::vector<ImageView<float32> > patch_list = get_ortho_patches(x);
    std::vector<ImageView<float32> > cost_list = get_cost_patches(patch_list);

    // This can be precomputed
    Matrix<double, 3, 2> de_d0;
    de_d0(0, 0) = -cos(m_lonlat[1]) * sin(m_lonlat[0]);
    de_d0(1, 0) = cos(m_lonlat[1]) * cos(m_lonlat[0]);
    de_d0(2, 0) = 0;

    de_d0(0, 1) = -sin(m_lonlat[1]) * cos(m_lonlat[0]);
    de_d0(1, 1) = -sin(m_lonlat[1]) * sin(m_lonlat[0]);
    de_d0(2, 1) = cos(m_lonlat[1]);

    // This can also can be precomputed.
    // TODO: Remove inverses
    Matrix2x2 S_inv = inverse(submatrix(m_georef.transform(), 0, 0, 2, 2));

    Vector3 e = cartography::lonlat_to_normal(m_lonlat);
    Vector3 x_ = x[3] * e;
    Vector3 n = math::SubVector<Vector4 const>(x, 0, 3);

    Vector4 result(0, 0, 0, 0);
    for (int i = 0; i < m_num_patches; i++) {
      Matrix<double, 3, 4> P(m_camera_list[i].camera_matrix());

      Vector3 p1(P(0, 0), P(0, 1), P(0, 2));
      Vector3 p2(P(1, 0), P(1, 1), P(1, 2));
      Vector3 p3(P(2, 0), P(2, 1), P(2, 2));

      Matrix3x3 Q1 = x_ * transpose(p1) + P(0, 3) * identity_matrix(3);
      Matrix3x3 Q2 = x_ * transpose(p2) + P(1, 3) * identity_matrix(3);
      Matrix3x3 Q3 = x_ * transpose(p3) + P(2, 3) * identity_matrix(3);

      double u1 = dot_prod(n, Q1 * e);
      double u2 = dot_prod(n, Q2 * e);
      double u3 = dot_prod(n, Q3 * e);

      Matrix3x3 A1 = u3 * Q1 - u1 * Q3;
      Matrix3x3 A2 = u3 * Q2 - u2 * Q3;

      Matrix<double, 2, 3> C_half;
      select_row(C_half, 0) = transpose(A1) * n;
      select_row(C_half, 1) = transpose(A2) * n;

      Matrix2x2 C_inv = inverse(C_half * de_d0);

      Matrix<double, 2, 3> tmp;
      select_row(tmp, 0) = select_col(A1, 0);
      select_row(tmp, 1) = select_col(A2, 0);
      Vector2 dz_dnx = -S_inv * C_inv * tmp * e;
      select_row(tmp, 0) = select_col(A1, 1);
      select_row(tmp, 1) = select_col(A2, 1);
      Vector2 dz_dny = -S_inv * C_inv * tmp * e;
      select_row(tmp, 0) = select_col(A1, 2);
      select_row(tmp, 1) = select_col(A2, 2);
      Vector2 dz_dnz = -S_inv * C_inv * tmp * e;
      select_row(tmp, 0) = u3 * p1 - u1 * p3;
      select_row(tmp, 1) = u3 * p2 - u2 * p3;
      Vector2 dz_dr = -S_inv * C_inv * tmp * e * dot_prod(n, e);

      Vector2 cost_sum(sum_of_pixel_values(m_kernel_deriv_x * cost_list[i]),
                       sum_of_pixel_values(m_kernel_deriv_y * cost_list[i]));

      result +=  Vector4(dot_prod(cost_sum, dz_dnx),
                         dot_prod(cost_sum, dz_dny),
                         dot_prod(cost_sum, dz_dnz),
                         dot_prod(cost_sum, dz_dr));
    }

    return result;
  }

  std::vector<ImageView<float32> > get_ortho_patches(domain_type const& x) const {
    // convert to gen_synth_scene plane type
    double radius = x[3] + m_georef.datum().radius(m_lonlat[0], m_lonlat[1]);
    Vector3 pt = cartography::lon_lat_radius_to_xyz(Vector3(m_lonlat[0], m_lonlat[1], radius));
    Vector3 normal = math::SubVector<Vector4 const>(x, 0, 3);
    Vector4 plane(normal[0], normal[1], normal[2], dot_prod(normal, pt));
   
    // project onto patch_list 
    std::vector<ImageView<float32> > patch_list(m_num_patches);

    for (unsigned i = 0; i < m_num_patches; i++) {
      boost::shared_ptr<camera::CameraModel> cam(new camera::PinholeModel(m_camera_list[i]));
      patch_list[i] = cartography::orthoproject(plane_dem_view(m_georef, plane, 
                                                               whole_kern, whole_kern),
                                                m_georef, m_image_list[i], cam,
                                                BilinearInterpolation(), ZeroEdgeExtension());
    }

    return patch_list;
  }

  std::vector<ImageView<float32> > 
  get_cost_patches(std::vector<ImageView<float32> > const& patch_list) const {
    std::vector<ImageView<float32> > result(m_num_patches);

    ImageView<float32> albedo(patch_list[0].cols(), patch_list[0].rows());

    BOOST_FOREACH(ImageView<float32> patch, patch_list) {
      albedo += patch;
    }

    albedo /= patch_list.size();

    for (int i = 0; i < m_num_patches; i++) {
      result[i] = patch_list[i] * log(albedo / patch_list[i]) + patch_list[i] - albedo;
    }

    return result;
  }

  unsigned dimension() const {
    return 4;
  }
};

}} // namespace multiview, vw

#endif
