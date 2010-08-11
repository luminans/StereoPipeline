#ifndef __ASP_MULTIVIEWTK_GEOMETRYOPTIMIZER_H__
#define __ASP_MULTIVIEWTK_GEOMETRYOPTIMIZER_H__

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

#include <asp/MultiviewTK/GeometryOptimizer.h>
#include <asp/MultiviewTK/MatrixImageMath.h>

namespace vw {
namespace multiview {

template <class ImageT>  
struct GeometryOptimizer {
  // Size of the correlation window
  static const int half_win = 8;
  static const int whole_win = half_win * 2 + 1;

  // Size of the smoothing kernel
  static const int half_kern = 2;
  static const int whole_kern = half_kern * 2 + 1;

  // Final size of correlation window, with padding
  static const int half_win_pad = half_win + half_kern;
  static const int whole_win_pad = half_win_pad * 2 + 1;

  typedef double result_type;

  // Domain and result type are planes. In this case, the last entry is
  // the radius in sphereical coords (unlike the stuff in gen_synth_scene)
  typedef Vector4 domain_type;
  typedef Vector4 gradient_type;

  const Vector2 m_lonlat;

  cartography::GeoReference m_georef;

  std::vector<ImageT> m_image_list;
  std::vector<camera::PinholeModel> m_camera_list;

  const unsigned m_num_patches;

  ImageView<double> m_window;

  GeometryOptimizer(int32 x, int32 y, cartography::GeoReference const& georef,
                    std::vector<ImageT> const& image_list, 
                    std::vector<camera::PinholeModel> const& camera_list) :
    m_lonlat(georef.pixel_to_lonlat(Vector2(x, y))), 
    m_georef(georef) , 
    m_image_list(image_list), m_camera_list(camera_list),
    m_num_patches(image_list.size()),
    m_window(generate_gaussian_derivative_kernel(whole_win_pad / 6.0, 0,
                                                 whole_win_pad / 6.0, 0,
                                                 M_PI / 2, whole_win_pad))
  {
    translate_georef(m_georef, Vector2(x - half_win_pad, y - half_win_pad));
  }
     
  result_type operator()(domain_type const& x) const {

    std::vector<ImageView<float32> > patch_list(get_ortho_patches(x));
    ImageView<float32> albedo(whole_win_pad, whole_win_pad);

    BOOST_FOREACH(ImageView<float32> patch, patch_list) {
      albedo += patch;
    }

    albedo /= m_num_patches;

    result_type result = 0;

    for (unsigned i = 0; i < m_num_patches; i++) {
      ImageView<float32> e = patch_list[i] - albedo;
      result += sum_of_pixel_values(crop(m_window * gaussian_filter(e * e, whole_kern / 6.0),
                                         half_kern, half_kern, whole_kern, whole_kern));
      
    }

    return result;
  } 

  gradient_type gradient(domain_type const& x) const {
    static int iter = 0;
    iter++;
    Vector3 e = cartography::lonlat_to_normal(m_lonlat);
    std::cout << x[3] / dot_prod(subvector(x, 0, 3), e) - m_georef.datum().radius(m_lonlat[0], m_lonlat[1]) << std::endl;

    std::vector<ImageView<float32> > patch_list(get_ortho_patches(x));
    ImageView<float32> albedo(whole_win_pad, whole_win_pad);

    BOOST_FOREACH(ImageView<float32> patch, patch_list) {
      albedo += patch;
    }

    albedo /= m_num_patches;

    ImageView<Vector2> grad_albedo = derivative_gaussian_filter(albedo, whole_kern / 6.0, 
                                                                whole_kern / 6.0, M_PI / 2, whole_kern);

    gradient_type result(0, 0, 0, 0);
    for (unsigned i = 0; i < m_num_patches; i++) {
      ImageView<Vector2> grad_l = grad_albedo * gaussian_filter(patch_list[i] - albedo, whole_kern / 6.0);

      result += sum_of_pixel_values(crop(m_window * 
        plane_jacobian_view(m_lonlat, m_georef, m_camera_list[i], x, whole_win_pad, whole_win_pad) * grad_l,
        half_kern, half_kern, whole_win, whole_win)); 

      std::stringstream ss;
      ss << "out/" << iter << "-" << i << ".tif";
      write_image(ss.str(), channel_cast_rescale<uint8>(patch_list[i]));

    }


    return result;
  }

  std::vector<ImageView<float32> > get_ortho_patches(domain_type const& x) const {
    // project onto patch_list 
    std::vector<ImageView<float32> > patch_list(m_num_patches);

    for (unsigned i = 0; i < m_num_patches; i++) {
      boost::shared_ptr<camera::CameraModel> cam(new camera::PinholeModel(m_camera_list[i]));
      patch_list[i] = cartography::orthoproject(plane_dem_view(m_georef, x, 
                                                               whole_win_pad, whole_win_pad),
                                                m_georef, m_image_list[i], cam,
                                                BilinearInterpolation(), ZeroEdgeExtension());
    }

    return patch_list;
  }

  unsigned dimension() const {
    return 4;
  }
};

}} // namespace multiview, vw

#endif
