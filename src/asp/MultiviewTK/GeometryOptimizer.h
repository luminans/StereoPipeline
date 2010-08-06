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

#include <asp/MultiviewTK/gen_synth_scene.h>
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

  GeometryOptimizer(int32 x, int32 y, cartography::GeoReference const& georef,
                    std::vector<ImageT> const& image_list, 
                    std::vector<camera::PinholeModel> const& camera_list) :
    m_lonlat(georef.pixel_to_lonlat(Vector2(x, y))), 
    m_georef(get_crop_georef(georef, BBox2i(x - half_kern, y - half_kern, whole_kern, whole_kern))), 
    m_image_list(image_list), m_camera_list(camera_list),
    m_num_patches(image_list.size()) {}
     
  result_type operator()(domain_type const& x) const {
    return result_type(0);
  } 

  gradient_type gradient(domain_type const& x) const {
    return gradient_type();
  }

  std::vector<ImageView<float32> > get_ortho_patches(domain_type const& x) const {
    // convert to gen_synth_scene plane type
    double radius = x[3] + m_georef.datum().radius(m_lonlat[0], m_lonlat[1]);
    Vector3 pt = cartography::lon_lat_radius_to_xyz(Vector3(m_lonlat[0], m_lonlat[1], radius));
    Vector3 normal = math::SubVector<Vector4>(x, 0, 3);
    Vector4 plane(normal[0], normal[1], normal[2], dot_prod(normal, pt));
   
    // project onto patch_list 
    std::vector<ImageView<float32> > patch_list(m_num_patches);

    for (unsigned i = 0; i < m_num_patches; i++) {
      patch_list[i] = cartography::orthoproject(plane_dem_view(m_georef, plane, 
                                                               whole_kern, whole_kern),
                                                m_georef, m_image_list[i], &m_camera_list[i],
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
