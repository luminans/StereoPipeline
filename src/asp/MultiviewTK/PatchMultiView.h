#ifndef __ASP_MULTIVIEWTK_REFMULTIVIEW_H__
#define __ASP_MULTIVIEWTK_REFMULTIVIEW_H__

#include <boost/foreach.hpp>

#include <vw/Cartography/GeoReference.h>
#include <vw/Camera/PinholeModel.h>

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Statistics.h>

namespace vw {
namespace multiview {

template <class DemT, class ImageT>
class PatchMultiView : public ImageViewBase<PatchMultiView<DemT, ImageT> > {
  DemT m_dem;
  cartography::GeoReference m_georef;
  std::vector<ImageT> m_image_list;
  std::vector<camera::PinholeModel> m_camera_list;
  public:
    typedef PixelMask<float32> pixel_type;
    typedef PixelMask<float32> const& result_type; // Should this be const or const&?
    typedef ProceduralPixelAccessor<PatchMultiView> pixel_accessor;

    PatchMultiView(DemT const& dem,
                   cartography::GeoReference georef,
                   std::vector<ImageT> const& image_list,
                   std::vector<camera::PinholeModel> const& camera_list) :
      m_dem(dem), m_georef(georef), m_image_list(image_list), m_camera_list(camera_list) {}

    inline int32 cols() const { return m_dem.cols(); }
    inline int32 rows() const { return m_dem.rows(); }
    inline int32 planes() const { return 1; }

    inline pixel_accessor origin() const { return pixel_accessor(*this); }

    inline result_type operator()(int32 col, int32 row, int32 /*plane*/ = 0) const {
      if (!is_valid(m_dem(col, row)))
        return result_type();

      // Deterine initial surface pi by looking at a window around m_dem
      // Initial bi, ci = 1
      //
      // Loop until converge:
      //   Create planar DEM patch using surface pi
      //   Project all images onto this patch
      //   Optimize bi, ci, a
      //   Calculate gradient for surface pi
      //   Search along gradient for min, update surface pi

      return result_type(); 
    }

    typedef PatchMultiView<typename DemT::prerasterize_type, ImageT> prerasterize_type;
    inline prerasterize_type prerasterize(BBox2i bbox) const {
      // TODO: compute the bounding box on each ImageT that we need to prerasterize
      // Right now we just hope they're lightweight
      return prerasterize_type(m_dem.prerasterize(bbox), m_georef, m_image_list, m_camera_list); 
    }

    template <class DestT>
    inline void rasterize(DestT const& dest, BBox2i bbox) const {
      vw::rasterize(prerasterize(bbox), dest, bbox);
    }  
};

template <class DemT, class ImageT>
PatchMultiView<DemT, ImageT> 
patch_multiview(ImageViewBase<DemT> const& dem,
                cartography::GeoReference const& georef,
                std::vector<ImageT> const& image_list,
                std::vector<camera::PinholeModel> const& camera_list) {
  return PatchMultiView<DemT, ImageT>(dem.impl(), georef, image_list, camera_list);
}

}} //namesace vw, multiview

#endif
