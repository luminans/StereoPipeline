// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <asp/MultiviewTK/multiview.h>
#include <asp/MultiviewTK/GeometryOptimizer.h>
#include <vw/Math/ConjugateGradient.h>

using std::cout;
using std::endl;
using std::string;

namespace fs = boost::filesystem;
using namespace vw;
using namespace vw::camera;
using namespace vw::cartography;
using namespace vw::multiview;

int main( int argc, char *argv[] ) {
  Options opts = parse_opts(argc, argv);

  std::vector<DiskImageView<float32> > image_list;
  std::vector<PinholeModel> camera_list;

  BOOST_FOREACH(string image, opts.image_names) {
    image_list.push_back(DiskImageView<float32>(image));
  }

  BOOST_FOREACH(string camera, opts.camera_names) {
    camera_list.push_back(PinholeModel(camera));
  }

  ImageViewRef<PixelMask<float32> > dem = 
    create_mask(crop(DiskImageView<float32>(opts.dem_name), opts.bbox), 
                opts.nodata_value);

  // TODO: All DEMs should really already have a transparency 
  // channel so nodata-values aren't needed... Fix point2dem!
  if (!opts.nodata_value) {
    dem = validate_mask(dem);
  }

  GeoReference georef;
  read_georeference(georef, opts.dem_name);
  translate_georef(georef, opts.bbox.min());

  int col = 500, row = 500;
  GeometryOptimizer<DiskImageView<float32> > go(col, row, georef, image_list,
                                                            camera_list);

  Vector4 plane = approx_dem_plane(col, row, dem, georef);

  vw_log().console_log().rule_set().add_rule(40, "math");

  Vector4 result = math::conjugate_gradient(go, plane, math::ArmijoStepSize(.1), 500);

  return 0;
}
