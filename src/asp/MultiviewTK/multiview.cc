// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include "multiview.h"

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

  std::vector<DiskImageView<PixelMask<float32> > > image_list;
  std::vector<PinholeModel> camera_list;

  BOOST_FOREACH(string image, opts.image_names) {
    image_list.push_back(DiskImageView<PixelMask<float32> >(image));
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

  GeoReference georef = get_crop_georef(opts.dem_name, opts.bbox);

  ImageViewRef<PixelMask<float32> > refined_dem =
    patch_multiview(dem, georef, image_list, camera_list);

  return 0;
}
