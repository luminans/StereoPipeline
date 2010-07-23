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

  ImageViewRef<float32> dem = crop(DiskImageView<float32>(opts.dem_name),
                                   opts.bbox);
  GeoReference georef = get_crop_georef(opts.dem_name, opts.bbox);

  for (unsigned i = 0; i < image_list.size(); i++) {
    std::stringstream output_name;
    output_name << opts.output_prefix << "-" << i << ".tif";

    boost::shared_ptr<CameraModel> cam_ptr(new PinholeModel(camera_list[i]));

    ImageViewRef<PixelMask<float32> > result_fp = orthoproject(dem, georef, 
                                                   image_list[i], cam_ptr,
                                                   BilinearInterpolation(), ZeroEdgeExtension());
    ImageViewRef<PixelMask<PixelGray<uint8> > > result = 
      pixel_cast_rescale<PixelMask<PixelGray<uint8> > >(normalize(result_fp));

    {    
      DiskImageResourceGDAL rsrc(output_name.str(), result.format(), Vector2i(256, 256));
      write_georeference(rsrc, georef);
      block_write_image(rsrc, result, TerminalProgressCallback("vw", output_name.str() + ": "));
    }
  }

  return 0;
}
