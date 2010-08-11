#ifndef __ASP_MULTIVIEWTK_MULTIVIEW_H__
#define __ASP_MULTIVIEWTK_MULTIVIEW_H__

#include <asp/MultiviewTK/PatchMultiView.h>
#include <asp/MultiviewTK/gen_synth_scene.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include <vw/FileIO.h>
#include <vw/Image.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <vw/Camera.h>


struct Options {
  vw::BBox2i bbox;
  vw::float32 nodata_value;
  std::string dem_name;
  std::vector<std::string> image_names;
  std::vector<std::string> camera_names;
  std::string output_prefix;
};

std::ostream& operator<<(std::ostream& os, const Options& opts) {
  os << "bbox: " << opts.bbox << std::endl;
  os << "dem_name: " << opts.dem_name << std::endl;
  os << "nodata_value: " << opts.nodata_value << std::endl;

  os << "image_names:";
  BOOST_FOREACH(std::string name, opts.image_names) {
    os << " " << name;
  }
  os << std::endl;

  os << "camera_names:";
  BOOST_FOREACH(std::string name, opts.camera_names) {
    os << " " << name;
  }
  os << std::endl;

  os << "output_prefix: " << opts.output_prefix << std::endl;

  return os;
}

Options parse_opts(int argc, char *argv[]) {
  namespace po = boost::program_options;
  namespace fs = boost::filesystem;

  Options opts;
  std::vector<int> bbox_bounds;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Display this help message")
    ("bbox", po::value<std::vector<int> >(&bbox_bounds)->multitoken(), "xoffset yoffset width height")
    ("nodata-value", po::value<vw::float32>(&opts.nodata_value)->default_value(0), "nodata value for DEM")
    ("output-prefix,o", po::value<std::string>(&opts.output_prefix), "Explicitly specify the output prefix")
    ("dem", po::value<std::string>(&opts.dem_name), "Explicitly specify the dem to refine")
    ("image", po::value<std::vector<std::string> >(&opts.image_names)->multitoken(), "Explicitly specify the images")
    ("camera", po::value<std::vector<std::string> >(&opts.camera_names)->multitoken(), "Specify the camera files (defaults to [image_name].pinhole")
    ;

  po::positional_options_description p;
  p.add("output-prefix", 1);
  p.add("dem", 1);
  p.add("image", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << "Usage: " << argv[0] << " output-prefix dem image1 image2 [image3 ...]" << std::endl;
    std::cout << desc << std::endl;
    exit(1);
  }

  if (!vm.count("dem") || !vm.count("output-prefix") ||
      opts.image_names.size() < 2) {
    std::cout << "Usage: " << argv[0] << " output-prefix dem image1 image2 [image3 ...]" << std::endl;
    std::cout << desc << std::endl;
    exit(1);
  }

  if (vm.count("bbox")) {
    if (bbox_bounds.size() == 4) {
      opts.bbox = vw::BBox2i(bbox_bounds[0], bbox_bounds[1], bbox_bounds[2], bbox_bounds[3]);
    } else {
      std::cout << "Error: --bbox must specify 4 arguments: xoffset yoffset width height" << std::endl;
      exit(1);
    }
  } else {
    vw::DiskImageView<vw::float32> dem(opts.dem_name);
    opts.bbox = bounding_box(dem);
  }

  if (!vm.count("camera")) {
    assert(opts.camera_names.empty());
    BOOST_FOREACH(std::string image, opts.image_names) {
      opts.camera_names.push_back(fs::path(image).replace_extension(".pinhole").string());
    }
  } else {
    if (opts.image_names.size() != opts.camera_names.size()) {
      std::cout << "Error: Please specify exactly one camera model per image" << std::endl;
      exit(1);
    }
  }

  return opts;
}

void
translate_georef(vw::cartography::GeoReference& georef, vw::Vector2i const& offset) {
  vw::Matrix3x3 affine = georef.transform();
  vw::Vector2 lonlat_offset = georef.pixel_to_lonlat(offset);
  affine(0, 2) = lonlat_offset.x();
  affine(1, 2) = lonlat_offset.y();
  georef.set_transform(affine);
}

template <class DemT>
vw::Vector4 approx_dem_plane(vw::int32 x, vw::int32 y, DemT dem, vw::cartography::GeoReference georef) {
  using namespace vw;
  std::vector<Vector2i> img_pts(3);
  std::vector<Vector3> pts(3);
  img_pts[0] = Vector2(x, y);
  img_pts[1] = Vector2(x - 1, y - 1);
  img_pts[2] = Vector2(x - 1, y);

  for (unsigned i = 0; i < 3; i++) {
    Vector2 lonlat = georef.pixel_to_lonlat(img_pts[i]);
    double rad = georef.datum().radius(lonlat[0], lonlat[1]) + 
                 dem(img_pts[i][0], img_pts[i][1]);
    pts[i] = cartography::lon_lat_radius_to_xyz(Vector3(lonlat[0], lonlat[1], rad));
  }

  return gen_plane(pts[0], pts[1], pts[2]);
}

#endif
