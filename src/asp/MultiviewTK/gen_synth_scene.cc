// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <asp/MultiviewTK/gen_synth_scene.h>

using std::cout;
using std::endl;
using std::string;

using namespace vw;
using namespace vw::camera;
using namespace vw::cartography;

int main( int argc, char *argv[] ) {
  const int dem_width = 1024, dem_height = 1024;
  //TODO: boost args
  string output_folder = ".";
  GeoReference georef = gen_dem_georef();
  std::vector<PinholeModel> camera_list = gen_camera_list();

  // Create DEM-initial.tif
  Vector4 dem_initial_plane = gen_plane(georef, -2605, dem_width, dem_height);
  write_image(output_folder + "/DEM-initial.tif", 
              pixel_cast<float32>(plane_dem_view(georef,
              dem_initial_plane, dem_width, dem_height)));

  // Create DEM-ground.tif
  Vector4 dem_ground_plane = gen_plane(georef, -2605, 888, dem_width, dem_height);
  write_image(output_folder + "/DEM-ground.tif", 
              pixel_cast<float32>(plane_dem_view(georef,
              dem_ground_plane, dem_width, dem_height)));

  // Create DRG-ground.tif
  boost::rand48 gen;
  write_image(output_folder + "/DRG-ground.tif",
              pixel_cast<float32>(uniform_noise_view(gen, dem_width, dem_height)));
  DiskImageView<float32> drg_ground(output_folder + "/DRG-ground.tif");

  /*
  for (unsigned i = 0; i < camera_list.size(); i++) {
    stringstream ss;
    ss << output_folder << "/" << i << ".tif";
    write_image(ss.str(), 
                gen_orbital_image(drg_ground, 
                                  camera_list[i], 
                                  dem_ground_plane
                                  1800, 1800));

  } */

  return 0;
}
