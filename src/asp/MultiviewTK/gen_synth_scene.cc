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
  GeoReference georef = gen_dem_georef();
  std::vector<PinholeModel> camera_list = gen_camera_list();

  return 0;
}
