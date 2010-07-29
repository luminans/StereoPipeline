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
  // 512x1024 DEM georef
  double georef_affine_data[] = { 0.00133653359715, 0.00000000000000, 56.1392457178,
                                  0.00000000000000,-0.00133653359715, 10.1949968063,
                                  0.00000000000000, 0.00000000000000, 1.00000000000 };
  GeoReference georef(Datum("D_MOON"), MatrixProxy<double>(georef_affine_data, 3, 3));

  Options opts = parse_opts(argc, argv);

  std::vector<PinholeModel> camera_list_cmp;

  BOOST_FOREACH(string camera, opts.camera_names) {
    camera_list_cmp.push_back(PinholeModel(camera));
  }

  cout << std::setprecision(12) << endl;

  std::vector<PinholeModel> camera_list;

  camera_list.push_back(PinholeModel(
    Vector3(966089.223462, 1557938.52831, 282405.060851),
    Quat(-0.0470851319085,0.358002222657,0.665010992829,-0.653741369612).rotation_matrix(),
    3802.7, 3802.7, 2757.875, 814.875,
    Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1), NullLensDistortion()
    ));

  camera_list.push_back(PinholeModel(
    Vector3(996045.240112, 1535573.08189, 299032.543447),
    Quat(-0.0540368023579, 0.364928919676, 0.665264280335, -0.649099641723).rotation_matrix(),
    3802.7, 3802.7, 1543.875, 823.875,
    Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1), NullLensDistortion()
    ));

  camera_list.push_back(PinholeModel(
    Vector3(1025532.88151, 1512446.05054, 315515.076309),
    Quat(-0.0588100780917, 0.371592737959, 0.667265027691, -0.642820032849).rotation_matrix(),
    3802.7, 3802.7, 363.875, 863.875,
    Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1), NullLensDistortion()
    ));

  camera_list.push_back(PinholeModel(
    Vector3(1054501.51827, 1488593.63697, 331834.757985),
    Quat(-0.0639146875174, 0.378441350432, 0.667943479689, -0.637611609793).rotation_matrix(),
    3802.7, 3802.7, -836.125, 863.875,
    Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1), NullLensDistortion()
    ));

  for (int i = 0; i < 4; i++) {
    cout << camera_list_cmp[i].point_to_pixel(camera_list[i].camera_center() + camera_list[i].pixel_to_vector(Vector2(100, 200))) << endl;
  }

  return 0;
}
