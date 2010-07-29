#ifndef __ASP_MULTIVIEWTK_GENSYNTHSCENE_H__
#define __ASP_MULTIVIEWTK_GENSYNTHSCENE_H__

#include <boost/foreach.hpp>

#include <vw/FileIO.h>
#include <vw/Image.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <vw/Camera.h>

vw::cartography::GeoReference gen_dem_georef() {
  using namespace vw;
  using namespace vw::cartography;
  // 512x1024 DEM georef
  // Original DEM mean: -2604.93 stddev:888.459
  double georef_affine_data[] = { 0.00133653359715, 0.00000000000000, 56.1392457178,
                                  0.00000000000000,-0.00133653359715, 10.1949968063,
                                  0.00000000000000, 0.00000000000000, 1.00000000000 };
  return GeoReference(Datum("D_MOON"), MatrixProxy<double>(georef_affine_data, 3, 3));
}

std::vector<vw::camera::PinholeModel> gen_camera_list() {
  using namespace vw;
  using namespace vw::camera;
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

  return camera_list;
}

#endif
