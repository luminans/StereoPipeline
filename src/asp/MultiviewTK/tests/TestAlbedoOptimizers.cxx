// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <boost/random/linear_congruential.hpp>

#include <asp/MultiviewTK/AlbedoOptimizers.h>

#include <vw/Image.h>
#include <vw/FileIO.h>

using namespace vw;
using namespace vw::multiview;

TEST(AlbedoOptimizers, AlbedoOptimizerGaussian) {
  static const int NUM_PATCHES = 4;
  float32 b_ground[] = { 3, 4, 5, 6 };
  float32 c_ground[] = { 4, 8, 1, 2 };

  boost::rand48 gen(10);
  ImageView<float32> albedo_ground = uniform_noise_view(gen, 30, 30);

  std::vector<ImageView<float32> > patch_list(NUM_PATCHES);
  for (int i = 0; i < NUM_PATCHES; i++) {
    patch_list[i] = b_ground[i] + c_ground[i] * albedo_ground;
  }

  vw_log().console_log().rule_set().add_rule(40, "math");
  AlbedoOptimizerData result;
  result = math::conjugate_gradient(AlbedoOptimizerGaussian(patch_list),
                                    AlbedoOptimizerData(patch_list),
                                    math::ArmijoStepSize(0.01),
                                    100, 1e-8);

  std::cout << "c: " << result.c / result.c[0] << std::endl;
  std::cout << "b: " << result.b - result.b[0] * result.c / result.c[0] << std::endl;

  write_image("albedo_ground.tif", channel_cast_rescale<uint8>(normalize(albedo_ground)));
  write_image("albedo_est.tif", channel_cast_rescale<uint8>(normalize(result.albedo)));
}
