// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <boost/random/linear_congruential.hpp>

#include <asp/MultiviewTK/AlbedoOptimizerGaussian.h>

#include <vw/Image.h>
#include <vw/FileIO.h>

using namespace vw;
using namespace vw::multiview;

TEST(AlbedoOptimizerGaussian, TestOptimize) {
  typedef ImageView<float32> PatchT;

  static const int NUM_PATCHES = 4;
  float32 b_ground[] = { 3, 4, 5, 6 };
  float32 c_ground[] = { 4, 4.4, 4.2, 4.8 };

  boost::rand48 gen(10);
  PatchT albedo_ground = gaussian_noise_view(gen, 0.5, 0.25, 10, 10);

  std::vector<PatchT> patch(NUM_PATCHES);
  for (int i = 0; i < NUM_PATCHES; i++) {
    patch[i] = b_ground[i] + c_ground[i] * albedo_ground;
  }

  AlbedoOptimizerGaussian optimizer(patch);
  PatchT albedo_est = optimizer.optimize_albedo();

  write_image("albedo_ground.tif", channel_cast_rescale<uint8>(normalize(albedo_ground)));
  write_image("albedo_est.tif", channel_cast_rescale<uint8>(normalize(albedo_est)));
}
