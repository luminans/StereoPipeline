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

class AlbedoOptimizerTest : public ::testing::Test {
protected:
  AlbedoOptimizerTest() {
    unsigned int num_patches = 4;
    float32 b_ground_data[] = { 3, 4, 5, 6 };
    float32 c_ground_data[] = { 4, 8, 1, 2 };

    b_ground = Vector<float32>(num_patches, b_ground_data);
    c_ground = Vector<float32>(num_patches, c_ground_data);

    boost::rand48 gen(10);
    albedo_ground = uniform_noise_view(gen, 30, 30);

    for (unsigned i = 0; i < num_patches; i++) {
      patch_list.push_back(b_ground[i] + c_ground[i] * albedo_ground);
    }

    // Normalize everything to b_0 = 0, c_0 = 1
    albedo_ground = b_ground[0] + c_ground[0] * copy(albedo_ground);
    b_ground = b_ground - b_ground[0] * c_ground / c_ground[0];
    c_ground = c_ground / c_ground[0];

    //vw_log().console_log().rule_set().add_rule(40, "math");
  }

  std::vector<ImageView<float32> > patch_list;
  ImageView<float32> albedo_ground;
  Vector<float32> b_ground;
  Vector<float32> c_ground;
};

TEST_F(AlbedoOptimizerTest, AlbedoOptimizerGaussian) {
  AlbedoOptimizerData result;
  result = math::conjugate_gradient(AlbedoOptimizerGaussian(patch_list),
                                    AlbedoOptimizerData(patch_list),
                                    math::ArmijoStepSize(0.01),
                                    100, 1e-8);

  // Normalize to b_0 = 0, c_0 = 1
  result.albedo = result.b[0] + result.c[0] * copy(result.albedo);
  result.b = result.b - result.b[0] * result.c / result.c[0];
  result.c = result.c / result.c[0];

  float32 avg_albedo_err = sum_of_pixel_values(abs(albedo_ground - result.albedo)) /
                           albedo_ground.cols() / albedo_ground.rows();

  EXPECT_LT(avg_albedo_err, 0.001);
  EXPECT_VECTOR_NEAR(b_ground, result.b, 0.01);
  EXPECT_VECTOR_NEAR(c_ground, result.c, 0.01);
}
