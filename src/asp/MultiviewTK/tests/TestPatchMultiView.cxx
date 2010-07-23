// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <gtest/gtest.h>
#include <test/Helpers.h>

#include <asp/MultiviewTK/PatchMultiView.h>

#include <vw/Image.h>

using namespace vw;
using namespace vw::multiview;

class PatchMultiViewTest : public ::testing::Test {
protected:
  PatchMultiViewTest() {}

  virtual void SetUp() {
  }

  // Variables here
};

TEST_F(PatchMultiViewTest, DoSomething) {
  EXPECT_EQ(0, 0);
}
