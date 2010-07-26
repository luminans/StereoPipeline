#ifndef __ASP_MULTIVIEWTK_ALBEDOOPTIMIZERGAUSSIAN_H__
#define __ASP_MULTIVIEWTK_ALBEDOOPTIMIZERGAUSSIAN_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Algorithms.h>

#include <boost/foreach.hpp>

namespace vw {
namespace multiview { 

std::vector<ImageView<float32> >
objective_gaussian(std::vector<ImageView<float32> > const& patch,
                   ImageView<float32> const& albedo,
                   std::vector<float32> const& b,
                   std::vector<float32> const& c) {
  std::vector<ImageView<float32> > obj(patch.size());
  for (unsigned i = 0; i < patch.size(); i++) {
    obj[i] = (b[i] + c[i] * albedo - patch[i]) *
             (b[i] + c[i] * albedo - patch[i]);
  }
  return obj;
}

template <class PatchT>
ImageView<float32> find_albedo_gaussian(std::vector<PatchT> const& patch) {
  static const unsigned MAX_ITER = 10;
  static const float32 CONV_TOL = 1e-6;

  unsigned num_patches = patch.size();

  // Initialize albedo = 1, b_k = 0, c_k = 1
  ImageView<float32> albedo(patch[0].cols(), patch[0].rows());
  fill(albedo, 1);
  std::vector<float32> b(num_patches, 0);
  std::vector<float32> c(num_patches, 1);

  ImageView<float32> albedo_grad(patch[0].cols(), patch[0].rows());
  std::vector<float32> b_grad(num_patches);
  std::vector<float32> c_grad(num_patches);

  float32 old_obj_sum;

  float32 obj_sum = 0;
  BOOST_FOREACH(ImageView<float32> obj, objective_gaussian(patch, albedo, b, c)) {
    obj_sum += sum_of_pixel_values(obj);
  }

  for (unsigned i = 0; i < MAX_ITER; i++) {
    old_obj_sum = obj_sum;

    // Calculate grads
    fill(albedo_grad, 0);
    for (unsigned k = 0; k < num_patches; k++) {
      ImageView<float32> temp = b[k] + c[k] * albedo - patch[k];
      b_grad[k] = sum_of_pixel_values(temp);
      c_grad[k] = sum_of_pixel_values(albedo * temp);
      albedo_grad += c[k] * temp;
    }

    // Apply grads
    for (unsigned k = 0; k < num_patches; k++) {
      b[k] += -b_grad[k] * 0.001;
      c[k] += -c_grad[k] * 0.001;
      albedo += -albedo_grad * 0.001;
    }

    // Recalculate obj;
    obj_sum = 0;
    BOOST_FOREACH(ImageView<float32> obj, objective_gaussian(patch, albedo, b, c)) {
      obj_sum += sum_of_pixel_values(obj);
    }

    std::cout << "prev obj: " << old_obj_sum << " curr obj: " << obj_sum << std::endl;
    std::cout << "improvement: " << (old_obj_sum - obj_sum) / obj_sum << std::endl;
  }

  return albedo;
}



}} //namespace vw, multiview

#endif
