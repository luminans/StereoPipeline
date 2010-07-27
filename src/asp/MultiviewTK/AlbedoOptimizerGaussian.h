#ifndef __ASP_MULTIVIEWTK_ALBEDOOPTIMIZERGAUSSIAN_H__
#define __ASP_MULTIVIEWTK_ALBEDOOPTIMIZERGAUSSIAN_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Manipulation.h>

#include <boost/foreach.hpp>

#include <vw/FileIO.h>

namespace vw {
namespace multiview { 

template <class PatchT, class AlbedoT>
std::vector<ImageView<float32> >
objective_gaussian(std::vector<PatchT> const& patch,
                   AlbedoT const& albedo,
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
  static const unsigned MAX_ITER = 20;
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

    // Find a good step size
    float32 step_size = 1;
    ImageView<float32> tmp_albedo;
    std::vector<float32> tmp_b(num_patches), tmp_c(num_patches);
    float32 tmp_obj_sum;
    for (unsigned j = 0; j < MAX_ITER; j++) {
      tmp_albedo = albedo - albedo_grad * step_size;
      for (unsigned k = 0; k < num_patches; k++) {
        tmp_b[k] = b[k] - b_grad[k] * step_size;
        tmp_c[k] = c[k] - c_grad[k] * step_size;
      }

      // Normalize so b_0 = 0, c_0 = 1
      tmp_albedo = tmp_b[0] + tmp_c[0] * copy(tmp_albedo);
      for (unsigned k = 0; k < num_patches; k++) {
        tmp_b[k] = tmp_b[k] - tmp_c[k] / tmp_c[0] * tmp_b[0];
        tmp_c[k] = tmp_c[k] / tmp_c[0];
      } 

      // Recalculate obj;
      tmp_obj_sum = 0;
      BOOST_FOREACH(ImageView<float32> obj, objective_gaussian(patch, tmp_albedo, tmp_b, tmp_c)) {
        tmp_obj_sum += sum_of_pixel_values(obj);
      }

      if (tmp_obj_sum < old_obj_sum) {
        break;
      } else {
        step_size /= 10;
      }
    }

    std::cout << "step size: " << step_size << std::endl;
    std::cout << "b: ";
    BOOST_FOREACH(float32 b, tmp_b) {
      std::cout << " " << b;
    }
    std::cout << std::endl;
    std::cout << "c: ";
    BOOST_FOREACH(float32 c, tmp_c) {
      std::cout << " " << c;
    }
    std::cout << std::endl;

    // Apply grads
    albedo = tmp_albedo;
    std::stringstream s;
    s << "albedo/albedo-" << i << ".tif";
    write_image(s.str(), channel_cast_rescale<uint8>(normalize(albedo)));
    b = tmp_b;
    c = tmp_c;
    obj_sum = tmp_obj_sum;

    std::cout << "prev obj: " << old_obj_sum << " curr obj: " << obj_sum << std::endl;
    std::cout << "improvement: " << (old_obj_sum - obj_sum) / obj_sum << std::endl;
  }

  return albedo;
}



}} //namespace vw, multiview

#endif
