#ifndef __ASP_MULTIVIEWTK_ALBEDOOPTIMIZERS_H__
#define __ASP_MULTIVIEWTK_ALBEDOOPTIMIZERS_H__

#include <vw/Image/ImageViewBase.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/ImageMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/Statistics.h>
#include <vw/Image/Algorithms.h>
#include <vw/Image/Manipulation.h>

#include <boost/foreach.hpp>

#include <vw/FileIO.h>

#include <vw/Math/ConjugateGradient.h>

namespace vw {
namespace multiview { 

struct AlbedoOptimizerData {
  ImageView<float32> albedo;
  Vector<float32> b;
  Vector<float32> c;

  AlbedoOptimizerData() {}

  template <class PatchT>
  AlbedoOptimizerData(std::vector<PatchT> patch_list) :
    albedo(patch_list[0].cols(), patch_list[0].rows()), 
    b(patch_list.size()), c(patch_list.size()) 
  {
    fill(albedo, 1);
    std::fill(b.begin(), b.end(), 0);
    std::fill(c.begin(), c.end(), 1);    
  }

  template <class AlbedoT>
  AlbedoOptimizerData(AlbedoT ca, Vector<float32> cb, Vector<float32> cc) :
    albedo(ca), b(cb), c(cc) {}

  AlbedoOptimizerData operator+ (AlbedoOptimizerData const& x) const {
    return AlbedoOptimizerData(albedo + x.albedo, b + x.b, c + x.c);
  }
};

template <class ScalarT>
typename boost::enable_if<IsScalar<ScalarT>, AlbedoOptimizerData>::type
operator*(ScalarT s, AlbedoOptimizerData const& x) {
  return AlbedoOptimizerData(s * x.albedo, s * x.b, s * x.c);
}

float32 dot_prod(AlbedoOptimizerData const& x, AlbedoOptimizerData const& y) {
  float32 result = sum_of_pixel_values(x.albedo * y.albedo);
  result += dot_prod(x.b, y.b);
  result += dot_prod(x.c, y.c);
  return result;
}

struct AlbedoOptimizerGaussian {
  typedef float32 result_type;
  typedef AlbedoOptimizerData domain_type;
  typedef AlbedoOptimizerData gradient_type;

  std::vector<ImageView<float32> > m_patch_list;
  const unsigned m_num_patches;

  template <class PatchT>
  AlbedoOptimizerGaussian(std::vector<PatchT> const& patch_list) :
    m_patch_list(patch_list.size()), m_num_patches(patch_list.size()) 
  {
    for (unsigned i = 0; i < m_num_patches; i++) {
      m_patch_list[i] = patch_list[i];
    } 
  } 

  result_type operator()(domain_type const& x) const {
    result_type result = 0;
    for (unsigned i = 0; i < m_num_patches; i++) {
      result += sum_of_pixel_values((x.b[i] + x.c[i] * x.albedo - m_patch_list[i]) *
                                    (x.b[i] + x.c[i] * x.albedo - m_patch_list[i]));
    }
    return result;
  }

  gradient_type gradient(domain_type const& x) const {
    ImageView<float32> albedo_grad(m_patch_list[0].cols(), m_patch_list[0].rows());
    Vector<float32> b_grad(m_num_patches);
    Vector<float32> c_grad(m_num_patches);
    fill(albedo_grad, 0);
    for (unsigned i = 0; i < m_num_patches; i++) {
      ImageView<float32> temp = x.b[i] + x.c[i] * x.albedo - m_patch_list[i];
      albedo_grad += x.c[i] * temp;
      b_grad[i] = sum_of_pixel_values(temp);
      c_grad[i] = sum_of_pixel_values(x.albedo * temp);
    }
    return gradient_type(albedo_grad, b_grad, c_grad);
  }

  // It looks like dimension can be any non-zero value...
  unsigned dimension() const { 
    return m_num_patches * m_patch_list[0].cols() * m_patch_list[0].rows(); 
  }
};

}} //namespace vw, multiview

#endif
