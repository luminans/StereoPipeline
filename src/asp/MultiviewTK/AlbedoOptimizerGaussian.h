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

class AlbedoOptimizerGaussian {
  typedef ImageView<float32> PatchT;
  typedef ImageView<float32> ObjT;

  std::vector<PatchT> m_patch;
  std::vector<ObjT> m_obj;
  std::vector<float32> m_b;
  std::vector<float32> m_c;
  unsigned m_num_patches;

  PatchT m_albedo;

  float32 m_obj_sum;

  void calc_objective() {
    m_obj_sum = 0;
    for (unsigned i = 0; i < m_num_patches; i++) {
      m_obj[i] = (m_b[i] + m_c[i] * m_albedo - m_patch[i]) *
                 (m_b[i] + m_c[i] * m_albedo - m_patch[i]);
      m_obj_sum += sum_of_pixel_values(m_obj[i]);
    }
  }

  public:
    AlbedoOptimizerGaussian(std::vector<PatchT> const& patch) :
      m_patch(patch), m_obj(patch.size()), m_b(patch.size(), 0), m_c(patch.size(), 1),
      m_num_patches(patch.size()) {
      // VW automatically fills with 0's
      PatchT albedo_num(m_patch[0].cols(), m_patch[0].rows());
      float32 albedo_denom = 0;
      for (unsigned i = 0; i < m_num_patches; i++) {
        albedo_num += m_c[i] * (m_patch[i] - m_b[i]);
        albedo_denom += m_c[i] * m_c[i];
      }

      m_albedo = albedo_num / albedo_denom;

      calc_objective();
    }

    PatchT optimize_albedo() {
      static const unsigned MAX_ITER = 10;
      static const float32 CONV_TOL = 1e-6;
 
      PatchT albedo_grad(m_patch[0].cols(), m_patch[0].rows()); 
      std::vector<float32> b_grad(m_num_patches);
      std::vector<float32> c_grad(m_num_patches);

      float32 old_obj_sum;

      for (unsigned i = 0; i < MAX_ITER; i++) {
        old_obj_sum = m_obj_sum;

        // Calculate grads
        fill(albedo_grad, 0);
        for (unsigned k = 0; k < m_num_patches; k++) {
          PatchT temp = m_b[k] + m_c[k] * m_albedo - m_patch[k];
          b_grad[k] = sum_of_pixel_values(temp);
          c_grad[k] = sum_of_pixel_values(m_albedo * temp);
          albedo_grad += m_c[k] * temp;
        }

        // Apply grads
        for (unsigned k = 0; k < m_num_patches; k++) {
          m_b[k] += -b_grad[k] * 0.0001;
          m_c[k] += -c_grad[k] * 0.0001;
          m_albedo += -albedo_grad * 0.0001;
        }

        // Recalculate obj;
        calc_objective();

        std::cout << "prev obj: " << old_obj_sum << " curr obj: " << m_obj_sum << std::endl;
      }

      return m_albedo;
    }
};  

}} //namespace vw, multiview

#endif
