/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/synthesis.hpp"

namespace hmap
{

namespace
{

// --- helpers

// Radial polynomial mask of the reference implementation: 1 at the patch
// center, decreasing to 0 before the corners, ((1 - (1 - 1/size) r^2)^+)^2.
Array build_patch_mask(int size)
{
  Array mask(glm::ivec2(size, size));

  float radius = 0.5f * (float)(size - 1);
  float offset = 1.f - 1.f / (float)size;

  for (int j = 0; j < size; j++)
    for (int i = 0; i < size; i++)
    {
      float x = ((float)i - radius) / radius;
      float y = ((float)j - radius) / radius;
      float v = std::max(0.f, 1.f - offset * (x * x + y * y));
      mask(i, j) = v * v;
    }

  return mask;
}

// Block-mean pooling, so that the mean of a low-resolution patch equals the
// mean of the paired high-resolution patch.
Array downsample_block_mean(const Array &array, int factor)
{
  if (factor == 1) return array;

  Array out(glm::ivec2(array.shape.x / factor, array.shape.y / factor));

  for (int j = 0; j < out.shape.y; j++)
    for (int i = 0; i < out.shape.x; i++)
    {
      float sum = 0.f;
      for (int q = 0; q < factor; q++)
        for (int p = 0; p < factor; p++)
          sum += array(i * factor + p, j * factor + q);
      out(i, j) = sum / (float)(factor * factor);
    }

  return out;
}

// Clamp-to-edge padding by `nbuffer` cells on every side.
Array pad_replicate(const Array &array, int nbuffer)
{
  Array out(
      glm::ivec2(array.shape.x + 2 * nbuffer, array.shape.y + 2 * nbuffer));

  for (int j = 0; j < out.shape.y; j++)
  {
    int jc = std::clamp(j - nbuffer, 0, array.shape.y - 1);
    for (int i = 0; i < out.shape.x; i++)
    {
      int ic = std::clamp(i - nbuffer, 0, array.shape.x - 1);
      out(i, j) = array(ic, jc);
    }
  }

  return out;
}

// --- dictionary

struct PatchDictionary
{
  int                dim_lr = 0;
  int                dim_hr = 0;
  int                n_atoms = 0;
  std::vector<float> atoms_lr; // n_atoms x dim_lr, unit L2 norm
  std::vector<float> atoms_hr; // n_atoms x dim_hr, scaled by the LR norm
};

PatchDictionary build_dictionary(const Array &exemplar_hr,
                                 const Array &exemplar_lr,
                                 int          factor,
                                 int          patch_size,
                                 int          stride,
                                 const Array &mask_lr,
                                 const Array &mask_hr)
{
  PatchDictionary dict;

  const int n = patch_size;
  const int nh = patch_size * factor;

  dict.dim_lr = n * n;
  dict.dim_hr = nh * nh;

  const int d1 = (exemplar_lr.shape.x - n) / stride + 1;
  const int d2 = (exemplar_lr.shape.y - n) / stride + 1;

  // patches with no variation (relative to the exemplar range) are useless
  // as atoms and would blow up the normalization
  const float norm_min = 1e-6f * (exemplar_lr.max() - exemplar_lr.min());

  dict.atoms_lr.reserve((size_t)d1 * d2 * dict.dim_lr);
  dict.atoms_hr.reserve((size_t)d1 * d2 * dict.dim_hr);

  std::vector<float> v_lr(dict.dim_lr);
  std::vector<float> v_hr(dict.dim_hr);

  for (int jp = 0; jp < d2; jp++)
    for (int ip = 0; ip < d1; ip++)
    {
      const int i0 = ip * stride;
      const int j0 = jp * stride;

      // low-resolution atom: mean removed, masked, unit norm
      double mean = 0.0;
      for (int q = 0; q < n; q++)
        for (int p = 0; p < n; p++)
          mean += exemplar_lr(i0 + p, j0 + q);
      mean /= (double)dict.dim_lr;

      double norm2 = 0.0;
      for (int q = 0; q < n; q++)
        for (int p = 0; p < n; p++)
        {
          float v = (float)(exemplar_lr(i0 + p, j0 + q) - mean) * mask_lr(p, q);
          v_lr[p * n + q] = v;
          norm2 += (double)v * v;
        }

      const float norm = (float)std::sqrt(norm2);
      if (norm <= norm_min) continue;

      // paired high-resolution atom, scaled by the low-resolution norm so
      // that the coefficients computed on the LR atoms apply directly
      const int I0 = i0 * factor;
      const int J0 = j0 * factor;

      double mean_hr = 0.0;
      for (int q = 0; q < nh; q++)
        for (int p = 0; p < nh; p++)
          mean_hr += exemplar_hr(I0 + p, J0 + q);
      mean_hr /= (double)dict.dim_hr;

      for (int q = 0; q < nh; q++)
        for (int p = 0; p < nh; p++)
          v_hr[p * nh + q] = (float)(exemplar_hr(I0 + p, J0 + q) - mean_hr) *
                             mask_hr(p, q) / norm;

      for (float v : v_lr)
        dict.atoms_lr.push_back(v / norm);
      dict.atoms_hr.insert(dict.atoms_hr.end(), v_hr.begin(), v_hr.end());
      dict.n_atoms++;
    }

  return dict;
}

// --- sparse coding (orthogonal matching pursuit)

struct SparseCode
{
  std::vector<int>   atoms;
  std::vector<float> coeffs;
};

double dot(const float *a, const float *b, int n)
{
  double s = 0.0;
  for (int k = 0; k < n; k++)
    s += (double)a[k] * b[k];
  return s;
}

// Solve the small symmetric positive definite system G c = b in place
// (Cholesky, G row-major k x k, solution returned in b). Returns false if G
// is not positive definite.
bool solve_spd(std::vector<double> &g, std::vector<double> &b, int k)
{
  for (int j = 0; j < k; j++)
  {
    double s = g[j * k + j];
    for (int p = 0; p < j; p++)
      s -= g[j * k + p] * g[j * k + p];
    if (s <= 1e-12) return false;

    double ljj = std::sqrt(s);
    g[j * k + j] = ljj;

    for (int i = j + 1; i < k; i++)
    {
      double t = g[i * k + j];
      for (int p = 0; p < j; p++)
        t -= g[i * k + p] * g[j * k + p];
      g[i * k + j] = t / ljj;
    }
  }

  // L y = b
  for (int i = 0; i < k; i++)
  {
    double t = b[i];
    for (int p = 0; p < i; p++)
      t -= g[i * k + p] * b[p];
    b[i] = t / g[i * k + i];
  }

  // L^T c = y
  for (int i = k - 1; i >= 0; i--)
  {
    double t = b[i];
    for (int p = i + 1; p < k; p++)
      t -= g[p * k + i] * b[p];
    b[i] = t / g[i * k + i];
  }

  return true;
}

// Greedy atom selection with a least-squares refit on the support after each
// selection. `residual`, `gram` and `rhs` are scratch buffers.
void omp_encode(const std::vector<float> &x,
                const PatchDictionary    &dict,
                int                       sparsity,
                SparseCode               &code,
                std::vector<float>       &residual,
                std::vector<double>      &gram,
                std::vector<double>      &rhs)
{
  code.atoms.clear();
  code.coeffs.clear();

  const int    dim = dict.dim_lr;
  const double x_norm2 = dot(x.data(), x.data(), dim);

  // flat patch: nothing to code, only the mean will be used
  if (x_norm2 <= 0.0) return;

  residual = x;

  for (int t = 0; t < sparsity; t++)
  {
    // atom most correlated with the current residual
    int    best = -1;
    double best_abs = 0.0;

    for (int a = 0; a < dict.n_atoms; a++)
    {
      if (std::find(code.atoms.begin(), code.atoms.end(), a) !=
          code.atoms.end())
        continue;

      double c = dot(residual.data(), &dict.atoms_lr[(size_t)a * dim], dim);
      if (std::abs(c) > best_abs)
      {
        best_abs = std::abs(c);
        best = a;
      }
    }

    if (best < 0 || best_abs <= 1e-6 * std::sqrt(x_norm2)) break;

    code.atoms.push_back(best);
    const int k = (int)code.atoms.size();

    // least-squares coefficients on the support (k is tiny, the Gram matrix
    // is simply rebuilt)
    gram.assign((size_t)k * k, 0.0);
    rhs.assign(k, 0.0);

    for (int i = 0; i < k; i++)
    {
      const float *di = &dict.atoms_lr[(size_t)code.atoms[i] * dim];
      rhs[i] = dot(di, x.data(), dim);
      for (int j = 0; j <= i; j++)
      {
        const float *dj = &dict.atoms_lr[(size_t)code.atoms[j] * dim];
        double       g = dot(di, dj, dim);
        gram[i * k + j] = g;
        gram[j * k + i] = g;
      }
    }

    if (!solve_spd(gram, rhs, k))
    {
      // the new atom is (numerically) dependent on the support: drop it and
      // keep the previous fit
      code.atoms.pop_back();
      break;
    }

    code.coeffs.resize(k);
    for (int i = 0; i < k; i++)
      code.coeffs[i] = (float)rhs[i];

    // residual = x - sum_i c_i d_i
    residual = x;
    for (int i = 0; i < k; i++)
    {
      const float *di = &dict.atoms_lr[(size_t)code.atoms[i] * dim];
      const float  c = code.coeffs[i];
      for (int p = 0; p < dim; p++)
        residual[p] -= c * di[p];
    }

    if (dot(residual.data(), residual.data(), dim) <= 1e-10 * x_norm2) break;
  }
}

} // namespace

// --- super-resolution

Array terrain_super_resolution(const Array &array,
                               const Array &exemplar,
                               int          factor,
                               int          patch_size,
                               int          analysis_stride,
                               int          synthesis_stride,
                               int          sparsity)
{
  if (!validate_non_empty(array) || !validate_non_empty(exemplar))
    return Array();

  if (factor < 1)
  {
    hmap::log::warn("terrain_super_resolution: factor must be >= 1 (got {})",
                    factor);
    return Array();
  }

  if (patch_size < 2)
  {
    hmap::log::warn(
        "terrain_super_resolution: patch_size must be >= 2 (got {})",
        patch_size);
    return Array();
  }

  if (analysis_stride < 1 || analysis_stride > patch_size ||
      synthesis_stride < 1 || synthesis_stride > patch_size)
  {
    hmap::log::warn("terrain_super_resolution: strides must be in [1, "
                    "patch_size = {}] (got analysis_stride = {}, "
                    "synthesis_stride = {})",
                    patch_size,
                    analysis_stride,
                    synthesis_stride);
    return Array();
  }

  if (sparsity < 1)
  {
    hmap::log::warn("terrain_super_resolution: sparsity must be >= 1 (got {})",
                    sparsity);
    return Array();
  }

  const int n = patch_size;
  const int nh = patch_size * factor;

  // --- dictionary from the exemplar

  Array exemplar_lr = downsample_block_mean(exemplar, factor);

  if (exemplar_lr.shape.x < n || exemplar_lr.shape.y < n)
  {
    hmap::log::warn("terrain_super_resolution: exemplar downsampled by {} is "
                    "{}x{}, smaller than the patch size {}",
                    factor,
                    exemplar_lr.shape.x,
                    exemplar_lr.shape.y,
                    n);
    return Array();
  }

  Array mask_lr = build_patch_mask(n);
  Array mask_hr = build_patch_mask(nh);

  PatchDictionary dict = build_dictionary(exemplar,
                                          exemplar_lr,
                                          factor,
                                          n,
                                          analysis_stride,
                                          mask_lr,
                                          mask_hr);

  if (dict.n_atoms == 0)
  {
    hmap::log::warn("terrain_super_resolution: no usable patch in the "
                    "exemplar (flat exemplar?)");
    return Array();
  }

  hmap::log::trace("terrain_super_resolution: {} atoms", dict.n_atoms);

  // --- sparse coding of the input patches (on a padded input, so that the
  // --- borders are covered by full patches)

  Array padded = pad_replicate(array, n);

  const int e1 = (padded.shape.x - n) / synthesis_stride + 1;
  const int e2 = (padded.shape.y - n) / synthesis_stride + 1;
  const int n_patches = e1 * e2;

  std::vector<SparseCode> codes(n_patches);
  std::vector<float>      means(n_patches);

#pragma omp parallel
  {
    std::vector<float>  x(dict.dim_lr);
    std::vector<float>  residual(dict.dim_lr);
    std::vector<double> gram;
    std::vector<double> rhs;

#pragma omp for schedule(dynamic, 16)
    for (int k = 0; k < n_patches; k++)
    {
      const int i0 = (k % e1) * synthesis_stride;
      const int j0 = (k / e1) * synthesis_stride;

      double mean = 0.0;
      for (int q = 0; q < n; q++)
        for (int p = 0; p < n; p++)
          mean += padded(i0 + p, j0 + q);
      mean /= (double)dict.dim_lr;

      for (int q = 0; q < n; q++)
        for (int p = 0; p < n; p++)
          x[p * n + q] = (float)(padded(i0 + p, j0 + q) - mean) * mask_lr(p, q);

      means[k] = (float)mean;
      omp_encode(x, dict, sparsity, codes[k], residual, gram, rhs);
    }
  }

  // --- high-resolution synthesis: mask-weighted blend of the coded patches

  Array acc(padded.shape * factor);
  Array weight(padded.shape * factor);

  std::vector<float> v_hr(dict.dim_hr);

  for (int k = 0; k < n_patches; k++)
  {
    const int I0 = (k % e1) * synthesis_stride * factor;
    const int J0 = (k / e1) * synthesis_stride * factor;

    for (int q = 0; q < nh; q++)
      for (int p = 0; p < nh; p++)
        v_hr[p * nh + q] = means[k] * mask_hr(p, q);

    for (size_t l = 0; l < codes[k].atoms.size(); l++)
    {
      const float
          *atom = &dict.atoms_hr[(size_t)codes[k].atoms[l] * dict.dim_hr];
      const float c = codes[k].coeffs[l];
      for (int idx = 0; idx < dict.dim_hr; idx++)
        v_hr[idx] += c * atom[idx];
    }

    for (int q = 0; q < nh; q++)
      for (int p = 0; p < nh; p++)
      {
        acc(I0 + p, J0 + q) += v_hr[p * nh + q];
        weight(I0 + p, J0 + q) += mask_hr(p, q);
      }
  }

  // --- normalize and crop the padding

  glm::ivec2 shape_out = array.shape * factor;
  Array      out(shape_out);

  int n_holes = 0;

  for (int j = 0; j < shape_out.y; j++)
    for (int i = 0; i < shape_out.x; i++)
    {
      float w = weight(i + nh, j + nh);
      if (w > 0.f)
        out(i, j) = acc(i + nh, j + nh) / w;
      else
        n_holes++;
    }

  if (n_holes > 0)
  {
    // only happens with strides too large for the mask footprint
    Array fallback = array.resample_to_shape_bicubic(shape_out);

    for (int j = 0; j < shape_out.y; j++)
      for (int i = 0; i < shape_out.x; i++)
        if (weight(i + nh, j + nh) <= 0.f) out(i, j) = fallback(i, j);

    hmap::log::warn("terrain_super_resolution: {} pixels not covered by any "
                    "patch (synthesis_stride too large for patch_size), "
                    "filled with the upsampled input",
                    n_holes);
  }

  return out;
}

} // namespace hmap
