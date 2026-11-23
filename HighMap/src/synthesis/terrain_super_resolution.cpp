/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/blending.hpp"
#include "highmap/coord_frame.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/synthesis.hpp"
#include "highmap/virtual_array.hpp"

namespace hmap
{

namespace
{

// --- Helpers

// radial polynomial mask of the reference implementation: 1 at the patch
// center, decreasing to 0 before the corners, ((1 - (1 - 1/size) r^2)^+)^2
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

// block-mean pooling, so that the mean of a low-resolution patch equals the
// mean of the paired high-resolution patch
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

// clamp-to-edge padding by `nbuffer` cells on every side
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

// --- Dictionary

struct PatchDictionary
{
  int                dim_lr = 0;
  int                dim_hr = 0;
  int                n_atoms = 0;
  std::vector<float> atoms_lr; // n_atoms x dim_lr, unit L2 norm
  std::vector<float> atoms_hr; // n_atoms x dim_hr, scaled by the LR norm
};

void append_to_dictionary(const Array     &exemplar_hr,
                          const Array     &exemplar_lr,
                          int              factor,
                          int              patch_size,
                          int              stride,
                          const Array     &mask_lr,
                          const Array     &mask_hr,
                          float            hr_weight,
                          int              dict_start,
                          int              budget,
                          int             &added_count,
                          size_t          &seen_count,
                          std::mt19937    &rng,
                          PatchDictionary &dict)
{
  const int n = patch_size;
  const int nh = patch_size * factor;

  if (exemplar_lr.shape.x < n || exemplar_lr.shape.y < n) return;
  if (budget <= 0) return;

  dict.dim_lr = n * n;
  dict.dim_hr = nh * nh;

  const int d1 = (exemplar_lr.shape.x - n) / stride + 1;
  const int d2 = (exemplar_lr.shape.y - n) / stride + 1;

  // patches with no variation (relative to the exemplar range) are useless
  // as atoms and would blow up the normalization
  const float norm_min = 1e-6f * (exemplar_lr.max() - exemplar_lr.min());
  if (norm_min <= 0.f) return;

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

      seen_count++;

      int slot = -1;

      if (added_count < budget)
      {
        slot = dict_start + added_count;
        added_count++;
        dict.n_atoms++;
        dict.atoms_lr.resize((size_t)dict.n_atoms * dict.dim_lr);
        dict.atoms_hr.resize((size_t)dict.n_atoms * dict.dim_hr);
      }
      else
      {
        std::uniform_int_distribution<size_t> dist(0, seen_count - 1);
        size_t                                r = dist(rng);
        if (r >= (size_t)budget) continue;
        slot = dict_start + (int)r;
      }

      // paired high-resolution atom, scaled by the low-resolution norm and
      // source weight so that LR coefficients apply directly
      const int I0 = i0 * factor;
      const int J0 = j0 * factor;

      double mean_hr = 0.0;
      for (int q = 0; q < nh; q++)
        for (int p = 0; p < nh; p++)
          mean_hr += exemplar_hr(I0 + p, J0 + q);
      mean_hr /= (double)dict.dim_hr;

      const float hr_scale = hr_weight / norm;

      for (int q = 0; q < nh; q++)
        for (int p = 0; p < nh; p++)
          v_hr[p * nh + q] = (float)(exemplar_hr(I0 + p, J0 + q) - mean_hr) *
                             mask_hr(p, q) * hr_scale;

      float *dst_lr = &dict.atoms_lr[(size_t)slot * dict.dim_lr];
      float *dst_hr = &dict.atoms_hr[(size_t)slot * dict.dim_hr];

      for (int idx = 0; idx < dict.dim_lr; idx++)
        dst_lr[idx] = v_lr[idx] / norm;

      for (int idx = 0; idx < dict.dim_hr; idx++)
        dst_hr[idx] = v_hr[idx];
    }
}

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
  const int max_possible = std::max(0, d1) * std::max(0, d2);

  if (max_possible == 0) return dict;

  dict.atoms_lr.reserve((size_t)max_possible * dict.dim_lr);
  dict.atoms_hr.reserve((size_t)max_possible * dict.dim_hr);

  int          added = 0;
  size_t       seen = 0;
  std::mt19937 rng(0);

  append_to_dictionary(exemplar_hr,
                       exemplar_lr,
                       factor,
                       patch_size,
                       stride,
                       mask_lr,
                       mask_hr,
                       1.f,
                       0,
                       max_possible,
                       added,
                       seen,
                       rng,
                       dict);

  return dict;
}

void extract_source_chunk_atoms(const Array     &src_chunk,
                                glm::vec2        scale,
                                int              factor,
                                int              patch_size,
                                int              analysis_stride,
                                const Array     &mask_lr,
                                const Array     &mask_hr,
                                float            weight,
                                int              dict_start,
                                int              budget,
                                int             &added_count,
                                size_t          &seen_count,
                                std::mt19937    &rng,
                                PatchDictionary &dict)
{
  const int nh = patch_size * factor;

  // cap per-window resampled resolution to bound temporary memory (<= 1 MB)
  const int max_resampled_dim = std::max(512, 2 * nh);

  const int win_x = std::max(
      2,
      std::min(src_chunk.shape.x,
               int(std::floor(float(max_resampled_dim) /
                              std::max(1e-6f, scale.x)))));
  const int win_y = std::max(
      2,
      std::min(src_chunk.shape.y,
               int(std::floor(float(max_resampled_dim) /
                              std::max(1e-6f, scale.y)))));

  for (int j0 = 0; j0 < src_chunk.shape.y; j0 += win_y)
    for (int i0 = 0; i0 < src_chunk.shape.x; i0 += win_x)
    {
      const int i1 = std::min(src_chunk.shape.x, i0 + win_x);
      const int j1 = std::min(src_chunk.shape.y, j0 + win_y);

      const int sub_nx = i1 - i0;
      const int sub_ny = j1 - j0;

      const int matched_nx = int(std::round(float(sub_nx) * scale.x));
      const int matched_ny = int(std::round(float(sub_ny) * scale.y));

      if (matched_nx < nh || matched_ny < nh) continue;

      Array sub = (i0 == 0 && j0 == 0 && i1 == src_chunk.shape.x &&
                   j1 == src_chunk.shape.y)
                      ? src_chunk
                      : src_chunk.extract_slice(glm::ivec4(i0, i1, j0, j1));

      Array exemplar_hr = (sub.shape.x == matched_nx &&
                           sub.shape.y == matched_ny)
                              ? sub
                              : sub.resample_to_shape_bilinear(
                                    glm::ivec2(matched_nx, matched_ny));

      Array exemplar_lr = downsample_block_mean(exemplar_hr, factor);

      append_to_dictionary(exemplar_hr,
                           exemplar_lr,
                           factor,
                           patch_size,
                           analysis_stride,
                           mask_lr,
                           mask_hr,
                           weight,
                           dict_start,
                           budget,
                           added_count,
                           seen_count,
                           rng,
                           dict);
    }
}

// --- Sparse coding (orthogonal matching pursuit)

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

// --- High-resolution synthesis from a prebuilt dictionary

Array synthesize_from_dictionary(const Array           &array,
                                 const PatchDictionary &dict,
                                 int                    factor,
                                 int                    patch_size,
                                 int                    synthesis_stride,
                                 int                    sparsity,
                                 const Array           &mask_lr,
                                 const Array           &mask_hr)
{
  const int n = patch_size;
  const int nh = patch_size * factor;

  // sparse coding of the input patches (on a padded input, so that the
  // borders are covered by full patches)
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

  // high-resolution synthesis: mask-weighted blend of the coded patches
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

  // normalize and crop the padding
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

void apply_feature_transfer_on_buffer(Array                 &buf,
                                      const PatchDictionary &dict,
                                      int                    ir,
                                      float                  amplitude,
                                      bool         target_prefiltering,
                                      int          factor,
                                      int          patch_size,
                                      int          synthesis_stride,
                                      int          sparsity,
                                      const Array &mask_lr,
                                      const Array &mask_hr)
{
  if (buf.shape.x < factor || buf.shape.y < factor) return;

  Array buf_lr = downsample_block_mean(buf, factor);

  Array synth = synthesize_from_dictionary(buf_lr,
                                           dict,
                                           factor,
                                           patch_size,
                                           synthesis_stride,
                                           sparsity,
                                           mask_lr,
                                           mask_hr);

  if (synth.shape != buf.shape)
    synth = synth.resample_to_shape_bilinear(buf.shape);

  const int max_ir = std::max(1, std::min(buf.shape.x, buf.shape.y) / 2);
  const int effective_ir = std::clamp(ir, 1, max_ir);

  buf = transfer(synth, buf, effective_ir, amplitude, target_prefiltering);
}

bool validate_feature_transfer_params(glm::ivec2 target_shape,
                                      glm::vec2  target_size,
                                      size_t     n_sources,
                                      size_t     n_frames,
                                      int        ir,
                                      int        factor,
                                      int        patch_size,
                                      int        analysis_stride,
                                      int        synthesis_stride,
                                      int        sparsity,
                                      int        max_atoms)
{
  if (target_shape.x <= 0 || target_shape.y <= 0) return false;
  if (target_size.x <= 0.f || target_size.y <= 0.f) return false;
  if (n_sources == 0 || n_sources != n_frames) return false;
  if (ir < 1) return false;
  if (factor < 1) return false;
  if (patch_size < 2) return false;
  if (analysis_stride < 1 || analysis_stride > patch_size) return false;
  if (synthesis_stride < 1 || synthesis_stride > patch_size) return false;
  if (sparsity < 1) return false;
  if (max_atoms < 1) return false;

  return true;
}

} // namespace

// --- Super-resolution

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

  // --- Dictionary from the exemplar

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

  return synthesize_from_dictionary(array,
                                    dict,
                                    factor,
                                    patch_size,
                                    synthesis_stride,
                                    sparsity,
                                    mask_lr,
                                    mask_hr);
}

// --- Terrain feature transfer

void transfer_terrain_features(
    Array                                 &target,
    const CoordFrame                      &frame_target,
    const std::vector<const Array *>      &sources,
    const std::vector<const CoordFrame *> &frames_source,
    int                                    ir,
    float                                  amplitude,
    bool                                   target_prefiltering,
    const std::vector<float>              &source_weights,
    int                                    factor,
    int                                    patch_size,
    int                                    analysis_stride,
    int                                    synthesis_stride,
    int                                    sparsity,
    int                                    max_atoms,
    std::uint32_t                          seed)
{
  const glm::vec2 target_size = frame_target.get_size();

  if (!validate_feature_transfer_params(target.shape,
                                        target_size,
                                        sources.size(),
                                        frames_source.size(),
                                        ir,
                                        factor,
                                        patch_size,
                                        analysis_stride,
                                        synthesis_stride,
                                        sparsity,
                                        max_atoms))
    return;

  if (amplitude == 0.f && !target_prefiltering) return;

  // count valid sources to divide the max_atoms budget
  std::vector<size_t> active_indices;
  active_indices.reserve(sources.size());

  for (size_t k = 0; k < sources.size(); ++k)
  {
    if (!sources[k] || !frames_source[k]) continue;
    if (sources[k]->shape.x <= 0 || sources[k]->shape.y <= 0) continue;

    const glm::vec2 s_size = frames_source[k]->get_size();
    if (s_size.x <= 0.f || s_size.y <= 0.f) continue;

    const float w = k < source_weights.size() ? source_weights[k] : 1.f;
    if (w <= 0.f) continue;

    active_indices.push_back(k);
  }

  if (active_indices.empty()) return;

  const glm::vec2 delta_target(target_size.x / float(target.shape.x),
                               target_size.y / float(target.shape.y));

  const int   n = patch_size;
  const int   nh = patch_size * factor;
  const Array mask_lr = build_patch_mask(n);
  const Array mask_hr = build_patch_mask(nh);

  PatchDictionary dict;
  dict.dim_lr = n * n;
  dict.dim_hr = nh * nh;

  const int n_active = int(active_indices.size());

  for (int a_idx = 0; a_idx < n_active; ++a_idx)
  {
    const size_t      k = active_indices[a_idx];
    const Array      &src = *sources[k];
    const CoordFrame &frame_src = *frames_source[k];
    const float       w = k < source_weights.size() ? source_weights[k] : 1.f;

    const int budget_start = (a_idx * max_atoms) / n_active;
    const int budget_end = ((a_idx + 1) * max_atoms) / n_active;
    const int budget_k = std::max(1, budget_end - budget_start);

    const glm::vec2 src_size = frame_src.get_size();
    const glm::vec2 delta_src(src_size.x / float(src.shape.x),
                              src_size.y / float(src.shape.y));

    const glm::vec2 scale(delta_src.x / delta_target.x,
                          delta_src.y / delta_target.y);

    const int    dict_start = dict.n_atoms;
    int          added_k = 0;
    size_t       seen_k = 0;
    std::mt19937 rng(seed + std::uint32_t(k) * 2654435761u);

    extract_source_chunk_atoms(src,
                               scale,
                               factor,
                               patch_size,
                               analysis_stride,
                               mask_lr,
                               mask_hr,
                               w,
                               dict_start,
                               budget_k,
                               added_k,
                               seen_k,
                               rng,
                               dict);
  }

  if (dict.n_atoms == 0)
  {
    hmap::log::warn(
        "transfer_terrain_features: no usable patch in source layers");
    return;
  }

  apply_feature_transfer_on_buffer(target,
                                   dict,
                                   ir,
                                   amplitude,
                                   target_prefiltering,
                                   factor,
                                   patch_size,
                                   synthesis_stride,
                                   sparsity,
                                   mask_lr,
                                   mask_hr);
}

void transfer_terrain_features(Array            &target,
                               const CoordFrame &frame_target,
                               const Array      &source,
                               const CoordFrame &frame_source,
                               int               ir,
                               float             amplitude,
                               bool              target_prefiltering,
                               int               factor,
                               int               patch_size,
                               int               analysis_stride,
                               int               synthesis_stride,
                               int               sparsity,
                               int               max_atoms,
                               std::uint32_t     seed)
{
  std::vector<const Array *>      sources = {&source};
  std::vector<const CoordFrame *> frames_source = {&frame_source};

  transfer_terrain_features(target,
                            frame_target,
                            sources,
                            frames_source,
                            ir,
                            amplitude,
                            target_prefiltering,
                            {},
                            factor,
                            patch_size,
                            analysis_stride,
                            synthesis_stride,
                            sparsity,
                            max_atoms,
                            seed);
}

void transfer_terrain_features(
    VirtualArray                            &target,
    const CoordFrame                        &frame_target,
    const std::vector<const VirtualArray *> &sources,
    const std::vector<const CoordFrame *>   &frames_source,
    int                                      ir,
    float                                    amplitude,
    bool                                     target_prefiltering,
    const std::vector<float>                &source_weights,
    int                                      factor,
    int                                      patch_size,
    int                                      analysis_stride,
    int                                      synthesis_stride,
    int                                      sparsity,
    int                                      max_atoms,
    std::uint32_t                            seed,
    const ComputeMode                       &cm)
{
  const glm::vec2 target_size = frame_target.get_size();

  if (!validate_feature_transfer_params(target.shape,
                                        target_size,
                                        sources.size(),
                                        frames_source.size(),
                                        ir,
                                        factor,
                                        patch_size,
                                        analysis_stride,
                                        synthesis_stride,
                                        sparsity,
                                        max_atoms))
    return;

  if (amplitude == 0.f && !target_prefiltering) return;

  std::vector<size_t> active_indices;
  active_indices.reserve(sources.size());

  for (size_t k = 0; k < sources.size(); ++k)
  {
    if (!sources[k] || !frames_source[k]) continue;
    if (sources[k]->shape.x <= 0 || sources[k]->shape.y <= 0) continue;

    const glm::vec2 s_size = frames_source[k]->get_size();
    if (s_size.x <= 0.f || s_size.y <= 0.f) continue;

    const float w = k < source_weights.size() ? source_weights[k] : 1.f;
    if (w <= 0.f) continue;

    active_indices.push_back(k);
  }

  if (active_indices.empty()) return;

  const glm::vec2 delta_target(target_size.x / float(target.shape.x),
                               target_size.y / float(target.shape.y));

  const int   n = patch_size;
  const int   nh = patch_size * factor;
  const Array mask_lr = build_patch_mask(n);
  const Array mask_hr = build_patch_mask(nh);

  PatchDictionary dict;
  dict.dim_lr = n * n;
  dict.dim_hr = nh * nh;

  const int n_active = int(active_indices.size());
  const int max_direct_dim = std::max(512, 2 * nh);

  // --- Streamed dictionary extraction from source VirtualArrays

  for (int a_idx = 0; a_idx < n_active; ++a_idx)
  {
    const size_t        k = active_indices[a_idx];
    const VirtualArray &src = *sources[k];
    const CoordFrame   &frame_src = *frames_source[k];
    const float         w = k < source_weights.size() ? source_weights[k] : 1.f;

    const int budget_start = (a_idx * max_atoms) / n_active;
    const int budget_end = ((a_idx + 1) * max_atoms) / n_active;
    const int budget_k = std::max(1, budget_end - budget_start);

    const glm::vec2 src_size = frame_src.get_size();
    const glm::vec2 delta_src(src_size.x / float(src.shape.x),
                              src_size.y / float(src.shape.y));

    const glm::vec2 scale(delta_src.x / delta_target.x,
                          delta_src.y / delta_target.y);

    const int    dict_start = dict.n_atoms;
    int          added_k = 0;
    size_t       seen_k = 0;
    std::mt19937 rng(seed + std::uint32_t(k) * 2654435761u);

    if (src.shape.x <= max_direct_dim && src.shape.y <= max_direct_dim)
    {
      Array src_arr = src.to_array(ComputeMode{});

      extract_source_chunk_atoms(src_arr,
                                 scale,
                                 factor,
                                 patch_size,
                                 analysis_stride,
                                 mask_lr,
                                 mask_hr,
                                 w,
                                 dict_start,
                                 budget_k,
                                 added_k,
                                 seen_k,
                                 rng,
                                 dict);
    }
    else
    {
      ComputeMode seq_cm;
      seq_cm.mode = ForEachMode::VA_SEQUENTIAL;
      seq_cm.trim_storage = cm.trim_storage;

      for_each_tile(
          src,
          [&](const Array &in_tile, const TileRegion &)
          {
            extract_source_chunk_atoms(in_tile,
                                       scale,
                                       factor,
                                       patch_size,
                                       analysis_stride,
                                       mask_lr,
                                       mask_hr,
                                       w,
                                       dict_start,
                                       budget_k,
                                       added_k,
                                       seen_k,
                                       rng,
                                       dict);
          },
          seq_cm);
    }
  }

  if (dict.n_atoms == 0)
  {
    hmap::log::warn(
        "transfer_terrain_features: no usable patch in source layers");
    return;
  }

  // --- Per-tile feature synthesis and transfer on target

  for_each_tile(
      target,
      [&](Array &tile, const TileRegion &region, ICommunicator &comm)
      {
        apply_feature_transfer_on_buffer(tile,
                                         dict,
                                         ir,
                                         amplitude,
                                         target_prefiltering,
                                         factor,
                                         patch_size,
                                         synthesis_stride,
                                         sparsity,
                                         mask_lr,
                                         mask_hr);

        if (comm.size() > 1)
          comm.sync_overlap_buffers(tile, region, SyncOperation::SmoothBlend);
      },
      cm);
}

void transfer_terrain_features(VirtualArray       &target,
                               const CoordFrame   &frame_target,
                               const VirtualArray &source,
                               const CoordFrame   &frame_source,
                               int                 ir,
                               float               amplitude,
                               bool                target_prefiltering,
                               int                 factor,
                               int                 patch_size,
                               int                 analysis_stride,
                               int                 synthesis_stride,
                               int                 sparsity,
                               int                 max_atoms,
                               std::uint32_t       seed,
                               const ComputeMode  &cm)
{
  std::vector<const VirtualArray *> sources = {&source};
  std::vector<const CoordFrame *>   frames_source = {&frame_source};

  transfer_terrain_features(target,
                            frame_target,
                            sources,
                            frames_source,
                            ir,
                            amplitude,
                            target_prefiltering,
                            {},
                            factor,
                            patch_size,
                            analysis_stride,
                            synthesis_stride,
                            sparsity,
                            max_atoms,
                            seed,
                            cm);
}

} // namespace hmap
