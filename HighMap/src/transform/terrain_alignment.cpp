/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/transform/terrain_alignment.hpp"

namespace hmap
{

namespace
{

float sample_array_bilinear(const Array &arr, float rx, float ry)
{
  if (arr.shape.x <= 1 || arr.shape.y <= 1) return arr(0, 0);

  float xi = std::clamp(rx, 0.f, 1.f) * float(arr.shape.x - 1);
  float yi = std::clamp(ry, 0.f, 1.f) * float(arr.shape.y - 1);

  int i = std::min(int(xi), arr.shape.x - 2);
  int j = std::min(int(yi), arr.shape.y - 2);

  float u = xi - float(i);
  float v = yi - float(j);

  return arr.get_value_bilinear_at(i, j, u, v);
}

} // namespace

// --- Alignment helpers

namespace
{

struct AlignmentSample
{
  float x;
  float y;
  float diff;
};

TerrainAlignment solve_weighted_offset(
    const std::vector<AlignmentSample> &samples,
    const std::vector<float>           &weights)
{
  double sum_w = 0.0;
  double sum_wd = 0.0;

  for (size_t k = 0; k < samples.size(); ++k)
  {
    double w = double(weights[k]);
    sum_w += w;
    sum_wd += w * double(samples[k].diff);
  }

  TerrainAlignment res;
  if (sum_w > 1e-12) res.offset = float(sum_wd / sum_w);
  return res;
}

TerrainAlignment solve_weighted_plane(
    const std::vector<AlignmentSample> &samples,
    const std::vector<float>           &weights)
{
  double a00 = 0.0, a01 = 0.0, a02 = 0.0;
  double a11 = 0.0, a12 = 0.0, a22 = 0.0;
  double b0 = 0.0, b1 = 0.0, b2 = 0.0;

  for (size_t k = 0; k < samples.size(); ++k)
  {
    double w = double(weights[k]);
    double x = double(samples[k].x);
    double y = double(samples[k].y);
    double d = double(samples[k].diff);

    a00 += w;
    a01 += w * x;
    a02 += w * y;
    a11 += w * x * x;
    a12 += w * x * y;
    a22 += w * y * y;

    b0 += w * d;
    b1 += w * x * d;
    b2 += w * y * d;
  }

  // symmetric 3x3 determinant via cofactors
  double c00 = a11 * a22 - a12 * a12;
  double c01 = a02 * a12 - a01 * a22;
  double c02 = a01 * a12 - a02 * a11;
  double det = a00 * c00 + a01 * c01 + a02 * c02;

  if (std::abs(det) < 1e-12) return solve_weighted_offset(samples, weights);

  double c11 = a00 * a22 - a02 * a02;
  double c12 = a01 * a02 - a00 * a12;
  double c22 = a00 * a11 - a01 * a01;

  TerrainAlignment res;
  res.offset = float((c00 * b0 + c01 * b1 + c02 * b2) / det);
  res.slope_x = float((c01 * b0 + c11 * b1 + c12 * b2) / det);
  res.slope_y = float((c02 * b0 + c12 * b1 + c22 * b2) / det);
  return res;
}

} // namespace

// --- Alignment Public API

TerrainAlignment compute_terrain_alignment(const Array      &lower,
                                           const CoordFrame &frame_lower,
                                           const Array      &upper,
                                           const CoordFrame &frame_upper,
                                           TerrainAlignMode  mode,
                                           float             huber_delta,
                                           int               max_iterations)
{
  if (!frame_lower.intersects(frame_upper)) return TerrainAlignment{};

  const AffineTransform2D xform = frame_upper.compute_relative_transform(
      frame_lower);

  const float inv_nx = upper.shape.x > 1 ? 1.f / float(upper.shape.x - 1) : 0.f;
  const float inv_ny = upper.shape.y > 1 ? 1.f / float(upper.shape.y - 1) : 0.f;

  std::vector<AlignmentSample> samples;
  samples.reserve(size_t(upper.shape.x) * size_t(upper.shape.y));

  for (int j = 0; j < upper.shape.y; ++j)
  {
    const float ry_upper = float(j) * inv_ny;

    for (int i = 0; i < upper.shape.x; ++i)
    {
      const float     rx_upper = float(i) * inv_nx;
      const glm::vec2 r_lower = xform.apply(rx_upper, ry_upper);

      if (!frame_lower.is_unit_point_within(r_lower.x, r_lower.y)) continue;

      const float lower_val = sample_array_bilinear(lower,
                                                    r_lower.x,
                                                    r_lower.y);
      samples.push_back(
          {rx_upper - 0.5f, ry_upper - 0.5f, lower_val - upper(i, j)});
    }
  }

  if (samples.empty()) return TerrainAlignment{};

  std::vector<float> weights(samples.size(), 1.f);

  const bool is_plane = (mode == TerrainAlignMode::TAM_PLANE ||
                         mode == TerrainAlignMode::TAM_ROBUST_PLANE);
  const bool is_robust = (mode == TerrainAlignMode::TAM_ROBUST_OFFSET ||
                          mode == TerrainAlignMode::TAM_ROBUST_PLANE);

  TerrainAlignment fit = is_plane ? solve_weighted_plane(samples, weights)
                                  : solve_weighted_offset(samples, weights);

  if (is_robust && huber_delta > 0.f)
  {
    for (int iter = 0; iter < max_iterations; ++iter)
    {
      for (size_t k = 0; k < samples.size(); ++k)
      {
        const float pred = fit.offset + fit.slope_x * samples[k].x +
                           fit.slope_y * samples[k].y;
        const float r = std::abs(samples[k].diff - pred);
        weights[k] = (r <= huber_delta) ? 1.f : (huber_delta / r);
      }

      fit = is_plane ? solve_weighted_plane(samples, weights)
                     : solve_weighted_offset(samples, weights);
    }
  }

  return fit;
}

TerrainAlignment compute_terrain_alignment(const VirtualArray &lower,
                                           const CoordFrame   &frame_lower,
                                           const VirtualArray &upper,
                                           const CoordFrame   &frame_upper,
                                           TerrainAlignMode    mode,
                                           float               huber_delta,
                                           int                 max_iterations,
                                           int                 coarse_size)
{
  if (!frame_lower.intersects(frame_upper)) return TerrainAlignment{};

  const glm::ivec2 lower_shape(
      std::min(lower.shape.x, std::max(2, coarse_size)),
      std::min(lower.shape.y, std::max(2, coarse_size)));
  const glm::ivec2 upper_shape(
      std::min(upper.shape.x, std::max(2, coarse_size)),
      std::min(upper.shape.y, std::max(2, coarse_size)));

  Array lower_coarse = lower.to_array(lower_shape, ComputeMode{});
  Array upper_coarse = upper.to_array(upper_shape, ComputeMode{});

  return compute_terrain_alignment(lower_coarse,
                                   frame_lower,
                                   upper_coarse,
                                   frame_upper,
                                   mode,
                                   huber_delta,
                                   max_iterations);
}

void apply_terrain_alignment(Array &upper, const TerrainAlignment &alignment)
{
  const float inv_nx = upper.shape.x > 1 ? 1.f / float(upper.shape.x - 1) : 0.f;
  const float inv_ny = upper.shape.y > 1 ? 1.f / float(upper.shape.y - 1) : 0.f;

  for (int j = 0; j < upper.shape.y; ++j)
  {
    const float y = float(j) * inv_ny - 0.5f;

    for (int i = 0; i < upper.shape.x; ++i)
    {
      const float x = float(i) * inv_nx - 0.5f;
      upper(i, j) += alignment.offset + alignment.slope_x * x +
                     alignment.slope_y * y;
    }
  }
}

void apply_terrain_alignment(VirtualArray           &upper,
                             const TerrainAlignment &alignment,
                             const ComputeMode      &cm)
{
  const float inv_nx = upper.shape.x > 1 ? 1.f / float(upper.shape.x - 1) : 0.f;
  const float inv_ny = upper.shape.y > 1 ? 1.f / float(upper.shape.y - 1) : 0.f;

  for_each_tile(
      upper,
      [&](Array &tile, const TileRegion &region)
      {
        const glm::ivec2 offset = region.global_offset();

        for (int j = 0; j < tile.shape.y; ++j)
        {
          const float y = float(offset.y + j) * inv_ny - 0.5f;

          for (int i = 0; i < tile.shape.x; ++i)
          {
            const float x = float(offset.x + i) * inv_nx - 0.5f;
            tile(i, j) += alignment.offset + alignment.slope_x * x +
                          alignment.slope_y * y;
          }
        }
      },
      cm);
}

void align_terrains(const Array      &lower,
                    const CoordFrame &frame_lower,
                    Array            &upper,
                    const CoordFrame &frame_upper,
                    TerrainAlignMode  mode,
                    float             huber_delta,
                    int               max_iterations)
{
  TerrainAlignment alignment = compute_terrain_alignment(lower,
                                                         frame_lower,
                                                         upper,
                                                         frame_upper,
                                                         mode,
                                                         huber_delta,
                                                         max_iterations);
  apply_terrain_alignment(upper, alignment);
}

void align_terrains(const VirtualArray &lower,
                    const CoordFrame   &frame_lower,
                    VirtualArray       &upper,
                    const CoordFrame   &frame_upper,
                    TerrainAlignMode    mode,
                    float               huber_delta,
                    int                 max_iterations,
                    int                 coarse_size,
                    const ComputeMode  &cm)
{
  TerrainAlignment alignment = compute_terrain_alignment(lower,
                                                         frame_lower,
                                                         upper,
                                                         frame_upper,
                                                         mode,
                                                         huber_delta,
                                                         max_iterations,
                                                         coarse_size);
  apply_terrain_alignment(upper, alignment, cm);
}

} // namespace hmap
