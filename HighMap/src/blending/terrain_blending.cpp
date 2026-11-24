/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/blending.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/math/core.hpp"
#include "highmap/range.hpp"

namespace hmap
{

namespace
{

// --- Helpers

float blend_value(float            lower_val,
                  float            upper_val,
                  float            rx_upper,
                  float            ry_upper,
                  TerrainBlendMode mode,
                  float            param)
{
  switch (mode)
  {
  case TerrainBlendMode::TBM_ADD: return lower_val + upper_val;

  case TerrainBlendMode::TBM_AVERAGE: return 0.5f * (lower_val + upper_val);

  case TerrainBlendMode::TBM_DISTANCE:
  {
    float d = std::min({rx_upper, 1.f - rx_upper, ry_upper, 1.f - ry_upper});
    float w = param > 0.f ? smoothstep3(std::clamp(d / param, 0.f, 1.f)) : 1.f;
    return lerp(lower_val, upper_val, w);
  }

  case TerrainBlendMode::TBM_MAX: return std::max(lower_val, upper_val);

  case TerrainBlendMode::TBM_MAX_SMOOTH:
    return maximum_smooth(lower_val, upper_val, param);

  case TerrainBlendMode::TBM_MIN: return std::min(lower_val, upper_val);

  case TerrainBlendMode::TBM_MIN_SMOOTH:
    return minimum_smooth(lower_val, upper_val, param);

  case TerrainBlendMode::TBM_SUBTRACT: return lower_val - upper_val;

  default: break;
  }

  return lower_val;
}

float edge_distance_weight(float rx,
                           float ry,
                           float feather_width,
                           bool  quintic = false)
{
  if (feather_width <= 0.f) return 1.f;

  float d = std::min({rx, 1.f - rx, ry, 1.f - ry});
  float t = std::clamp(d / feather_width, 0.f, 1.f);
  return quintic ? smoothstep5(t) : smoothstep3(t);
}

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

struct BoundarySample
{
  float x;
  float y;
  float diff;
};

template <typename LowerSampler, typename UpperSampler>
std::vector<BoundarySample> sample_overlap_boundary(
    const CoordFrame &frame_lower,
    const CoordFrame &frame_upper,
    LowerSampler    &&sample_lower,
    UpperSampler    &&sample_upper,
    int               boundary_samples)
{
  const int n = std::max(4, boundary_samples);

  const AffineTransform2D xform_u2l = frame_upper.compute_relative_transform(
      frame_lower);

  std::vector<BoundarySample> samples;
  samples.reserve(size_t(4) * size_t(n));

  auto add_pt = [&](float rx, float ry)
  {
    const glm::vec2 r_lower = xform_u2l.apply(rx, ry);
    if (!frame_lower.is_unit_point_within(r_lower.x, r_lower.y)) return;

    const float lower_val = sample_lower(r_lower.x, r_lower.y);
    const float upper_val = sample_upper(rx, ry);
    samples.push_back({rx, ry, lower_val - upper_val});
  };

  // bottom edge (y = 0, x: 0 -> 1)
  for (int k = 0; k < n; ++k)
    add_pt(float(k) / float(n), 0.f);

  // right edge (x = 1, y: 0 -> 1)
  for (int k = 0; k < n; ++k)
    add_pt(1.f, float(k) / float(n));

  // top edge (y = 1, x: 1 -> 0)
  for (int k = 0; k < n; ++k)
    add_pt(1.f - float(k) / float(n), 1.f);

  // left edge (x = 0, y: 1 -> 0)
  for (int k = 0; k < n; ++k)
    add_pt(0.f, 1.f - float(k) / float(n));

  return samples;
}

float evaluate_mvc_membrane(const std::vector<BoundarySample> &samples,
                            float                              rx,
                            float                              ry)
{
  if (samples.empty()) return 0.f;
  if (samples.size() == 1) return samples[0].diff;

  const size_t n = samples.size();
  double       sum_w = 0.0;
  double       sum_wd = 0.0;

  for (size_t k = 0; k < n; ++k)
  {
    const float dx = samples[k].x - rx;
    const float dy = samples[k].y - ry;
    const float dist = std::hypot(dx, dy);

    if (dist < 1e-6f) return samples[k].diff;

    const BoundarySample &prev = samples[(k + n - 1) % n];
    const BoundarySample &next = samples[(k + 1) % n];

    const float dx_prev = prev.x - rx;
    const float dy_prev = prev.y - ry;
    const float dx_next = next.x - rx;
    const float dy_next = next.y - ry;

    const float d_prev = std::hypot(dx_prev, dy_prev);
    const float d_next = std::hypot(dx_next, dy_next);

    if (d_prev < 1e-6f) return prev.diff;
    if (d_next < 1e-6f) return next.diff;

    auto half_tan =
        [](float ax, float ay, float bx, float by, float la, float lb) -> double
    {
      const double dot = double(ax) * double(bx) + double(ay) * double(by);
      const double cross = double(ax) * double(by) - double(ay) * double(bx);
      const double denom = double(la) * double(lb) + dot;
      if (std::abs(denom) < 1e-12) return 0.0;
      return std::abs(cross) / denom;
    };

    const double w = (half_tan(dx_prev, dy_prev, dx, dy, d_prev, dist) +
                      half_tan(dx, dy, dx_next, dy_next, dist, d_next)) /
                     double(dist);

    sum_w += w;
    sum_wd += w * double(samples[k].diff);
  }

  return sum_w > 1e-12 ? float(sum_wd / sum_w) : 0.f;
}

void solve_poisson_defect_red_black(Array                      &u,
                                    const Array                &rhs,
                                    const std::vector<uint8_t> &is_interior,
                                    int                         iterations,
                                    float                       omega)
{
  const int nx = u.shape.x;
  const int ny = u.shape.y;

  for (int it = 0; it < iterations; ++it)
  {
    for (int color = 0; color < 2; ++color)
    {
      for (int j = 1; j < ny - 1; ++j)
      {
        for (int i = 1; i < nx - 1; ++i)
        {
          if (((i + j) & 1) != color) continue;

          const size_t idx = size_t(j) * size_t(nx) + size_t(i);
          if (!is_interior[idx]) continue;

          const float sum_nb = u(i - 1, j) + u(i + 1, j) + u(i, j - 1) +
                               u(i, j + 1);
          const float u_star = 0.25f * (sum_nb - rhs(i, j));
          u(i, j) = u(i, j) + omega * (u_star - u(i, j));
        }
      }
    }
  }
}

void build_poisson_patch_inputs(const Array                &lower,
                                const Array                &upper_proj,
                                const std::vector<uint8_t> &in_mask,
                                const Array                &feather_weights,
                                bool                        mixed_gradients,
                                Array                      &u,
                                Array                      &rhs,
                                std::vector<uint8_t>       &is_interior)
{
  const int nx = lower.shape.x;
  const int ny = lower.shape.y;

  u = Array(lower.shape, 0.f);
  rhs = Array(lower.shape, 0.f);
  is_interior.assign(size_t(nx) * size_t(ny), 0);

  for (int j = 0; j < ny; ++j)
  {
    for (int i = 0; i < nx; ++i)
    {
      const size_t idx = size_t(j) * size_t(nx) + size_t(i);
      if (!in_mask[idx]) continue;

      const bool interior = (i > 0 && i < nx - 1 && j > 0 && j < ny - 1) &&
                            in_mask[idx - 1] && in_mask[idx + 1] &&
                            in_mask[idx - size_t(nx)] &&
                            in_mask[idx + size_t(nx)];

      if (interior)
      {
        is_interior[idx] = 1;

        if (mixed_gradients)
        {
          // divergence of (w * grad(upper) + (1 - w) * grad(lower)) -
          // laplacian(upper) = div((1 - w) * (grad(lower) - grad(upper)))
          auto diff_at = [&](int x, int y)
          { return lower(x, y) - upper_proj(x, y); };

          const float w_c = 1.f - feather_weights(i, j);
          const float w_e = 0.5f * (w_c + (1.f - feather_weights(i + 1, j)));
          const float w_w = 0.5f * (w_c + (1.f - feather_weights(i - 1, j)));
          const float w_n = 0.5f * (w_c + (1.f - feather_weights(i, j + 1)));
          const float w_s = 0.5f * (w_c + (1.f - feather_weights(i, j - 1)));

          const float d_c = diff_at(i, j);
          rhs(i, j) = w_e * (diff_at(i + 1, j) - d_c) +
                      w_w * (diff_at(i - 1, j) - d_c) +
                      w_n * (diff_at(i, j + 1) - d_c) +
                      w_s * (diff_at(i, j - 1) - d_c);
        }
      }
      else
      {
        // Dirichlet boundary condition: u = lower - upper
        u(i, j) = lower(i, j) - upper_proj(i, j);
      }
    }
  }
}

} // namespace

// --- Public API

void blend_terrains(Array            &lower,
                    const CoordFrame &frame_lower,
                    const Array      &upper,
                    const CoordFrame &frame_upper,
                    TerrainBlendMode  mode,
                    float             param)
{
  switch (mode)
  {
  case TerrainBlendMode::TBM_HEIGHT:
    blend_terrains_height(lower, frame_lower, upper, frame_upper, param);
    return;

  case TerrainBlendMode::TBM_MULTISCALE:
    blend_terrains_multiscale(lower, frame_lower, upper, frame_upper, param);
    return;

  case TerrainBlendMode::TBM_POISSON:
    blend_terrains_poisson(lower,
                           frame_lower,
                           upper,
                           frame_upper,
                           400,
                           param,
                           param > 0.f);
    return;

  case TerrainBlendMode::TBM_SEAMLESS:
    blend_terrains_seamless(lower, frame_lower, upper, frame_upper, param);
    return;

  case TerrainBlendMode::TBM_SLOPE:
    blend_terrains_slope(lower, frame_lower, upper, frame_upper, param);
    return;

  default: break;
  }

  if (!frame_lower.intersects(frame_upper)) return;

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;

  for (int j = 0; j < lower.shape.y; ++j)
  {
    const float ry_lower = float(j) * inv_ny;

    for (int i = 0; i < lower.shape.x; ++i)
    {
      const float     rx_lower = float(i) * inv_nx;
      const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

      if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y)) continue;

      const float upper_val = sample_array_bilinear(upper,
                                                    r_upper.x,
                                                    r_upper.y);
      lower(i, j) = blend_value(lower(i, j),
                                upper_val,
                                r_upper.x,
                                r_upper.y,
                                mode,
                                param);
    }
  }
}

void blend_terrains(VirtualArray       &lower,
                    const CoordFrame   &frame_lower,
                    const VirtualArray &upper,
                    const CoordFrame   &frame_upper,
                    TerrainBlendMode    mode,
                    float               param,
                    const ComputeMode  &cm)
{
  switch (mode)
  {
  case TerrainBlendMode::TBM_HEIGHT:
    blend_terrains_height(lower,
                          frame_lower,
                          upper,
                          frame_upper,
                          param,
                          0.f,
                          0.2f,
                          cm);
    return;

  case TerrainBlendMode::TBM_MULTISCALE:
    blend_terrains_multiscale(lower,
                              frame_lower,
                              upper,
                              frame_upper,
                              param,
                              3,
                              16,
                              1.f,
                              128,
                              cm);
    return;

  case TerrainBlendMode::TBM_POISSON:
    blend_terrains_poisson(lower,
                           frame_lower,
                           upper,
                           frame_upper,
                           400,
                           param,
                           param > 0.f,
                           128,
                           25,
                           cm);
    return;

  case TerrainBlendMode::TBM_SEAMLESS:
    blend_terrains_seamless(lower,
                            frame_lower,
                            upper,
                            frame_upper,
                            param,
                            64,
                            cm);
    return;

  case TerrainBlendMode::TBM_SLOPE:
    blend_terrains_slope(lower,
                         frame_lower,
                         upper,
                         frame_upper,
                         param,
                         0.1f,
                         2,
                         cm);
    return;

  default: break;
  }

  if (!frame_lower.intersects(frame_upper)) return;

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const glm::vec4 upper_bbox = {0.f, 1.f, 0.f, 1.f};
  const float     upper_dx = upper_bbox.y - upper_bbox.x;
  const float     upper_dy = upper_bbox.w - upper_bbox.z;

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;

  for_each_tile(
      lower,
      [&](Array &tile, const UvTile &region)
      {
        const glm::ivec2 offset = region.global_offset();

        for (int j = 0; j < tile.shape.y; ++j)
        {
          const int   gj = offset.y + j;
          const float ry_lower = float(gj) * inv_ny;

          for (int i = 0; i < tile.shape.x; ++i)
          {
            const int       gi = offset.x + i;
            const float     rx_lower = float(gi) * inv_nx;
            const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

            if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y))
              continue;

            const float ux = upper_bbox.x + r_upper.x * upper_dx;
            const float uy = upper_bbox.z + r_upper.y * upper_dy;
            const float upper_val = upper.get_bilinear(ux, uy);

            tile(i, j) = blend_value(tile(i, j),
                                     upper_val,
                                     r_upper.x,
                                     r_upper.y,
                                     mode,
                                     param);
          }
        }
      },
      cm);
}

void blend_terrains_height(Array            &lower,
                           const CoordFrame &frame_lower,
                           const Array      &upper,
                           const CoordFrame &frame_upper,
                           float             feather_width,
                           float             height_offset,
                           float             height_range)
{
  if (!frame_lower.intersects(frame_upper)) return;

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;
  const float inv_range = 1.f / std::max(height_range, 1e-6f);

  for (int j = 0; j < lower.shape.y; ++j)
  {
    const float ry_lower = float(j) * inv_ny;

    for (int i = 0; i < lower.shape.x; ++i)
    {
      const float     rx_lower = float(i) * inv_nx;
      const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

      if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y)) continue;

      const float upper_val = sample_array_bilinear(upper,
                                                    r_upper.x,
                                                    r_upper.y);
      const float lower_val = lower(i, j);

      const float w_d = edge_distance_weight(r_upper.x,
                                             r_upper.y,
                                             feather_width);
      const float t_h = std::clamp(
          ((upper_val - lower_val) - height_offset) * inv_range + 0.5f,
          0.f,
          1.f);
      const float w_h = smoothstep3(t_h);

      lower(i, j) = lerp(lower_val, upper_val, w_d * w_h);
    }
  }
}

void blend_terrains_height(VirtualArray       &lower,
                           const CoordFrame   &frame_lower,
                           const VirtualArray &upper,
                           const CoordFrame   &frame_upper,
                           float               feather_width,
                           float               height_offset,
                           float               height_range,
                           const ComputeMode  &cm)
{
  if (!frame_lower.intersects(frame_upper)) return;

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const glm::vec4 upper_bbox = {0.f, 1.f, 0.f, 1.f};
  const float     upper_dx = upper_bbox.y - upper_bbox.x;
  const float     upper_dy = upper_bbox.w - upper_bbox.z;

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;
  const float inv_range = 1.f / std::max(height_range, 1e-6f);

  for_each_tile(
      lower,
      [&](Array &tile, const UvTile &region)
      {
        const glm::ivec2 offset = region.global_offset();

        for (int j = 0; j < tile.shape.y; ++j)
        {
          const float ry_lower = float(offset.y + j) * inv_ny;

          for (int i = 0; i < tile.shape.x; ++i)
          {
            const float     rx_lower = float(offset.x + i) * inv_nx;
            const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

            if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y))
              continue;

            const float ux = upper_bbox.x + r_upper.x * upper_dx;
            const float uy = upper_bbox.z + r_upper.y * upper_dy;
            const float upper_val = upper.get_bilinear(ux, uy);
            const float lower_val = tile(i, j);

            const float w_d = edge_distance_weight(r_upper.x,
                                                   r_upper.y,
                                                   feather_width);
            const float t_h = std::clamp(
                ((upper_val - lower_val) - height_offset) * inv_range + 0.5f,
                0.f,
                1.f);
            const float w_h = smoothstep3(t_h);

            tile(i, j) = lerp(lower_val, upper_val, w_d * w_h);
          }
        }
      },
      cm);
}

void blend_terrains_multiscale(Array            &lower,
                               const CoordFrame &frame_lower,
                               const Array      &upper,
                               const CoordFrame &frame_upper,
                               float             feather_width,
                               int               n_levels,
                               int               base_radius,
                               float             detail_gain)
{
  if (!frame_lower.intersects(frame_upper)) return;

  const int levels = std::max(1, n_levels);
  const int nx = lower.shape.x;
  const int ny = lower.shape.y;

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const float inv_nx = nx > 1 ? 1.f / float(nx - 1) : 0.f;
  const float inv_ny = ny > 1 ? 1.f / float(ny - 1) : 0.f;

  Array                upper_proj = lower;
  std::vector<uint8_t> in_mask(size_t(nx) * size_t(ny), 0);
  Array                r_ux(lower.shape, 0.f);
  Array                r_uy(lower.shape, 0.f);

  for (int j = 0; j < ny; ++j)
  {
    const float ry_lower = float(j) * inv_ny;

    for (int i = 0; i < nx; ++i)
    {
      const float     rx_lower = float(i) * inv_nx;
      const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

      if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y)) continue;

      in_mask[size_t(j) * size_t(nx) + size_t(i)] = 1;
      upper_proj(i, j) = sample_array_bilinear(upper, r_upper.x, r_upper.y);
      r_ux(i, j) = r_upper.x;
      r_uy(i, j) = r_upper.y;
    }
  }

  // decompose into low-pass base and high-pass detail bands
  std::vector<Array> low_bands;
  std::vector<Array> high_bands;
  low_bands.reserve(size_t(levels) + 1);
  high_bands.reserve(size_t(levels) + 1);

  low_bands.push_back(lower);
  high_bands.push_back(upper_proj);

  for (int l = 0; l < levels; ++l)
  {
    const int ir = std::max(1, base_radius >> (levels - 1 - l));
    Array     l_next = low_bands.back();
    Array     h_next = high_bands.back();
    smooth_cpulse(l_next, ir);
    smooth_cpulse(h_next, ir);
    low_bands.push_back(std::move(l_next));
    high_bands.push_back(std::move(h_next));
  }

  // blend coarsest low-pass base with wide feathering
  Array out(lower.shape, 0.f);
  for (int j = 0; j < ny; ++j)
  {
    for (int i = 0; i < nx; ++i)
    {
      if (!in_mask[size_t(j) * size_t(nx) + size_t(i)])
      {
        out(i, j) = lower(i, j);
        continue;
      }

      const float w_base = edge_distance_weight(r_ux(i, j),
                                                r_uy(i, j),
                                                feather_width,
                                                true);
      out(i, j) = lerp(low_bands.back()(i, j), high_bands.back()(i, j), w_base);
    }
  }

  // add detail bands with progressively narrower feathering
  for (int l = levels - 1; l >= 0; --l)
  {
    const float band_feather = std::max(0.05f,
                                        feather_width *
                                            std::pow(0.65f, float(levels - l)));

    for (int j = 0; j < ny; ++j)
    {
      for (int i = 0; i < nx; ++i)
      {
        if (!in_mask[size_t(j) * size_t(nx) + size_t(i)]) continue;

        const float d_low = low_bands[l](i, j) - low_bands[l + 1](i, j);
        const float d_high = (high_bands[l](i, j) - high_bands[l + 1](i, j)) *
                             detail_gain;
        const float w_band = edge_distance_weight(r_ux(i, j),
                                                  r_uy(i, j),
                                                  band_feather,
                                                  true);

        out(i, j) += lerp(d_low, d_high, w_band);
      }
    }
  }

  lower = std::move(out);
}

void blend_terrains_multiscale(VirtualArray       &lower,
                               const CoordFrame   &frame_lower,
                               const VirtualArray &upper,
                               const CoordFrame   &frame_upper,
                               float               feather_width,
                               int                 n_levels,
                               int                 base_radius,
                               float               detail_gain,
                               int                 coarse_size,
                               const ComputeMode  &cm)
{
  if (!frame_lower.intersects(frame_upper)) return;

  const int levels = std::max(1, n_levels);

  // Level 1: coarse global proxy for low-frequency base
  const glm::ivec2 coarse_shape(
      std::min(lower.shape.x, std::max(4, coarse_size)),
      std::min(lower.shape.y, std::max(4, coarse_size)));

  Array lower_coarse = lower.to_array(coarse_shape, ComputeMode{});
  Array upper_coarse = upper.to_array(
      glm::ivec2(std::min(upper.shape.x, std::max(4, coarse_size)),
                 std::min(upper.shape.y, std::max(4, coarse_size))),
      ComputeMode{});

  const float scale_x = float(coarse_shape.x) /
                        float(std::max(1, lower.shape.x));
  const int coarse_radius = std::max(1, int(std::round(base_radius * scale_x)));

  Array lower_base_coarse = lower_coarse;
  smooth_cpulse(lower_base_coarse, coarse_radius);

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  // project upper_coarse onto lower_coarse grid
  Array upper_proj_coarse = lower_coarse;
  for (int j = 0; j < coarse_shape.y; ++j)
  {
    const float ry = coarse_shape.y > 1 ? float(j) / float(coarse_shape.y - 1)
                                        : 0.f;
    for (int i = 0; i < coarse_shape.x; ++i)
    {
      const float rx = coarse_shape.x > 1 ? float(i) / float(coarse_shape.x - 1)
                                          : 0.f;
      const glm::vec2 r_u = xform.apply(rx, ry);
      if (frame_upper.is_unit_point_within(r_u.x, r_u.y))
        upper_proj_coarse(i, j) = sample_array_bilinear(upper_coarse,
                                                        r_u.x,
                                                        r_u.y);
    }
  }

  Array upper_base_coarse = upper_proj_coarse;
  smooth_cpulse(upper_base_coarse, coarse_radius);

  // Level 2: tile-local detail blending + halo sync
  const glm::vec4 upper_bbox = {0.f, 1.f, 0.f, 1.f};
  const float     upper_dx = upper_bbox.y - upper_bbox.x;
  const float     upper_dy = upper_bbox.w - upper_bbox.z;

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;

  for_each_tile(
      lower,
      [&](Array &tile, const UvTile &region, ICommunicator &comm)
      {
        const glm::ivec2 offset = region.global_offset();
        const int        tx = tile.shape.x;
        const int        ty = tile.shape.y;

        Array                upper_proj = tile;
        std::vector<uint8_t> in_mask(size_t(tx) * size_t(ty), 0);
        Array                r_ux(tile.shape, 0.f);
        Array                r_uy(tile.shape, 0.f);

        bool any_inside = false;
        for (int j = 0; j < ty; ++j)
        {
          const float ry_lower = float(offset.y + j) * inv_ny;

          for (int i = 0; i < tx; ++i)
          {
            const float     rx_lower = float(offset.x + i) * inv_nx;
            const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

            if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y))
              continue;

            any_inside = true;
            in_mask[size_t(j) * size_t(tx) + size_t(i)] = 1;
            const float ux = upper_bbox.x + r_upper.x * upper_dx;
            const float uy = upper_bbox.z + r_upper.y * upper_dy;
            upper_proj(i, j) = upper.get_bilinear(ux, uy);
            r_ux(i, j) = r_upper.x;
            r_uy(i, j) = r_upper.y;
          }
        }

        if (any_inside)
        {
          // fine local high-pass decomposition
          const int fine_ir = std::max(
              1,
              std::min(base_radius, std::min(tx, ty) / 4));
          Array low_fine = tile;
          Array high_fine = upper_proj;
          smooth_cpulse(low_fine, fine_ir);
          smooth_cpulse(high_fine, fine_ir);

          const float fine_feather = feather_width *
                                     std::pow(0.5f, float(levels));

          for (int j = 0; j < ty; ++j)
          {
            const float ry_lower = float(offset.y + j) * inv_ny;

            for (int i = 0; i < tx; ++i)
            {
              if (!in_mask[size_t(j) * size_t(tx) + size_t(i)]) continue;

              const float rx_lower = float(offset.x + i) * inv_nx;

              const float base_l = sample_array_bilinear(lower_base_coarse,
                                                         rx_lower,
                                                         ry_lower);
              const float base_u = sample_array_bilinear(upper_base_coarse,
                                                         rx_lower,
                                                         ry_lower);

              const float w_base = edge_distance_weight(r_ux(i, j),
                                                        r_uy(i, j),
                                                        feather_width);
              const float blended_base = lerp(base_l, base_u, w_base);

              const float d_low = tile(i, j) - low_fine(i, j);
              const float d_high = (upper_proj(i, j) - high_fine(i, j)) *
                                   detail_gain;
              const float w_fine = edge_distance_weight(r_ux(i, j),
                                                        r_uy(i, j),
                                                        fine_feather);

              tile(i, j) = blended_base + lerp(d_low, d_high, w_fine);
            }
          }
        }

        if (comm.size() > 1)
          comm.sync_overlap_buffers(tile, region, SyncOperation::SmoothBlend);
      },
      cm);
}

void blend_terrains_poisson(Array            &lower,
                            const CoordFrame &frame_lower,
                            const Array      &upper,
                            const CoordFrame &frame_upper,
                            int               iterations,
                            float             feather_width,
                            bool              mixed_gradients)
{
  if (!frame_lower.intersects(frame_upper)) return;

  const int nx = lower.shape.x;
  const int ny = lower.shape.y;

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const float inv_nx = nx > 1 ? 1.f / float(nx - 1) : 0.f;
  const float inv_ny = ny > 1 ? 1.f / float(ny - 1) : 0.f;

  Array                upper_proj = lower;
  Array                feather_w(lower.shape, 1.f);
  std::vector<uint8_t> in_mask(size_t(nx) * size_t(ny), 0);

  for (int j = 0; j < ny; ++j)
  {
    const float ry_lower = float(j) * inv_ny;

    for (int i = 0; i < nx; ++i)
    {
      const float     rx_lower = float(i) * inv_nx;
      const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

      if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y)) continue;

      in_mask[size_t(j) * size_t(nx) + size_t(i)] = 1;
      upper_proj(i, j) = sample_array_bilinear(upper, r_upper.x, r_upper.y);
      feather_w(i,
                j) = edge_distance_weight(r_upper.x, r_upper.y, feather_width);
    }
  }

  Array                u;
  Array                rhs;
  std::vector<uint8_t> is_interior;
  build_poisson_patch_inputs(lower,
                             upper_proj,
                             in_mask,
                             feather_w,
                             mixed_gradients,
                             u,
                             rhs,
                             is_interior);

  const float pi = static_cast<float>(M_PI);
  const float rho = 0.5f * (std::cos(pi / float(std::max(2, nx))) +
                            std::cos(pi / float(std::max(2, ny))));
  const float omega = 2.f / (1.f + std::sqrt(std::max(0.f, 1.f - rho * rho)));

  solve_poisson_defect_red_black(u, rhs, is_interior, iterations, omega);

  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      if (in_mask[size_t(j) * size_t(nx) + size_t(i)])
        lower(i, j) = upper_proj(i, j) + u(i, j);
}

void blend_terrains_poisson(VirtualArray       &lower,
                            const CoordFrame   &frame_lower,
                            const VirtualArray &upper,
                            const CoordFrame   &frame_upper,
                            int                 iterations,
                            float               feather_width,
                            bool                mixed_gradients,
                            int                 coarse_size,
                            int                 sync_interval,
                            const ComputeMode  &cm)
{
  if (!frame_lower.intersects(frame_upper)) return;

  // Level 1: coarse global Poisson defect solve
  const glm::ivec2 coarse_shape(
      std::min(lower.shape.x, std::max(4, coarse_size)),
      std::min(lower.shape.y, std::max(4, coarse_size)));

  Array lower_coarse = lower.to_array(coarse_shape, ComputeMode{});
  Array upper_coarse = upper.to_array(
      glm::ivec2(std::min(upper.shape.x, std::max(4, coarse_size)),
                 std::min(upper.shape.y, std::max(4, coarse_size))),
      ComputeMode{});

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  Array                upper_proj_c = lower_coarse;
  Array                feather_w_c(coarse_shape, 1.f);
  std::vector<uint8_t> in_mask_c(size_t(coarse_shape.x) *
                                     size_t(coarse_shape.y),
                                 0);

  for (int j = 0; j < coarse_shape.y; ++j)
  {
    const float ry = coarse_shape.y > 1 ? float(j) / float(coarse_shape.y - 1)
                                        : 0.f;
    for (int i = 0; i < coarse_shape.x; ++i)
    {
      const float rx = coarse_shape.x > 1 ? float(i) / float(coarse_shape.x - 1)
                                          : 0.f;
      const glm::vec2 r_u = xform.apply(rx, ry);
      if (!frame_upper.is_unit_point_within(r_u.x, r_u.y)) continue;

      in_mask_c[size_t(j) * size_t(coarse_shape.x) + size_t(i)] = 1;
      upper_proj_c(i, j) = sample_array_bilinear(upper_coarse, r_u.x, r_u.y);
      feather_w_c(i, j) = edge_distance_weight(r_u.x, r_u.y, feather_width);
    }
  }

  Array                u_coarse;
  Array                rhs_c;
  std::vector<uint8_t> is_interior_c;
  build_poisson_patch_inputs(lower_coarse,
                             upper_proj_c,
                             in_mask_c,
                             feather_w_c,
                             mixed_gradients,
                             u_coarse,
                             rhs_c,
                             is_interior_c);

  const float pi = static_cast<float>(M_PI);
  const float rho_c = 0.5f * (std::cos(pi / float(coarse_shape.x)) +
                              std::cos(pi / float(coarse_shape.y)));
  const float omega_c = 2.f /
                        (1.f + std::sqrt(std::max(0.f, 1.f - rho_c * rho_c)));

  solve_poisson_defect_red_black(u_coarse,
                                 rhs_c,
                                 is_interior_c,
                                 std::max(200, iterations),
                                 omega_c);

  // Level 2: tile-local Red-Black SOR refinement with halo sync
  const glm::vec4 upper_bbox = {0.f, 1.f, 0.f, 1.f};
  const float     upper_dx = upper_bbox.y - upper_bbox.x;
  const float     upper_dy = upper_bbox.w - upper_bbox.z;

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;

  for_each_tile(
      lower,
      [&](Array &tile, const UvTile &region, ICommunicator &comm)
      {
        const glm::ivec2 offset = region.global_offset();
        const int        tx = tile.shape.x;
        const int        ty = tile.shape.y;

        Array                upper_proj = tile;
        Array                feather_w(tile.shape, 1.f);
        std::vector<uint8_t> in_mask(size_t(tx) * size_t(ty), 0);

        bool any_inside = false;
        for (int j = 0; j < ty; ++j)
        {
          const float ry_lower = float(offset.y + j) * inv_ny;

          for (int i = 0; i < tx; ++i)
          {
            const float     rx_lower = float(offset.x + i) * inv_nx;
            const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

            if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y))
              continue;

            any_inside = true;
            in_mask[size_t(j) * size_t(tx) + size_t(i)] = 1;
            const float ux = upper_bbox.x + r_upper.x * upper_dx;
            const float uy = upper_bbox.z + r_upper.y * upper_dy;
            upper_proj(i, j) = upper.get_bilinear(ux, uy);
            feather_w(i, j) = edge_distance_weight(r_upper.x,
                                                   r_upper.y,
                                                   feather_width);
          }
        }

        Array                u_tile;
        Array                rhs_tile;
        std::vector<uint8_t> is_interior;
        build_poisson_patch_inputs(tile,
                                   upper_proj,
                                   in_mask,
                                   feather_w,
                                   mixed_gradients,
                                   u_tile,
                                   rhs_tile,
                                   is_interior);

        // initialize interior from coarse global solution
        for (int j = 0; j < ty; ++j)
        {
          const float ry_lower = float(offset.y + j) * inv_ny;
          for (int i = 0; i < tx; ++i)
          {
            const size_t idx = size_t(j) * size_t(tx) + size_t(i);
            if (is_interior[idx])
            {
              const float rx_lower = float(offset.x + i) * inv_nx;
              u_tile(i,
                     j) = sample_array_bilinear(u_coarse, rx_lower, ry_lower);
            }
          }
        }

        const float rho_t = 0.5f * (std::cos(pi / float(std::max(2, tx))) +
                                    std::cos(pi / float(std::max(2, ty))));
        const float omega_t = 2.f /
                              (1.f +
                               std::sqrt(std::max(0.f, 1.f - rho_t * rho_t)));

        const int batch_size = (comm.size() > 1 && sync_interval > 0)
                                   ? sync_interval
                                   : iterations;

        for (int it = 0; it < iterations; it += batch_size)
        {
          const int n_it = std::min(batch_size, iterations - it);
          if (any_inside)
            solve_poisson_defect_red_black(u_tile,
                                           rhs_tile,
                                           is_interior,
                                           n_it,
                                           omega_t);

          if (comm.size() > 1)
            comm.sync_overlap_buffers(u_tile, region, SyncOperation::Average);
        }

        if (any_inside)
        {
          for (int j = 0; j < ty; ++j)
            for (int i = 0; i < tx; ++i)
              if (in_mask[size_t(j) * size_t(tx) + size_t(i)])
                tile(i, j) = upper_proj(i, j) + u_tile(i, j);
        }

        if (comm.size() > 1)
          comm.sync_overlap_buffers(tile, region, SyncOperation::SmoothBlend);
      },
      cm);
}

void blend_terrains_seamless(Array            &lower,
                             const CoordFrame &frame_lower,
                             const Array      &upper,
                             const CoordFrame &frame_upper,
                             float             feather_width,
                             int               boundary_samples)
{
  if (!frame_lower.intersects(frame_upper)) return;

  const auto samples = sample_overlap_boundary(
      frame_lower,
      frame_upper,
      [&](float rx, float ry) { return sample_array_bilinear(lower, rx, ry); },
      [&](float rx, float ry) { return sample_array_bilinear(upper, rx, ry); },
      boundary_samples);

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;

  for (int j = 0; j < lower.shape.y; ++j)
  {
    const float ry_lower = float(j) * inv_ny;

    for (int i = 0; i < lower.shape.x; ++i)
    {
      const float     rx_lower = float(i) * inv_nx;
      const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

      if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y)) continue;

      const float upper_val = sample_array_bilinear(upper,
                                                    r_upper.x,
                                                    r_upper.y);
      const float phi = evaluate_mvc_membrane(samples, r_upper.x, r_upper.y);
      const float w = edge_distance_weight(r_upper.x,
                                           r_upper.y,
                                           feather_width,
                                           true);

      lower(i, j) = lerp(lower(i, j), upper_val + phi, w);
    }
  }
}

void blend_terrains_seamless(VirtualArray       &lower,
                             const CoordFrame   &frame_lower,
                             const VirtualArray &upper,
                             const CoordFrame   &frame_upper,
                             float               feather_width,
                             int                 boundary_samples,
                             const ComputeMode  &cm)
{
  if (!frame_lower.intersects(frame_upper)) return;

  const glm::vec4 lower_bbox = {0.f, 1.f, 0.f, 1.f};
  const float     lower_dx = lower_bbox.y - lower_bbox.x;
  const float     lower_dy = lower_bbox.w - lower_bbox.z;

  const glm::vec4 upper_bbox = {0.f, 1.f, 0.f, 1.f};
  const float     upper_dx = upper_bbox.y - upper_bbox.x;
  const float     upper_dy = upper_bbox.w - upper_bbox.z;

  const auto samples = sample_overlap_boundary(
      frame_lower,
      frame_upper,
      [&](float rx, float ry)
      {
        return lower.get_bilinear(lower_bbox.x + rx * lower_dx,
                                  lower_bbox.z + ry * lower_dy);
      },
      [&](float rx, float ry)
      {
        return upper.get_bilinear(upper_bbox.x + rx * upper_dx,
                                  upper_bbox.z + ry * upper_dy);
      },
      boundary_samples);

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;

  for_each_tile(
      lower,
      [&](Array &tile, const UvTile &region)
      {
        const glm::ivec2 offset = region.global_offset();

        for (int j = 0; j < tile.shape.y; ++j)
        {
          const float ry_lower = float(offset.y + j) * inv_ny;

          for (int i = 0; i < tile.shape.x; ++i)
          {
            const float     rx_lower = float(offset.x + i) * inv_nx;
            const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

            if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y))
              continue;

            const float ux = upper_bbox.x + r_upper.x * upper_dx;
            const float uy = upper_bbox.z + r_upper.y * upper_dy;
            const float upper_val = upper.get_bilinear(ux, uy);

            const float phi = evaluate_mvc_membrane(samples,
                                                    r_upper.x,
                                                    r_upper.y);
            const float w = edge_distance_weight(r_upper.x,
                                                 r_upper.y,
                                                 feather_width,
                                                 true);

            tile(i, j) = lerp(tile(i, j), upper_val + phi, w);
          }
        }
      },
      cm);
}

void blend_terrains_slope(Array            &lower,
                          const CoordFrame &frame_lower,
                          const Array      &upper,
                          const CoordFrame &frame_upper,
                          float             feather_width,
                          float             slope_scale,
                          int               ir_slope)
{
  if (!frame_lower.intersects(frame_upper)) return;

  const int nx = lower.shape.x;
  const int ny = lower.shape.y;

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const float inv_nx = nx > 1 ? 1.f / float(nx - 1) : 0.f;
  const float inv_ny = ny > 1 ? 1.f / float(ny - 1) : 0.f;

  Array                upper_proj = lower;
  std::vector<uint8_t> in_mask(size_t(nx) * size_t(ny), 0);
  Array                r_ux(lower.shape, 0.f);
  Array                r_uy(lower.shape, 0.f);

  for (int j = 0; j < ny; ++j)
  {
    const float ry_lower = float(j) * inv_ny;

    for (int i = 0; i < nx; ++i)
    {
      const float     rx_lower = float(i) * inv_nx;
      const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

      if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y)) continue;

      in_mask[size_t(j) * size_t(nx) + size_t(i)] = 1;
      upper_proj(i, j) = sample_array_bilinear(upper, r_upper.x, r_upper.y);
      r_ux(i, j) = r_upper.x;
      r_uy(i, j) = r_upper.y;
    }
  }

  Array gn_lower = gradient_norm(lower);
  Array gn_upper = gradient_norm(upper_proj);

  if (ir_slope > 0)
  {
    smooth_cpulse(gn_lower, ir_slope);
    smooth_cpulse(gn_upper, ir_slope);
  }

  const float inv_scale = 1.f / std::max(slope_scale, 1e-6f);

  for (int j = 0; j < ny; ++j)
  {
    for (int i = 0; i < nx; ++i)
    {
      if (!in_mask[size_t(j) * size_t(nx) + size_t(i)]) continue;

      const float w_d = edge_distance_weight(r_ux(i, j),
                                             r_uy(i, j),
                                             feather_width);
      const float t_s = std::clamp(
          (gn_upper(i, j) - gn_lower(i, j)) * inv_scale + 0.5f,
          0.f,
          1.f);
      const float w_s = smoothstep3(t_s);

      lower(i, j) = lerp(lower(i, j), upper_proj(i, j), w_d * w_s);
    }
  }
}

void blend_terrains_slope(VirtualArray       &lower,
                          const CoordFrame   &frame_lower,
                          const VirtualArray &upper,
                          const CoordFrame   &frame_upper,
                          float               feather_width,
                          float               slope_scale,
                          int                 ir_slope,
                          const ComputeMode  &cm)
{
  if (!frame_lower.intersects(frame_upper)) return;

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const glm::vec4 upper_bbox = {0.f, 1.f, 0.f, 1.f};
  const float     upper_dx = upper_bbox.y - upper_bbox.x;
  const float     upper_dy = upper_bbox.w - upper_bbox.z;

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;
  const float inv_scale = 1.f / std::max(slope_scale, 1e-6f);

  for_each_tile(
      lower,
      [&](Array &tile, const UvTile &region, ICommunicator &comm)
      {
        const glm::ivec2 offset = region.global_offset();
        const int        tx = tile.shape.x;
        const int        ty = tile.shape.y;

        Array                upper_proj = tile;
        std::vector<uint8_t> in_mask(size_t(tx) * size_t(ty), 0);
        Array                r_ux(tile.shape, 0.f);
        Array                r_uy(tile.shape, 0.f);

        bool any_inside = false;
        for (int j = 0; j < ty; ++j)
        {
          const float ry_lower = float(offset.y + j) * inv_ny;

          for (int i = 0; i < tx; ++i)
          {
            const float     rx_lower = float(offset.x + i) * inv_nx;
            const glm::vec2 r_upper = xform.apply(rx_lower, ry_lower);

            if (!frame_upper.is_unit_point_within(r_upper.x, r_upper.y))
              continue;

            any_inside = true;
            in_mask[size_t(j) * size_t(tx) + size_t(i)] = 1;
            const float ux = upper_bbox.x + r_upper.x * upper_dx;
            const float uy = upper_bbox.z + r_upper.y * upper_dy;
            upper_proj(i, j) = upper.get_bilinear(ux, uy);
            r_ux(i, j) = r_upper.x;
            r_uy(i, j) = r_upper.y;
          }
        }

        if (any_inside)
        {
          Array gn_lower = gradient_norm(tile);
          Array gn_upper = gradient_norm(upper_proj);

          if (ir_slope > 0)
          {
            smooth_cpulse(gn_lower, ir_slope);
            smooth_cpulse(gn_upper, ir_slope);
          }

          for (int j = 0; j < ty; ++j)
          {
            for (int i = 0; i < tx; ++i)
            {
              if (!in_mask[size_t(j) * size_t(tx) + size_t(i)]) continue;

              const float w_d = edge_distance_weight(r_ux(i, j),
                                                     r_uy(i, j),
                                                     feather_width);
              const float t_s = std::clamp(
                  (gn_upper(i, j) - gn_lower(i, j)) * inv_scale + 0.5f,
                  0.f,
                  1.f);
              const float w_s = smoothstep3(t_s);

              tile(i, j) = lerp(tile(i, j), upper_proj(i, j), w_d * w_s);
            }
          }
        }

        if (comm.size() > 1)
          comm.sync_overlap_buffers(tile, region, SyncOperation::SmoothBlend);
      },
      cm);
}

} // namespace hmap
