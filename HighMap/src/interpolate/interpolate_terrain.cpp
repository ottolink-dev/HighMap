/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "highmap/coord_frame.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/interpolate/interpolate_terrain.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/core.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

namespace
{

struct SourcePatch
{
  Array      data;
  glm::ivec2 origin{0, 0};
  bool       valid = false;
};

SourcePatch assemble_source_patch(const VirtualArray      &h_source,
                                  const glm::ivec2        &target_shape,
                                  const TileRegion        &target_region,
                                  const AffineTransform2D &tf)
{
  SourcePatch patch;

  const float inv_ntx = (target_shape.x > 1)
                            ? (1.f / static_cast<float>(target_shape.x - 1))
                            : 0.f;
  const float inv_nty = (target_shape.y > 1)
                            ? (1.f / static_cast<float>(target_shape.y - 1))
                            : 0.f;

  const glm::ivec2 ij0 = target_region.global_offset();

  const int it_corners[2] = {ij0.x, ij0.x + target_region.shape.x - 1};
  const int jt_corners[2] = {ij0.y, ij0.y + target_region.shape.y - 1};

  float min_xs = std::numeric_limits<float>::max();
  float max_xs = std::numeric_limits<float>::lowest();
  float min_ys = std::numeric_limits<float>::max();
  float max_ys = std::numeric_limits<float>::lowest();

  for (int cj = 0; cj < 2; ++cj)
    for (int ci = 0; ci < 2; ++ci)
    {
      const float     rx = static_cast<float>(it_corners[ci]) * inv_ntx;
      const float     ry = static_cast<float>(jt_corners[cj]) * inv_nty;
      const glm::vec2 uv = tf.apply(rx, ry);

      const float xs = uv.x * static_cast<float>(h_source.shape.x - 1);
      const float ys = uv.y * static_cast<float>(h_source.shape.y - 1);

      min_xs = std::min(min_xs, xs);
      max_xs = std::max(max_xs, xs);
      min_ys = std::min(min_ys, ys);
      max_ys = std::max(max_ys, ys);
    }

  // pad by 1 cell for bilinear / rounding neighborhood
  int is_min = static_cast<int>(std::floor(min_xs)) - 1;
  int is_max = static_cast<int>(std::ceil(max_xs)) + 1;
  int js_min = static_cast<int>(std::floor(min_ys)) - 1;
  int js_max = static_cast<int>(std::ceil(max_ys)) + 1;

  if (is_max < 0 || is_min >= h_source.shape.x || js_max < 0 ||
      js_min >= h_source.shape.y)
    return patch;

  is_min = std::clamp(is_min, 0, h_source.shape.x - 1);
  is_max = std::clamp(is_max, 0, h_source.shape.x - 1);
  js_min = std::clamp(js_min, 0, h_source.shape.y - 1);
  js_max = std::clamp(js_max, 0, h_source.shape.y - 1);

  const glm::ivec2 patch_shape(is_max - is_min + 1, js_max - js_min + 1);
  patch.data = Array(patch_shape);
  patch.origin = {is_min, js_min};
  patch.valid = true;

  const int stx_min = is_min / h_source.tile_shape.x;
  const int stx_max = is_max / h_source.tile_shape.x;
  const int sty_min = js_min / h_source.tile_shape.y;
  const int sty_max = js_max / h_source.tile_shape.y;

  for (int sty = sty_min; sty <= sty_max; ++sty)
    for (int stx = stx_min; stx <= stx_max; ++stx)
    {
      const TileRegion src_reg = h_source.tile_region_from_tile_coords(stx,
                                                                       sty);
      const Array     &src_tile = h_source.storage->get_tile(src_reg);
      const glm::ivec2 src_off = src_reg.global_offset();

      // copy only core cells of the source tile that fall inside patch
      const int core_gx0 = stx * h_source.tile_shape.x;
      const int core_gx1 = std::min(core_gx0 + h_source.tile_shape.x,
                                    h_source.shape.x);
      const int core_gy0 = sty * h_source.tile_shape.y;
      const int core_gy1 = std::min(core_gy0 + h_source.tile_shape.y,
                                    h_source.shape.y);

      const int gx_start = std::max(is_min, core_gx0);
      const int gx_end = std::min(is_max + 1, core_gx1);
      const int gy_start = std::max(js_min, core_gy0);
      const int gy_end = std::min(js_max + 1, core_gy1);

      for (int gy = gy_start; gy < gy_end; ++gy)
      {
        const int ly = gy - src_off.y;
        const int py = gy - js_min;
        for (int gx = gx_start; gx < gx_end; ++gx)
        {
          patch.data(gx - is_min, py) = src_tile(gx - src_off.x, ly);
        }
      }

      h_source.storage->release_tile(src_reg);
    }

  return patch;
}

inline float sample_patch(const SourcePatch &patch,
                          const glm::ivec2  &src_shape,
                          float              us,
                          float              vs,
                          bool               bilinear)
{
  const float fx = us * static_cast<float>(src_shape.x - 1) -
                   static_cast<float>(patch.origin.x);
  const float fy = vs * static_cast<float>(src_shape.y - 1) -
                   static_cast<float>(patch.origin.y);

  if (bilinear)
  {
    int i0 = static_cast<int>(std::floor(fx));
    int j0 = static_cast<int>(std::floor(fy));
    i0 = std::clamp(i0, 0, std::max(0, patch.data.shape.x - 2));
    j0 = std::clamp(j0, 0, std::max(0, patch.data.shape.y - 2));

    const float u = std::clamp(fx - static_cast<float>(i0), 0.f, 1.f);
    const float v = std::clamp(fy - static_cast<float>(j0), 0.f, 1.f);

    return patch.data.get_value_bilinear_at(i0, j0, u, v);
  }
  else
  {
    int i = std::clamp(static_cast<int>(std::round(fx)),
                       0,
                       patch.data.shape.x - 1);
    int j = std::clamp(static_cast<int>(std::round(fy)),
                       0,
                       patch.data.shape.y - 1);
    return patch.data(i, j);
  }
}

} // namespace

void flatten_terrain(VirtualArray       &h_source1,
                     const VirtualArray &h_source2,
                     const CoordFrame   &t_source1,
                     const CoordFrame   &t_source2,
                     const ComputeMode  &cm,
                     bool                bilinear)
{
  if (h_source1.empty() || h_source2.empty()) return;
  if (!t_source1.intersects(t_source2)) return;

  const AffineTransform2D tf = t_source1.compute_relative_transform(t_source2);
  const float             inv_nx = (h_source1.shape.x > 1)
                                       ? (1.f / static_cast<float>(h_source1.shape.x - 1))
                                       : 0.f;
  const float             inv_ny = (h_source1.shape.y > 1)
                                       ? (1.f / static_cast<float>(h_source1.shape.y - 1))
                                       : 0.f;

  auto lambda = [&](Array &tile, const TileRegion &region)
  {
    SourcePatch patch = assemble_source_patch(h_source2,
                                              h_source1.shape,
                                              region,
                                              tf);
    if (!patch.valid) return;

    const glm::ivec2 ij0 = region.global_offset();
    const float      dus = tf.ax * inv_nx;
    const float      dvs = tf.ay * inv_nx;

    for (int j = 0; j < region.shape.y; ++j)
    {
      const float rx0 = static_cast<float>(ij0.x) * inv_nx;
      const float ry = static_cast<float>(ij0.y + j) * inv_ny;
      float       us = tf.ax * rx0 + tf.bx * ry + tf.cx;
      float       vs = tf.ay * rx0 + tf.by * ry + tf.cy;

      for (int i = 0; i < region.shape.x; ++i)
      {
        if (t_source2.is_unit_point_within(us, vs))
        {
          const float v_source2 = sample_patch(patch,
                                               h_source2.shape,
                                               us,
                                               vs,
                                               bilinear);
          const float r = t_source2.normalized_shape_factor_rel(us, vs);
          tile(i, j) = lerp(tile(i, j), v_source2, r);
        }

        us += dus;
        vs += dvs;
      }
    }
  };

  for_each_tile(h_source1, lambda, cm);
}

void flatten_terrain(const std::vector<const VirtualArray *> &h_sources,
                     VirtualArray                            &h_target,
                     const std::vector<const CoordFrame *>   &t_sources,
                     const CoordFrame                        &t_target,
                     const ComputeMode                       &cm,
                     bool                                     bilinear)
{
  if (!validate_non_empty(h_sources, "h_sources") ||
      !validate_non_empty(t_sources, "t_sources"))
    return;

  if (h_sources.size() != t_sources.size())
  {
    hmap::log::warn("Size mismatch between h_sources and t_sources");
    return;
  }

  // interpolation only for the first layer
  hmap::interpolate_terrain(*h_sources[0],
                            h_target,
                            *t_sources[0],
                            t_target,
                            cm,
                            bilinear);

  // process in-place the remaining layers by flattening them on the
  // current state of the target layer: target <= target & source
  for (size_t k = 1; k < h_sources.size(); ++k)
    flatten_terrain(h_target,
                    *h_sources[k],
                    t_target,
                    *t_sources[k],
                    cm,
                    bilinear);
}

void interpolate_terrain(const VirtualArray &h_source,
                         VirtualArray       &h_target,
                         const CoordFrame   &t_source,
                         const CoordFrame   &t_target,
                         const ComputeMode  &cm,
                         bool                bilinear,
                         float               fill_value)
{
  if (h_source.empty() || h_target.empty()) return;

  if (!t_target.intersects(t_source))
  {
    h_target.fill(fill_value, cm);
    return;
  }

  const AffineTransform2D tf = t_target.compute_relative_transform(t_source);
  const float             inv_nx = (h_target.shape.x > 1)
                                       ? (1.f / static_cast<float>(h_target.shape.x - 1))
                                       : 0.f;
  const float             inv_ny = (h_target.shape.y > 1)
                                       ? (1.f / static_cast<float>(h_target.shape.y - 1))
                                       : 0.f;

  auto lambda = [&](Array &tile, const TileRegion &region)
  {
    SourcePatch patch = assemble_source_patch(h_source,
                                              h_target.shape,
                                              region,
                                              tf);
    if (!patch.valid)
    {
      tile = fill_value;
      return;
    }

    const glm::ivec2 ij0 = region.global_offset();
    const float      dus = tf.ax * inv_nx;
    const float      dvs = tf.ay * inv_nx;

    for (int j = 0; j < region.shape.y; ++j)
    {
      const float rx0 = static_cast<float>(ij0.x) * inv_nx;
      const float ry = static_cast<float>(ij0.y + j) * inv_ny;
      float       us = tf.ax * rx0 + tf.bx * ry + tf.cx;
      float       vs = tf.ay * rx0 + tf.by * ry + tf.cy;

      for (int i = 0; i < region.shape.x; ++i)
      {
        if (t_source.is_unit_point_within(us, vs))
          tile(i, j) = sample_patch(patch, h_source.shape, us, vs, bilinear);
        else
          tile(i, j) = fill_value;

        us += dus;
        vs += dvs;
      }
    }
  };

  for_each_tile(h_target, lambda, cm);
}

} // namespace hmap
