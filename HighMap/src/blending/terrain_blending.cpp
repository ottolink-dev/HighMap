/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

#include <algorithm>

#include "highmap/blending.hpp"
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
  }

  return lower_val;
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

} // namespace

// --- Public API

void blend_terrains(Array            &lower,
                    const CoordFrame &frame_lower,
                    const Array      &upper,
                    const CoordFrame &frame_upper,
                    TerrainBlendMode  mode,
                    float             param)
{
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
  if (!frame_lower.intersects(frame_upper)) return;

  const AffineTransform2D xform = frame_lower.compute_relative_transform(
      frame_upper);

  const glm::vec4 upper_bbox = upper.bbox;
  const float     upper_dx = upper_bbox.y - upper_bbox.x;
  const float     upper_dy = upper_bbox.w - upper_bbox.z;

  const float inv_nx = lower.shape.x > 1 ? 1.f / float(lower.shape.x - 1) : 0.f;
  const float inv_ny = lower.shape.y > 1 ? 1.f / float(lower.shape.y - 1) : 0.f;

  for_each_tile(
      lower,
      [&](Array &tile, const TileRegion &region)
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

} // namespace hmap
