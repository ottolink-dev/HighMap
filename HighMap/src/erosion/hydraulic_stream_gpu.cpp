/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/blending.hpp"
#include "highmap/erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/math/array.hpp"
#include "highmap/random.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

namespace
{

void bspline_weights(float t, float w[4])
{
  float t2 = t * t, t3 = t2 * t;
  w[0] = (1.f - 3.f * t + 3.f * t2 - t3) / 6.f;
  w[1] = (4.f - 6.f * t2 + 3.f * t3) / 6.f;
  w[2] = (1.f + 3.f * t + 3.f * t2 - 3.f * t3) / 6.f;
  w[3] = t3 / 6.f;
}

Array upsample_bicubic_warp(const Array  &src,
                            glm::ivec2    dst_shape,
                            float         warp,
                            std::uint32_t seed)
{
  if (warp <= 0.f) return src.resample_to_shape_bicubic(dst_shape);

  int sw = src.shape.x;
  int sh = src.shape.y;
  int dw = dst_shape.x;
  int dh = dst_shape.y;

  Array dst(dst_shape);

  const float rx = static_cast<float>(sw) / static_cast<float>(dw);
  const float ry = static_cast<float>(sh) / static_cast<float>(dh);

#pragma omp parallel for collapse(2) schedule(static)
  for (int y = 0; y < dh; ++y)
  {
    for (int x = 0; x < dw; ++x)
    {
      float fx = (static_cast<float>(x) + 0.5f) * rx - 0.5f;
      float fy = (static_cast<float>(y) + 0.5f) * ry - 0.5f;

      if (warp > 0.f)
      {
        float u = fx * 0.5f;
        float v = fy * 0.5f;
        fx += warp * vnoise(u + 13.7f, v + 3.1f, seed);
        fy += warp * vnoise(u + 71.3f, v + 47.9f, seed + 101u);
      }

      int   xi = static_cast<int>(std::floor(fx));
      int   yi = static_cast<int>(std::floor(fy));
      float wx[4], wy[4];
      bspline_weights(fx - static_cast<float>(xi), wx);
      bspline_weights(fy - static_cast<float>(yi), wy);

      float val = 0.f;
      for (int j = 0; j < 4; ++j)
      {
        int   row = std::min(sh - 1, std::max(0, yi - 1 + j));
        float r = 0.f;
        for (int i = 0; i < 4; ++i)
        {
          int col = std::min(sw - 1, std::max(0, xi - 1 + i));
          r += wx[i] * src(col, row);
        }
        val += wy[j] * r;
      }
      dst(x, y) = val;
    }
  }

  return dst;
}

} // namespace

void hydraulic_stream_log(Array &z,
                          float  c_erosion,
                          float  talus_ref,
                          int    deposition_ir,
                          float  deposition_scale_ratio,
                          float  gradient_power,
                          float  gradient_scaling_ratio,
                          int    gradient_prefilter_ir,
                          float  saturation_ratio,
                          Array *p_bedrock,
                          Array *p_moisture_map,
                          Array *p_erosion_map,
                          Array *p_deposition_map,
                          Array *p_flow_map)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;

  // keep a backup of the input if the erosion / deposition maps need
  // to be computed
  Array z_bckp = Array();
  if (p_erosion_map || p_deposition_map) z_bckp = z;

  // use flow accumulation to determine erosion intensity
  Array facc = flow_accumulation_dinf(z, talus_ref);
  facc = log10(facc);
  remap(facc);

  if (saturation_ratio < 1.f)
    saturate(facc, 0.f, saturation_ratio, 0.1f * saturation_ratio);

  // scale erosion with local gradient
  Array gn = hmap::gradient_norm(z);

  // TODO next operations in one local pass
  {
    gpu::smooth_cpulse(gn, gradient_prefilter_ir);
    remap(gn);
    gn = pow(gn, gradient_power);
    gn = smoothstep5_lower(gn);
    facc *= (1.f - gradient_scaling_ratio) + gradient_scaling_ratio * gn;
  }

  // preserve local peaks
  {
    Array re = gpu::relative_elevation_square_kernel(z, gradient_prefilter_ir);
    re = smoothstep3(1.f - re);
    facc *= re;
  }

  if (p_moisture_map)
    z -= (*p_moisture_map) * c_erosion * facc;
  else
    z -= c_erosion * facc;

  // mimic deposition
  Array zd = z;
  gpu::smooth_fill_holes(zd, deposition_ir);
  zd = gpu::blend_gradients(zd, z, deposition_ir);
  z = lerp(z, zd, deposition_scale_ratio);

  // enforce bedrock
  if (p_bedrock) z = maximum(*p_bedrock, z);

  // splatmaps
  if (p_erosion_map)
  {
    *p_erosion_map = z_bckp - z;
    clamp_min(*p_erosion_map, 0.f);
  }

  if (p_deposition_map)
  {
    *p_deposition_map = z - z_bckp;
    clamp_min(*p_deposition_map, 0.f);
  }

  if (p_flow_map) *p_flow_map = facc;
}

void hydraulic_stream_log(Array       &z,
                          float        c_erosion,
                          float        talus_ref,
                          const Array *p_mask,
                          int          deposition_ir,
                          float        deposition_scale_ratio,
                          float        gradient_power,
                          float        gradient_scaling_ratio,
                          int          gradient_prefilter_ir,
                          float        saturation_ratio,
                          Array       *p_bedrock,
                          Array       *p_moisture_map,
                          Array       *p_erosion_map,
                          Array       *p_deposition_map,
                          Array       *p_flow_map)
{
  if (!validate_non_empty(z)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;

  if (!p_mask)
    gpu::hydraulic_stream_log(z,
                              c_erosion,
                              talus_ref,
                              deposition_ir,
                              deposition_scale_ratio,
                              gradient_power,
                              gradient_scaling_ratio,
                              gradient_prefilter_ir,
                              saturation_ratio,
                              p_bedrock,
                              p_moisture_map,
                              p_erosion_map,
                              p_deposition_map,
                              p_flow_map);
  else
  {
    Array z_f = z;
    gpu::hydraulic_stream_log(z_f,
                              c_erosion,
                              talus_ref,
                              deposition_ir,
                              deposition_scale_ratio,
                              gradient_power,
                              gradient_scaling_ratio,
                              gradient_prefilter_ir,
                              saturation_ratio,
                              p_bedrock,
                              p_moisture_map,
                              p_erosion_map,
                              p_deposition_map,
                              p_flow_map);
    z = lerp(z, z_f, *(p_mask));
  }
}

void hydraulic_stream_log_multiscale(Array        &z,
                                     float         c_erosion,
                                     float         talus_ref,
                                     int           levels,
                                     int           deposition_ir,
                                     float         deposition_scale_ratio,
                                     float         gradient_power,
                                     float         gradient_scaling_ratio,
                                     int           gradient_prefilter_ir,
                                     float         saturation_ratio,
                                     const Array  *p_bedrock,
                                     const Array  *p_moisture_map,
                                     Array        *p_erosion_map,
                                     Array        *p_deposition_map,
                                     Array        *p_flow_map,
                                     float         mix,
                                     float         warp,
                                     std::uint32_t seed)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;

  int nlevels = std::max(1, levels);

  Array z_bckp = Array();
  if ((p_erosion_map != nullptr) || (p_deposition_map != nullptr)) z_bckp = z;

  // build halving pyramid ladder (coarsest first, final level == z.shape)
  std::vector<glm::ivec2> ladder(nlevels);
  for (int i = 0; i < nlevels; ++i)
  {
    int shift = nlevels - 1 - i;
    ladder[i] = {std::max(2, z.shape.x >> shift),
                 std::max(2, z.shape.y >> shift)};
  }

  // iterate coarse-to-fine across resolution levels
  for (int i = 0; i < nlevels; ++i)
  {
    if (i < nlevels - 1)
    {
      // --- Coarse resolution level: compute macro-scale valley incision delta

      Array z_coarse = z.resample_to_shape(ladder[i]);
      Array z_coarse_before = z_coarse;

      // resample optional input maps to coarse shape
      Array  level_bedrock, level_moisture;
      Array *p_lvl_bedrock = nullptr;
      Array *p_lvl_moisture = nullptr;

      if (p_bedrock)
      {
        level_bedrock = p_bedrock->resample_to_shape(ladder[i]);
        p_lvl_bedrock = &level_bedrock;
      }
      if (p_moisture_map)
      {
        level_moisture = p_moisture_map->resample_to_shape(ladder[i]);
        p_lvl_moisture = &level_moisture;
      }

      int level_deposition_ir = std::max(1,
                                         (deposition_ir * ladder[i].x) /
                                             z.shape.x);
      int level_gradient_prefilter_ir = std::max(
          1,
          (gradient_prefilter_ir * ladder[i].x) / z.shape.x);
      float level_talus_ref = talus_ref * static_cast<float>(z.shape.x) /
                              static_cast<float>(ladder[i].x);

      // run erosion on coarse grid without deposition
      gpu::hydraulic_stream_log(z_coarse,
                                c_erosion,
                                level_talus_ref,
                                level_deposition_ir,
                                /* deposition_scale_ratio */ 0.f,
                                gradient_power,
                                gradient_scaling_ratio,
                                level_gradient_prefilter_ir,
                                saturation_ratio,
                                p_lvl_bedrock,
                                p_lvl_moisture,
                                nullptr,
                                nullptr,
                                nullptr);

      // compute incision delta and upsample smoothly to full resolution with
      // optional domain warp
      Array delta_coarse = z_coarse_before - z_coarse;
      clamp_min(delta_coarse, 0.f);

      Array delta_full = upsample_bicubic_warp(
          delta_coarse,
          z.shape,
          warp,
          seed + static_cast<std::uint32_t>(i * 37));
      clamp_min(delta_full, 0.f);

      // carve coarse valley delta into full-resolution terrain
      z -= std::clamp(mix, 0.f, 1.f) * delta_full;
      if (p_bedrock) z = maximum(*p_bedrock, z);
    }
    else
    {
      // --- Final level: run erosion and deposition directly on full-detail
      // terrain

      Array  bedrock_copy, moisture_copy;
      Array *p_b_arg = nullptr;
      Array *p_m_arg = nullptr;
      if (p_bedrock)
      {
        bedrock_copy = *p_bedrock;
        p_b_arg = &bedrock_copy;
      }
      if (p_moisture_map)
      {
        moisture_copy = *p_moisture_map;
        p_m_arg = &moisture_copy;
      }

      gpu::hydraulic_stream_log(z,
                                c_erosion,
                                talus_ref,
                                deposition_ir,
                                deposition_scale_ratio,
                                gradient_power,
                                gradient_scaling_ratio,
                                gradient_prefilter_ir,
                                saturation_ratio,
                                p_b_arg,
                                p_m_arg,
                                nullptr,
                                nullptr,
                                p_flow_map);
    }
  }

  // splatmaps
  if (p_erosion_map)
  {
    *p_erosion_map = z_bckp - z;
    clamp_min(*p_erosion_map, 0.f);
  }

  if (p_deposition_map)
  {
    *p_deposition_map = z - z_bckp;
    clamp_min(*p_deposition_map, 0.f);
  }
}

void hydraulic_stream_log_multiscale(Array        &z,
                                     float         c_erosion,
                                     float         talus_ref,
                                     const Array  *p_mask,
                                     int           levels,
                                     int           deposition_ir,
                                     float         deposition_scale_ratio,
                                     float         gradient_power,
                                     float         gradient_scaling_ratio,
                                     int           gradient_prefilter_ir,
                                     float         saturation_ratio,
                                     const Array  *p_bedrock,
                                     const Array  *p_moisture_map,
                                     Array        *p_erosion_map,
                                     Array        *p_deposition_map,
                                     Array        *p_flow_map,
                                     float         mix,
                                     float         warp,
                                     std::uint32_t seed)
{
  if (!validate_non_empty(z)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;

  if (!p_mask)
    gpu::hydraulic_stream_log_multiscale(z,
                                         c_erosion,
                                         talus_ref,
                                         levels,
                                         deposition_ir,
                                         deposition_scale_ratio,
                                         gradient_power,
                                         gradient_scaling_ratio,
                                         gradient_prefilter_ir,
                                         saturation_ratio,
                                         p_bedrock,
                                         p_moisture_map,
                                         p_erosion_map,
                                         p_deposition_map,
                                         p_flow_map,
                                         mix,
                                         warp,
                                         seed);
  else
  {
    Array z_f = z;
    gpu::hydraulic_stream_log_multiscale(z_f,
                                         c_erosion,
                                         talus_ref,
                                         levels,
                                         deposition_ir,
                                         deposition_scale_ratio,
                                         gradient_power,
                                         gradient_scaling_ratio,
                                         gradient_prefilter_ir,
                                         saturation_ratio,
                                         p_bedrock,
                                         p_moisture_map,
                                         p_erosion_map,
                                         p_deposition_map,
                                         p_flow_map,
                                         mix,
                                         warp,
                                         seed);
    z = lerp(z, z_f, *(p_mask));
  }
}

} // namespace hmap::gpu
