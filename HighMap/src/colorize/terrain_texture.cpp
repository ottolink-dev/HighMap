/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/colorize.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/math.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"
#include "highmap/texture.hpp"

#include <numbers>

namespace hmap
{

Texture terrain_texture(const Array                &z,
                        const Array                *p_deposition,
                        const TerrainTextureParams &p,
                        const glm::vec4            &bbox,
                        float                       zmin,
                        float                       zmax,
                        const Array                *p_noise)
{
  if (!validate_non_empty(z)) return Texture();
  if (p_noise && !validate_same_shape(z, *p_noise)) return Texture();

  // resolve elevation range
  if (zmin > zmax)
  {
    const glm::vec2 z_range = z.range();
    zmin = z_range.x;
    zmax = z_range.y;
  }

  const glm::ivec2 shape = z.shape;
  const int        nx = shape.x;
  const int        ny = shape.y;

  // pixel radii computed from a single reference scale (region radius)
  const int ir_region = std::max(1, p.ir_region);
  const int ir_crest = std::max(1, int(ir_region / 5.f));
  const int ir_fine = std::max(1, int(ir_region / 20.f));
  const int ir_smooth = std::max(1, int(ir_region / 25.6f));

  const float zrange = std::max(zmax - zmin, 1e-6f);

  // --- Terrain analysis

  // slope and fall lines, read from a lightly smoothed surface
  Array zs = z;
  smooth_cpulse(zs, ir_smooth);

  Array gx, gy;
  Array slope = gradient_norm(zs, &gx, &gy);

  const Array fine = topographic_residual(z, ir_fine, true);
  const Array crest = topographic_residual(z, ir_crest, true);
  const Array region = topographic_residual(z, ir_region, true);

  // drainage on a log scale: 0 on the ridges, 1 in the main river
  Array flow = flow_accumulation_dinf(zs, 0.1f);
  {
    const float lmax = std::log(std::max(flow.max(), 2.f));
    for (float &v : flow.vector)
      v = std::log(std::max(v, 1.f)) / lmax;
  }
  Array moisture = flow;
  smooth_cpulse(moisture, ir_crest);

  // normalize deposition into [0, 1] relative to the 90th percentile
  // of non-zero cells; this prevents extreme outlier deposits from
  // dominating the scale while preserving sparsity where no sediment
  // accumulated
  Array deposit(shape, 0.f);
  if (p_deposition && p_deposition->shape == shape)
  {
    // reference scale: 90th percentile among active deposits
    // (ignoring zero cells)
    const float ref = p_deposition->range_percentile(0.f, 0.9f, 1024, true).y;

    if (ref > 1e-6f)
    {
      deposit = *p_deposition / ref;
      clamp(deposit, 0.f, 1.f);
    }
  }

  // --- Noise

  Array n_large;
  Array n_mid;
  Array n_fine;
  Array n_rill;
  Array n_grain;

  if (p_noise)
  {
    n_large = *p_noise;
    n_mid = *p_noise;
    n_fine = *p_noise;
    n_rill = *p_noise;
    n_grain = *p_noise;
  }
  else if (p.color_noise <= 0.f)
  {
    n_large = Array(shape, 0.f);
    n_mid = Array(shape, 0.f);
    n_fine = Array(shape, 0.f);
    n_rill = Array(shape, 0.f);
    n_grain = Array(shape, 0.f);
  }
  else
  {
    const auto noise = [&](float kw, int seed_offset, int octaves)
    {
      const float k = kw * p.noise_scale;
      Array       a = noise_fbm(NoiseType::PERLIN,
                          shape,
                                {k, k},
                          p.seed + seed_offset,
                          octaves,
                          0.f,
                          0.5f,
                          2.f,
                          nullptr,
                          nullptr,
                          nullptr,
                          bbox);
      normalize_and_clip(a);
      return a;
    };

    n_large = noise(3.f, 0, 4);
    n_mid = noise(12.f, 1, 5);
    n_fine = noise(64.f, 2, 3);
    n_rill = noise(160.f, 3, 2);
    n_grain = noise(320.f, 4, 2);
  }

  // --- Fall lines

  // bare high crests: their colour runs down the slopes below them
  Array source(shape);

#pragma omp parallel for
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      const float h = (z(i, j) - zmin) / zrange;
      source(i, j) = smoothstep3(0.6f * crest(i, j) + 0.6f * region(i, j),
                                 0.05f,
                                 0.35f) *
                     smoothstep3(h + 0.15f * n_mid(i, j), 0.3f, 0.75f) *
                     (0.75f + 0.25f * n_fine(i, j));
    }

  const int streak_steps = std::max(0, p.streak_length);
  const int rill_steps = std::max(1, int(ir_region / 5.f));

  Array streak(shape), rill(shape);

#pragma omp parallel for
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      // walk up the fall line, one cell per step when step is 1
      const auto walk = [&](int steps, float step, auto &&visit)
      {
        float x = (float)i, y = (float)j;
        for (int s = 1; s <= steps; ++s)
        {
          const int   cx = int(x + 0.5f);
          const int   cy = int(y + 0.5f);
          const float g = std::hypot(gx(cx, cy), gy(cx, cy));
          if (g < 1e-6f) return;

          x += step * gx(cx, cy) / g;
          y += step * gy(cx, cy) / g;
          if (x < 0.f || y < 0.f || x > nx - 1.f || y > ny - 1.f) return;

          visit(s, int(x + 0.5f), int(y + 0.5f));
        }
      };

      // fine noise averaged along the fall line: rills and gullies
      float sum = n_rill(i, j);
      int   count = 1;
      walk(rill_steps,
           1.f,
           [&](int, int ci, int cj)
           {
             sum += n_rill(ci, cj);
             ++count;
           });
      rill(i, j) = std::clamp(sum / std::sqrt((float)count), -1.f, 1.f);

      float best = source(i, j);
      walk(streak_steps,
           1.f,
           [&](int s, int ci, int cj)
           {
             best = std::max(best,
                             source(ci, cj) * (1.f - s / (streak_steps + 1.f)));
           });
      streak(i, j) = best;
    }

  // --- Colors

  Texture out(shape, 4);

  const float cn = p.color_noise;
  const float two_pi = 2.f * std::numbers::pi_v<float>;
  const float flow_width = 0.8f - 0.5f * p.flow_width;

#pragma omp parallel for
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      const float h = (z(i, j) - zmin) / zrange;
      const float s = slope(i, j) / p.rock_talus; // 1 where rock takes over
      const float c = crest(i, j);
      const float r = region(i, j);
      const float nl = n_large(i, j);
      const float nm = n_mid(i, j);
      const float nf = n_fine(i, j);
      const float ng = n_grain(i, j);
      const float rl = rill(i, j);

      // --- Grass & dry grass

      const float wet_val = 0.5f + 0.6f * moisture(i, j) - 0.35f * r -
                            0.2f * c + 0.15f * nl;
      const float dry_val = 1.f - wet_val + 0.35f * (h - 0.4f) + 0.12f * nm +
                            0.08f * nf;
      const float dry_mask = smoothstep3(dry_val, 0.42f, 0.58f);

      const float     grass_sh = cn * (0.3f * nm + 0.3f * nf + 0.2f * ng);
      const glm::vec3 c_grass = shade_color(p.grass, grass_sh);
      const glm::vec3 c_dry = shade_color(p.dry, grass_sh);

      glm::vec3 col = glm::mix(c_grass, c_dry, dry_mask);

      // --- Soil

      const float soil_val = s + 0.2f * nm + 0.1f * nf +
                             0.25f * std::max(c, 0.f);
      const float soil_mask = smoothstep3(soil_val, 0.35f, 0.6f);
      const float soil_sh = cn * (0.25f * rl + 0.3f * nm + 0.15f * ng);

      col = glm::mix(col, shade_color(p.soil, soil_sh), soil_mask);

      // --- Sediment

      const float flow_mask = smoothstep3(flow(i, j) + 0.05f * nf,
                                          flow_width,
                                          flow_width + 0.08f);
      const float dep_mask = smoothstep3(deposit(i, j) + 0.15f * nf,
                                         0.4f,
                                         0.8f);

      float sed_mask = std::max(p.flow * flow_mask, p.deposition * dep_mask);
      sed_mask *= 1.f - smoothstep3(s, 0.6f, 1.f);

      const float sed_sh = cn * (0.4f * nf + 0.2f * nm + 0.2f * ng);
      col = glm::mix(col, shade_color(p.sediment, sed_sh), sed_mask);

      // --- Rock & strata

      const float edge = p.rock_edge;
      const float rock_val = s + 0.25f * nm + 0.1f * nf +
                             0.3f * std::max(c, 0.f) * h;
      const float rock_mask = smoothstep3(rock_val, 1.f - edge, 1.f + edge);

      const float band = std::sin(two_pi * (14.f * h + 0.4f * nl + 0.1f * nm));
      const float rock_sh = 0.35f * p.strata * band +
                            cn * (0.2f * rl + 0.3f * nf + 0.15f * ng) +
                            0.25f * c;

      col = glm::mix(col, shade_color(p.rock, rock_sh), rock_mask);

      // --- Ridge streaks

      const float streak_val = streak(i, j) * (0.75f + 0.25f * (1.f + nf));
      const float streak_mask = smoothstep3(streak_val, 0.25f, 0.65f);
      const float slope_gate = smoothstep3(s, 0.15f, 0.4f);
      const float spill_mask = p.streaks * streak_mask * slope_gate;

      const float ridge_sh = cn * (0.3f * nf + 0.15f * rl + 0.15f * ng);
      col = glm::mix(col, shade_color(p.ridge, ridge_sh), spill_mask);

      // --- Snow

      const float snow_val = h - p.snow_line + 0.08f * nm + 0.04f * nf -
                             0.15f * std::max(0.f, s - 0.9f) - 0.05f * c;
      const float snow_mask = p.snow_amount *
                              smoothstep3(snow_val, -0.04f, 0.04f);
      const float snow_sh = 0.15f * nf - 0.3f * std::max(-c, 0.f);

      col = glm::mix(col, shade_color(p.snow, snow_sh), snow_mask);

      // --- Micro-relief & ambient shading

      const float cavity_factor = std::clamp(-0.7f * c - 0.4f * r, 0.f, 1.f);
      col *= 1.f - p.cavity * 0.35f * cavity_factor;
      col *= 1.f + p.detail * 0.22f * fine(i, j);
      col *= 1.f + cn * (0.05f * nm + 0.06f * nf + 0.05f * ng);
      col += cn * 0.015f * nl * glm::vec3(1.f, 0.6f, 0.1f);

      // --- Output RGBA

      out(i, j, 0) = std::clamp(col.r, 0.f, 1.f);
      out(i, j, 1) = std::clamp(col.g, 0.f, 1.f);
      out(i, j, 2) = std::clamp(col.b, 0.f, 1.f);
      out(i, j, 3) = 1.f;
    }

  return out;
}

} // namespace hmap
