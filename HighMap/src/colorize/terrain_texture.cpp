/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "cl_wrapper/run.hpp"

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

namespace hmap
{

Texture terrain_texture(const Array                &z,
                        const Array                *p_deposition,
                        const Array                *p_snow_depth,
                        const TerrainTextureParams &p,
                        const glm::vec4            &bbox,
                        const Array                *p_noise,
                        ICommunicator              &comm)
{
  if (!validate_non_empty(z)) return Texture();
  if (p_noise && !validate_same_shape(z, *p_noise)) return Texture();
  if (p_snow_depth && !validate_same_shape(z, *p_snow_depth)) return Texture();

  // resolve elevation range
  const glm::vec2 z_range = z.range(comm);
  const float     zmin = z_range.x;
  const float     zmax = z_range.y;

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
  gpu::smooth_cpulse(zs, ir_smooth);

  // gradient evaluated on CPU
  Array gx, gy;
  Array slope = gradient_norm(zs, &gx, &gy);

  // precalculate normalized directions on CPU to avoid hypot and division on
  // GPU
  Array dir_x(shape, 0.f);
  Array dir_y(shape, 0.f);

#pragma omp parallel for
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      const float g = std::sqrt(gx(i, j) * gx(i, j) + gy(i, j) * gy(i, j) +
                                1e-9f);
      dir_x(i, j) = gx(i, j) / g;
      dir_y(i, j) = gy(i, j) / g;
    }

  // topographic residuals via GPU
  const Array fine = gpu::topographic_residual(z, ir_fine, true);
  const Array crest = gpu::topographic_residual(z, ir_crest, true);
  const Array region = gpu::topographic_residual(z, ir_region, true);

  // drainage on a log scale: 0 on the ridges, 1 in the main river
  Array flow = flow_accumulation_dinf(zs, 0.1f, 0, {}, comm);
  {
    const float lmax = std::log(std::max(flow.max(comm), 2.f));
    for (float &v : flow.vector)
      v = std::log(std::max(v, 1.f)) / lmax;
  }
  Array moisture = flow;
  gpu::smooth_cpulse(moisture, ir_crest);

  // normalize deposition into [0, 1] relative to the maximum
  Array deposit(shape, 0.f);
  if (p_deposition && p_deposition->shape == shape)
  {
    const float ref = p_deposition->max(comm);
    if (ref > 1e-6f)
    {
      deposit = *p_deposition / ref;
      clamp(deposit, 0.f, 1.f);
    }
  }

  // normalize snow depth into [0, 1] relative to global maximum across tiles
  const bool has_snow_depth = (p_snow_depth != nullptr);
  Array      snow_depth(shape, 0.f);
  if (has_snow_depth)
  {
    const float max_snow_depth = p_snow_depth->max(comm);
    if (max_snow_depth > 1e-6f)
    {
      snow_depth = *p_snow_depth / max_snow_depth;
      clamp(snow_depth, 0.f, 1.f);
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
      Array       a = gpu::noise_fbm(NoiseType::PERLIN,
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

  // --- Bare high crests source layer

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

  // --- Fall-line walk (Kernel 1)

  const int streak_steps = std::max(0, p.streak_length);
  const int rill_steps = std::max(1, int(ir_region / 5.f));

  Array streak(shape);
  Array rill(shape);

  {
    auto run_walk = clwrapper::Run("terrain_fallline_walk");

    run_walk.bind_imagef("dir_x", dir_x.vector, nx, ny);
    run_walk.bind_imagef("dir_y", dir_y.vector, nx, ny);
    run_walk.bind_imagef("in_source", source.vector, nx, ny);
    run_walk.bind_imagef("in_n_rill", n_rill.vector, nx, ny);
    run_walk.bind_imagef("out_rill", rill.vector, nx, ny, true);
    run_walk.bind_imagef("out_streak", streak.vector, nx, ny, true);

    run_walk.bind_arguments(nx, ny, rill_steps, streak_steps);
    run_walk.execute({nx, ny});

    run_walk.read_imagef("out_rill");
    run_walk.read_imagef("out_streak");
  }

  // --- Shading & color compositing (Kernel 2)

  Texture out(shape, 4);

  {
    auto run_shade = clwrapper::Run("terrain_shade");

    run_shade.bind_imagef("in_z", z.vector, nx, ny);
    run_shade.bind_imagef("in_slope", slope.vector, nx, ny);
    run_shade.bind_imagef("in_crest", crest.vector, nx, ny);
    run_shade.bind_imagef("in_region", region.vector, nx, ny);
    run_shade.bind_imagef("in_fine", fine.vector, nx, ny);
    run_shade.bind_imagef("in_flow", flow.vector, nx, ny);
    run_shade.bind_imagef("in_deposit", deposit.vector, nx, ny);
    run_shade.bind_imagef("in_moisture", moisture.vector, nx, ny);
    run_shade.bind_imagef("in_streak", streak.vector, nx, ny);
    run_shade.bind_imagef("in_rill", rill.vector, nx, ny);
    run_shade.bind_imagef("in_nl", n_large.vector, nx, ny);
    run_shade.bind_imagef("in_nm", n_mid.vector, nx, ny);
    run_shade.bind_imagef("in_nf", n_fine.vector, nx, ny);
    run_shade.bind_imagef("in_ng", n_grain.vector, nx, ny);
    run_shade.bind_imagef("in_snow_depth", snow_depth.vector, nx, ny);

    run_shade.bind_buffer<float>("out_r", out[0].vector);
    run_shade.bind_buffer<float>("out_g", out[1].vector);
    run_shade.bind_buffer<float>("out_b", out[2].vector);
    run_shade.bind_buffer<float>("out_a", out[3].vector);

    const float flow_width = 0.8f - 0.5f * p.flow_width;

    cl_float3 col_grass = {{p.grass.r, p.grass.g, p.grass.b, 0.f}};
    cl_float3 col_dry = {{p.dry.r, p.dry.g, p.dry.b, 0.f}};
    cl_float3 col_soil = {{p.soil.r, p.soil.g, p.soil.b, 0.f}};
    cl_float3 col_sediment = {{p.sediment.r, p.sediment.g, p.sediment.b, 0.f}};
    cl_float3 col_rock = {{p.rock.r, p.rock.g, p.rock.b, 0.f}};
    cl_float3 col_ridge = {{p.ridge.r, p.ridge.g, p.ridge.b, 0.f}};
    cl_float3 col_snow = {{p.snow.r, p.snow.g, p.snow.b, 0.f}};

    const int use_snow_depth = has_snow_depth ? 1 : 0;

    run_shade.bind_arguments(nx,
                             ny,
                             use_snow_depth,
                             zmin,
                             zrange,
                             p.rock_talus,
                             p.rock_edge,
                             p.strata,
                             p.snow_line,
                             p.snow_amount,
                             p.color_noise,
                             flow_width,
                             p.flow,
                             p.deposition,
                             p.streaks,
                             p.cavity,
                             p.detail,
                             col_grass,
                             col_dry,
                             col_soil,
                             col_sediment,
                             col_rock,
                             col_ridge,
                             col_snow);

    run_shade.execute({nx, ny});

    run_shade.read_buffer("out_r");
    run_shade.read_buffer("out_g");
    run_shade.read_buffer("out_b");
    run_shade.read_buffer("out_a");
  }

  return out;
}

} // namespace hmap
