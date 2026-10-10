/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

/* Clean-room port of the flow-coupled particle erosion model from
 * erosiv/soillib (LGPL-3, Nicholas McDonald): persistent discharge/momentum
 * fields couple particles through the mean flow; a bank-stability debris
 * flow runs in the same solver loop against a separate sediment layer. */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/erosion/erosion_parameters.hpp"
#include "highmap/erosion/hydraulic_erosion.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"

namespace hmap::gpu
{

namespace detail
{

// Runs `steps` erosion iterations on the given model state (all arrays same
// shape). State lives in the host vectors between kernel invocations
// (CLWrapper Runs own their device buffers).
void mcdonald_run_steps(Array        &bed,
                        Array        &sed,
                        Array        &dis,
                        Array        &mx,
                        Array        &my,
                        const Array  *p_moisture_map,
                        int           steps,
                        std::uint32_t seed,
                        float         world_extent_km,
                        float         z_scale_km,
                        int           samples,
                        int           maxage,
                        float         lrate,
                        float         time_step,
                        float         rainfall,
                        float         evap_rate,
                        float         gravity,
                        float         viscosity,
                        float         bed_shear,
                        float         crit_slope,
                        float         settle_rate,
                        float         thermal_rate,
                        float         deposition_rate,
                        float         suspension_rate,
                        float         exit_slope)
{
  int nx = bed.shape.x;
  int ny = bed.shape.y;

  const float cell_m = world_extent_km * 1e3f / (float)nx;
  const float z_m = z_scale_km * 1e3f;

  Array tr_d(bed.shape), tr_mx(bed.shape), tr_my(bed.shape);

  for (int step = 0; step < steps; ++step)
  {
    // reset (host-side, replaces soillib's reset kernel — algebraically
    // identical since the tracks are re-uploaded fresh each step)
    std::fill(tr_d.vector.begin(), tr_d.vector.end(), 0.f);
    std::fill(tr_mx.vector.begin(), tr_mx.vector.end(), 0.f);
    std::fill(tr_my.vector.begin(), tr_my.vector.end(), 0.f);

    { // solve
      auto run = clwrapper::Run("mcdonald_solve");
      run.bind_buffer<float>("bed", bed.vector);
      run.bind_buffer<float>("sed", sed.vector);
      run.bind_buffer<float>("dis", dis.vector);
      run.bind_buffer<float>("mx", mx.vector);
      run.bind_buffer<float>("my", my.vector);
      run.bind_buffer<float>("tr_d", tr_d.vector);
      run.bind_buffer<float>("tr_mx", tr_mx.vector);
      run.bind_buffer<float>("tr_my", tr_my.vector);
      helper_bind_optional_buffer(run, "moisture_map", p_moisture_map);
      run.bind_arguments(nx,
                         ny,
                         samples,
                         seed,
                         (std::uint32_t)step,
                         cell_m,
                         z_m,
                         time_step,
                         rainfall,
                         evap_rate,
                         gravity,
                         viscosity,
                         bed_shear,
                         deposition_rate,
                         suspension_rate,
                         exit_slope,
                         maxage,
                         p_moisture_map ? 1 : 0);
      run.write_buffer("bed");
      run.write_buffer("sed");
      run.write_buffer("dis");
      run.write_buffer("mx");
      run.write_buffer("my");
      run.write_buffer("tr_d");
      run.write_buffer("tr_mx");
      run.write_buffer("tr_my");
      run.execute(samples);
      run.read_buffer("bed");
      run.read_buffer("sed");
      run.read_buffer("tr_d");
      run.read_buffer("tr_mx");
      run.read_buffer("tr_my");
    }

    { // filter (exponential mix of tracks into flow fields)
      auto run = clwrapper::Run("mcdonald_filter");
      run.bind_buffer<float>("dis", dis.vector);
      run.bind_buffer<float>("mx", mx.vector);
      run.bind_buffer<float>("my", my.vector);
      run.bind_buffer<float>("tr_d", tr_d.vector);
      run.bind_buffer<float>("tr_mx", tr_mx.vector);
      run.bind_buffer<float>("tr_my", tr_my.vector);
      run.bind_arguments(nx, ny, lrate);
      run.write_buffer("dis");
      run.write_buffer("mx");
      run.write_buffer("my");
      run.write_buffer("tr_d");
      run.write_buffer("tr_mx");
      run.write_buffer("tr_my");
      run.execute({nx, ny});
      run.read_buffer("dis");
      run.read_buffer("mx");
      run.read_buffer("my");
    }

    { // debris flow
      auto run = clwrapper::Run("mcdonald_debris");
      run.bind_buffer<float>("bed", bed.vector);
      run.bind_buffer<float>("sed", sed.vector);
      run.bind_arguments(nx,
                         ny,
                         samples,
                         seed,
                         (std::uint32_t)step,
                         cell_m,
                         z_m,
                         time_step,
                         gravity,
                         crit_slope,
                         settle_rate,
                         thermal_rate);
      run.write_buffer("bed");
      run.write_buffer("sed");
      run.execute(samples);
      run.read_buffer("bed");
      run.read_buffer("sed");
    }
  }
}

} // namespace detail

void hydraulic_mcdonald(Array                &z,
                        int                   steps,
                        std::uint32_t         seed,
                        const McDonaldParams &params,
                        const Array          *p_moisture_map,
                        Array                *p_sediment_map,
                        Array                *p_discharge_map)
{
  auto p = params.to_physical();
  hydraulic_mcdonald(z,
                     steps,
                     seed,
                     p_moisture_map,
                     p_sediment_map,
                     p_discharge_map,
                     p.world_extent_km,
                     p.z_scale_km,
                     p.samples,
                     p.maxage,
                     p.lrate,
                     p.time_step,
                     p.rainfall,
                     p.evap_rate,
                     p.gravity,
                     p.viscosity,
                     p.bed_shear,
                     p.crit_slope,
                     p.settle_rate,
                     p.thermal_rate,
                     p.deposition_rate,
                     p.suspension_rate,
                     p.exit_slope);
}

void hydraulic_mcdonald(Array        &z,
                        int           steps,
                        std::uint32_t seed,
                        const Array  *p_moisture_map,
                        Array        *p_sediment_map,
                        Array        *p_discharge_map,
                        float         world_extent_km,
                        float         z_scale_km,
                        int           samples,
                        int           maxage,
                        float         lrate,
                        float         time_step,
                        float         rainfall,
                        float         evap_rate,
                        float         gravity,
                        float         viscosity,
                        float         bed_shear,
                        float         crit_slope,
                        float         settle_rate,
                        float         thermal_rate,
                        float         deposition_rate,
                        float         suspension_rate,
                        float         exit_slope)
{
  if (!validate_non_empty(z)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;

  Array sed(z.shape), dis(z.shape), mx(z.shape), my(z.shape);

  detail::mcdonald_run_steps(z,
                             sed,
                             dis,
                             mx,
                             my,
                             p_moisture_map,
                             steps,
                             seed,
                             world_extent_km,
                             z_scale_km,
                             samples,
                             maxage,
                             lrate,
                             time_step,
                             rainfall,
                             evap_rate,
                             gravity,
                             viscosity,
                             bed_shear,
                             crit_slope,
                             settle_rate,
                             thermal_rate,
                             deposition_rate,
                             suspension_rate,
                             exit_slope);

  if (p_sediment_map) *p_sediment_map = sed;
  if (p_discharge_map) *p_discharge_map = dis;

  // total surface = bedrock + sediment
  for (size_t k = 0; k < z.vector.size(); ++k)
    z.vector[k] += sed.vector[k];
}

void hydraulic_mcdonald_multiscale(Array                  &z,
                                   std::uint32_t           seed,
                                   const std::vector<int> &steps_per_level,
                                   const McDonaldParams   &params,
                                   const Array            *p_moisture_map,
                                   Array                  *p_sediment_map,
                                   Array                  *p_discharge_map)
{
  auto p = params.to_physical();
  hydraulic_mcdonald_multiscale(z,
                                seed,
                                steps_per_level,
                                p_moisture_map,
                                p_sediment_map,
                                p_discharge_map,
                                p.world_extent_km,
                                p.z_scale_km,
                                p.samples,
                                p.maxage,
                                p.lrate,
                                p.time_step,
                                p.rainfall,
                                p.evap_rate,
                                p.gravity,
                                p.viscosity,
                                p.bed_shear,
                                p.crit_slope,
                                p.settle_rate,
                                p.thermal_rate,
                                p.deposition_rate,
                                p.suspension_rate,
                                p.exit_slope);
}

void hydraulic_mcdonald_multiscale(Array                  &z,
                                   std::uint32_t           seed,
                                   const std::vector<int> &steps_per_level,
                                   const Array            *p_moisture_map,
                                   Array                  *p_sediment_map,
                                   Array                  *p_discharge_map,
                                   float                   world_extent_km,
                                   float                   z_scale_km,
                                   int                     samples,
                                   int                     maxage,
                                   float                   lrate,
                                   float                   time_step,
                                   float                   rainfall,
                                   float                   evap_rate,
                                   float                   gravity,
                                   float                   viscosity,
                                   float                   bed_shear,
                                   float                   crit_slope,
                                   float                   settle_rate,
                                   float                   thermal_rate,
                                   float                   deposition_rate,
                                   float                   suspension_rate,
                                   float                   exit_slope)
{
  if (!validate_non_empty(z)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;

  int nlevels = (int)steps_per_level.size();
  if (nlevels == 0) return;

  // halving ladder, coarsest first; final level == z.shape
  std::vector<glm::ivec2> ladder(nlevels);
  for (int i = 0; i < nlevels; ++i)
  {
    int shift = nlevels - 1 - i;
    ladder[i] = {std::max(2, z.shape.x >> shift),
                 std::max(2, z.shape.y >> shift)};
  }

  Array bed = z.resample_to_shape(ladder[0]);
  Array sed(ladder[0]), dis(ladder[0]), mx(ladder[0]), my(ladder[0]);

  Array moist;
  if (p_moisture_map) moist = p_moisture_map->resample_to_shape(ladder[0]);

  for (int i = 0; i < nlevels; ++i)
  {
    if (i > 0)
    {
      bed = bed.resample_to_shape(ladder[i]);
      sed = sed.resample_to_shape(ladder[i]);
      dis = dis.resample_to_shape(ladder[i]);
      mx = mx.resample_to_shape(ladder[i]);
      my = my.resample_to_shape(ladder[i]);
      if (p_moisture_map) moist = p_moisture_map->resample_to_shape(ladder[i]);
    }

    detail::mcdonald_run_steps(bed,
                               sed,
                               dis,
                               mx,
                               my,
                               p_moisture_map ? &moist : nullptr,
                               steps_per_level[i],
                               seed + (std::uint32_t)i,
                               world_extent_km,
                               z_scale_km,
                               samples,
                               maxage,
                               lrate,
                               time_step,
                               rainfall,
                               evap_rate,
                               gravity,
                               viscosity,
                               bed_shear,
                               crit_slope,
                               settle_rate,
                               thermal_rate,
                               deposition_rate,
                               suspension_rate,
                               exit_slope);
  }

  if (p_sediment_map) *p_sediment_map = sed;
  if (p_discharge_map) *p_discharge_map = dis;

  z = bed;
  for (size_t k = 0; k < z.vector.size(); ++k)
    z.vector[k] += sed.vector[k];
}

} // namespace hmap::gpu
