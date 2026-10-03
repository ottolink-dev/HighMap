/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include <algorithm>
#include <cstdint>

namespace hmap
{

/**
 * @brief Parameters for the MISE (Multiscale Implicit Stream-power Erosion)
 * solver.
 *
 * Adapted from MISE by Leonhard (https://github.com/Leonhardmaster2).
 */
struct MiseParams
{
  // --- Fluvial incision

  float strength = 0.5f;   // overall erosion power (K x duration)
  float area_exp = 0.45f;  // drainage area exponent m
  float downcutting = 1.f; // cap on normalised drainage area A in (0, 1]

  // --- Hillslopes

  float thermal = 0.05f; // rate x duration of plain relaxation
  float debris = 200.f;  // extra rate per unit of upstream steep area
  float talus = 3.2f;    // critical slope
  float cliff = 6.8f;    // hard slope limit along flow lines (0 = none)

  // --- Deposition

  float deposition_slope = 1.8f;    // slope below which sediment settles
  float deposition_exp = 0.1f;      // decay of equilibrium slope with drainage
                                    // area
  float deposition_area = 1e-4f;    // normalised area where decay starts
  float deposition_rate = 1.f;      // share of available room filled per step
  float lake_fill = 1.f;            // share of sediment reaching closed
                                    // depression that stays in it
  float sediment_erodibility = 2.f; // erodibility multiplier for loose
                                    // sediment
  float sediment_talus = 1.4f;      // repose slope of loose sediment
  int   sediment_relax_iters = 2;   // sediment talus sweeps per step (0 =
                                    // off)

  // --- Flow routing

  float receiver_exp = 2.f;         // receiver drop probability exponent
  float slope_correction_max = 2.f; // clamp of |grad z| / directional
                                    // slope

  // --- Multi-scale schedule

  int   base_res = 128;      // coarsest level resolution across
  float share_ratio = 0.6f;  // strength ratio geometric decay across
                             // levels
  int   iters0 = 12;         // flow-routing updates at coarsest level
  float iters_decay = 0.75f; // iteration decay per level
  int   iters_min = 3;       // minimum iterations per level
  int   skip_finest = 0;     // number of finest levels skipped
  float warp = 0.35f;        // domain warp when upsampling
  float antialias = 0.5f;    // binomial blur factor at full resolution

  bool extrapolate_border = true; // extrapolate border cells from
                                  // interior
  std::uint32_t seed = 1;         // random seed
  int           threads = 0;      // thread count (0 = hardware
                                  // concurrency)
  bool exact_flood = false;       // exact priority queue instead of
                                  // buckets
};

} // namespace hmap

namespace hmap::gpu
{

/**
 * @brief Parameters for McDonald particle-based hydraulic erosion solver.
 */
struct McDonaldParams
{
  float strength = 0.5f;     // overall erosion power (scales suspension &
                             // thermal rates)
  float deposition = 0.5f;   // sediment retention vs transport
  float crit_slope = 0.57f;  // critical slope [m/m]
  float meandering = 0.5f;   // flow-coupling & momentum inertia
  float scale = 1.0f;        // domain extent multiplier (scales world_extent_km
                             // and base z_scale_km)
  float relief_scale = 1.0f; // relative vertical relief scale multiplier

  struct PhysicalParams
  {
    float world_extent_km;
    float z_scale_km;
    int   samples;
    int   maxage;
    float lrate;
    float time_step;
    float rainfall;
    float evap_rate;
    float gravity;
    float viscosity;
    float bed_shear;
    float crit_slope;
    float settle_rate;
    float thermal_rate;
    float deposition_rate;
    float suspension_rate;
    float exit_slope;
  };

  PhysicalParams to_physical() const
  {
    PhysicalParams p;
    float          s = std::max(1e-4f, this->scale);
    float          r = std::max(1e-4f, this->relief_scale);
    p.world_extent_km = 40.f * s;
    p.z_scale_km = 4.f * s * r;
    p.samples = 8192;
    p.maxage = 512;
    p.lrate = 0.1f + 0.8f * std::clamp(this->meandering, 0.f, 1.f);
    p.time_step = 10.f;
    p.rainfall = 1.f;
    p.evap_rate = 1e-9f;
    p.gravity = 9.81f;
    p.viscosity = 0.01f + 0.03f * std::clamp(this->meandering, 0.f, 1.f);
    p.bed_shear = 0.02f - 0.015f * std::clamp(this->meandering, 0.f, 1.f);
    p.crit_slope = std::max(1e-4f, this->crit_slope);
    p.settle_rate = 0.1f *
                    (0.2f + 1.6f * std::clamp(this->deposition, 0.f, 1.f));
    p.thermal_rate = 2.5e-3f *
                     (0.2f + 1.6f * std::clamp(this->strength, 0.f, 1.f));
    p.deposition_rate = 5e-3f *
                        (0.2f + 1.6f * std::clamp(this->deposition, 0.f, 1.f));
    p.suspension_rate = 2.5e-4f *
                        (0.2f + 1.6f * std::clamp(this->strength, 0.f, 1.f));
    p.exit_slope = 0.01f;
    return p;
  }
};

} // namespace hmap::gpu
