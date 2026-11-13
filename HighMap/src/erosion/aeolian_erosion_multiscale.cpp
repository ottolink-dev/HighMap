/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <cmath>

#include "highmap/array.hpp"
#include "highmap/erosion/aeolian_erosion.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap
{

void aeolian_erosion_multiscale(Array        &z,
                                int           steps,
                                std::uint32_t seed,
                                glm::vec2     wind,
                                const Array  *p_bedrock,
                                const Array  *p_vegetation,
                                const Array  *p_abrasion_mask,
                                Array        *p_sediment_map,
                                Array        *p_bedrock_map,
                                float         cell_size,
                                float         matter_to_move,
                                float         sand_repose_angle,
                                float         shadow_angle_min,
                                float         shadow_angle_max,
                                float         max_shadow_dist,
                                float         abrasion_rate,
                                float         repose_rate,
                                int           levels,
                                float         coarse_fraction,
                                int           min_level_size,
                                int           shadow_map_interval,
                                int           repose_interval)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_vegetation && !validate_same_shape(z, *p_vegetation)) return;
  if (p_abrasion_mask && !validate_same_shape(z, *p_abrasion_mask)) return;

  const glm::ivec2 shape = z.shape;

  if (cell_size <= 0.f) cell_size = 1.f / static_cast<float>(shape.x);

  // --- Full-resolution layers

  Array bed(shape, 0.f);
  if (p_bedrock) bed = *p_bedrock;

  Array sed = z - bed;
  clamp_min(sed, 0.f);

  // --- Number of usable levels (level 0 is the full resolution)

  int n_levels = std::max(1, levels);

  while (n_levels > 1)
  {
    const int factor = 1 << (n_levels - 1);

    if (std::min(shape.x, shape.y) / factor >= min_level_size) break;

    n_levels--;
  }

  // --- Steps distribution

  coarse_fraction = std::clamp(coarse_fraction, 0.f, 0.95f);

  const int n_coarse = n_levels - 1;
  const int steps_coarse_total = (n_coarse > 0)
                                     ? static_cast<int>(coarse_fraction *
                                                        static_cast<float>(
                                                            steps))
                                     : 0;
  const int steps_per_coarse = (n_coarse > 0) ? steps_coarse_total / n_coarse
                                              : 0;

  int steps_used = 0;

  // --- Loop from the coarsest level to the full resolution

  for (int level = n_levels - 1; level >= 0; --level)
  {
    const bool is_finest = (level == 0);
    const int  factor = 1 << level;

    const int level_steps = is_finest ? std::max(0, steps - steps_used)
                                      : steps_per_coarse;

    steps_used += is_finest ? 0 : level_steps;

    if (level_steps <= 0) continue;

    // seed offset so that levels do not reuse the same random sequence
    const std::uint32_t level_seed = seed + 7919u *
                                                static_cast<std::uint32_t>(
                                                    level);

    if (is_finest)
    {
      Array z_fine = bed + sed;

      aeolian_erosion(z_fine,
                      level_steps,
                      level_seed,
                      wind,
                      &bed,
                      p_vegetation,
                      p_abrasion_mask,
                      &sed,
                      &bed,
                      cell_size,
                      matter_to_move,
                      0,
                      sand_repose_angle,
                      shadow_angle_min,
                      shadow_angle_max,
                      max_shadow_dist,
                      abrasion_rate,
                      repose_rate,
                      0,
                      TileRegion::dummy(),
                      SingleArrayComm::instance(),
                      shadow_map_interval,
                      repose_interval);

      continue;
    }

    // --- Coarse level: downsample the layers

    const glm::ivec2 shape_l = {std::max(1, shape.x / factor),
                                std::max(1, shape.y / factor)};

    Array bed_l = bed.resample_to_shape(shape_l);
    Array sed_l = sed.resample_to_shape(shape_l);
    clamp_min(sed_l, 0.f);

    Array bed_l0 = bed_l;
    Array sed_l0 = sed_l;

    Array veg_l;
    Array mask_l;

    if (p_vegetation) veg_l = p_vegetation->resample_to_shape(shape_l);
    if (p_abrasion_mask) mask_l = p_abrasion_mask->resample_to_shape(shape_l);

    Array z_l = bed_l + sed_l;

    aeolian_erosion(z_l,
                    level_steps,
                    level_seed,
                    wind,
                    &bed_l,
                    p_vegetation ? &veg_l : nullptr,
                    p_abrasion_mask ? &mask_l : nullptr,
                    &sed_l,
                    &bed_l,
                    cell_size * static_cast<float>(factor),
                    matter_to_move,
                    0,
                    sand_repose_angle,
                    shadow_angle_min,
                    shadow_angle_max,
                    max_shadow_dist,
                    abrasion_rate,
                    repose_rate,
                    0,
                    TileRegion::dummy(),
                    SingleArrayComm::instance(),
                    shadow_map_interval,
                    repose_interval);

    // --- Upsample only the changes and apply them to the full resolution

    Array d_sed = (sed_l - sed_l0).resample_to_shape(shape);
    sed += d_sed;
    clamp_min(sed, 0.f);

    if (abrasion_rate > 0.f)
    {
      Array d_bed = (bed_l - bed_l0).resample_to_shape(shape);
      bed += d_bed;
    }
  }

  // --- Outputs

  z = bed + sed;

  if (p_sediment_map) *p_sediment_map = sed;
  if (p_bedrock_map) *p_bedrock_map = bed;
}

} // namespace hmap
