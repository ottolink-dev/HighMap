/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

#include "highmap/flora/forest_growth.hpp"
#include "highmap/internal/flora_utils.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math.hpp"

#include <unordered_map>

namespace hmap
{

Forest grow_forest_soft_core(const Forest               &forest,
                             const std::vector<Species> &species,
                             const InteractionMatrix    &repulsion_distances,
                             const InteractionMatrix    &repulsion_strengths,
                             size_t                      target_count,
                             const Array                &max_radius_scale,
                             float            max_radius_scale_strength,
                             uint32_t         seed,
                             const glm::vec4 &bbox)
{
  if (forest.empty()) return Forest();

  ForestScaleSampler sampler(max_radius_scale, max_radius_scale_strength, bbox);
  SpeciesLookup      lookup(species);

  size_t max_target = (target_count == 0) ? forest.size() : target_count;

  // build interaction distances matrix if empty from species
  InteractionMatrix eff_repulsion_dist = repulsion_distances;
  if (eff_repulsion_dist.size == 0 && !species.empty())
  {
    eff_repulsion_dist = InteractionMatrix::from_species(species, 2.0f);
  }

  InteractionMatrix eff_repulsion_str = repulsion_strengths;
  if (eff_repulsion_str.size == 0)
  {
    size_t num_sp = std::max<size_t>(1, species.size());
    eff_repulsion_str = InteractionMatrix::uniform(num_sp, 1.0f);
  }

  // determine maximum interaction distance
  float max_r = 0.0f;
  for (float r : eff_repulsion_dist.values)
  {
    max_r = std::max(max_r, r);
  }

  // if still 0, fallback to twice maximum tree radius
  if (max_r <= 1e-6f)
  {
    for (const auto &sp : species)
    {
      max_r = std::max(max_r, 2.0f * sp.radius_max);
    }
  }

  // prioritize candidates randomly
  std::mt19937                          rng(seed);
  std::uniform_real_distribution<float> dist_uni(0.0f, 1.0f);

  std::vector<size_t> indices(forest.size());
  std::iota(indices.begin(), indices.end(), 0);

  std::vector<float> priorities(forest.size());
  for (size_t i = 0; i < forest.size(); ++i)
  {
    priorities[i] = dist_uni(rng);
  }

  std::sort(indices.begin(),
            indices.end(),
            [&](size_t a, size_t b) { return priorities[a] < priorities[b]; });

  auto assign_final_radius = [&](const Tree &t) -> Tree
  {
    Tree    res = t;
    Species sp_i = lookup.get(t);
    float   s = sampler.sample(t.position.x, t.position.y);
    float   r_max_eff = lerp(sp_i.radius_min, sp_i.radius_max, s);
    res.radius = std::clamp(t.radius, sp_i.radius_min, r_max_eff);
    return res;
  };

  if (max_r <= 1e-6f)
  {
    std::vector<Tree> accepted;
    accepted.reserve(std::min(max_target, forest.size()));
    for (size_t i = 0; i < std::min(max_target, indices.size()); ++i)
    {
      accepted.push_back(assign_final_radius(forest[indices[i]]));
    }
    return Forest(std::move(accepted));
  }

  // --- Spatial Hash Grid for Neighbor Lookup

  float                                            cell_size = max_r;
  std::unordered_map<int64_t, std::vector<size_t>> grid;

  auto compute_cell = [&](float x, float y) -> std::pair<int, int>
  {
    int gx = static_cast<int>(std::floor((x - bbox.x) / cell_size));
    int gy = static_cast<int>(std::floor((y - bbox.z) / cell_size));
    return {gx, gy};
  };

  auto make_key = [](int gx, int gy) -> int64_t
  {
    return (static_cast<int64_t>(gx) << 32) ^
           (static_cast<int64_t>(gy) & 0xFFFFFFFF);
  };

  std::vector<Tree> accepted;
  accepted.reserve(std::min(max_target, forest.size()));

  for (size_t idx : indices)
  {
    const Tree &p = forest[idx];
    auto [gx, gy] = compute_cell(p.position.x, p.position.y);

    float survival_prob = 1.0f;

    for (int dy = -1; dy <= 1; ++dy)
    {
      for (int dx = -1; dx <= 1; ++dx)
      {
        int64_t key = make_key(gx + dx, gy + dy);
        auto    it = grid.find(key);
        if (it == grid.end()) continue;

        for (size_t accepted_idx : it->second)
        {
          const Tree &q = accepted[accepted_idx];

          float dist_sq = (p.position.x - q.position.x) *
                              (p.position.x - q.position.x) +
                          (p.position.y - q.position.y) *
                              (p.position.y - q.position.y);

          float r = eff_repulsion_dist.get(p.species_id, q.species_id);
          float theta = eff_repulsion_str.get(p.species_id, q.species_id);

          if (r > 1e-6f && theta > 1e-6f)
          {
            float r_sq = r * r;
            if (dist_sq < 9.0f * r_sq)
            {
              float factor = 1.0f - theta * std::exp(-dist_sq / r_sq);
              survival_prob *= std::max(0.0f, factor);
            }
          }
        }
      }
    }

    if (dist_uni(rng) < survival_prob)
    {
      accepted.push_back(assign_final_radius(p));
      grid[make_key(gx, gy)].push_back(accepted.size() - 1);

      if (accepted.size() >= max_target) break;
    }
  }

  return Forest(std::move(accepted));
}

} // namespace hmap
