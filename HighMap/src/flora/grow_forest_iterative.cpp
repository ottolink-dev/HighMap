/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "point_sampler/metrics.hpp"

#include "highmap/flora/forest_growth.hpp"
#include "highmap/internal/flora_utils.hpp"
#include "highmap/math.hpp"

namespace hmap
{

Forest grow_forest_iterative(const Forest               &forest,
                             const std::vector<Species> &species,
                             size_t                      iterations,
                             float                       growth_rate,
                             const InteractionMatrix    &competition_matrix,
                             const Array                &max_radius_scale,
                             float            max_radius_scale_strength,
                             bool             prune_unviable,
                             const glm::vec4 &bbox)
{
  if (forest.empty()) return Forest();

  ForestScaleSampler sampler(max_radius_scale, max_radius_scale_strength, bbox);
  SpeciesLookup      lookup(species);

  // if only 1 tree, clamp and return
  if (forest.size() == 1)
  {
    Tree    t = forest[0];
    Species sp_i = lookup.get(t);
    float   s = sampler.sample(t.position.x, t.position.y);
    float   r_max_eff = lerp(sp_i.radius_min, sp_i.radius_max, s);
    t.radius = std::clamp(t.radius, sp_i.radius_min, r_max_eff);
    return Forest(std::vector<Tree>{t});
  }

  std::vector<Tree> current_trees;
  current_trees.reserve(forest.size());
  for (const auto &t : forest)
  {
    current_trees.push_back(t);
  }

  float rate = std::clamp(growth_rate, 0.0f, 1.0f);
  if (rate <= 1e-6f || iterations == 0)
  {
    return Forest(std::move(current_trees));
  }

  // --- Iterative Growth & Competition Simulation Loop

  for (size_t iter = 0; iter < iterations; ++iter)
  {
    if (current_trees.empty()) break;

    // build 2d point list for spatial neighbor query
    std::vector<ps::Point<float, 2>> points;
    points.reserve(current_trees.size());
    for (const auto &tree : current_trees)
    {
      points.push_back({tree.position.x, tree.position.y});
    }

    // query nearest neighbors (up to 8 neighbors)
    size_t safe_k = std::min<size_t>(8, current_trees.size() - 1);
    auto   neighbors_idx = ps::nearest_neighbors_indices(points, safe_k);

    std::vector<Tree> next_trees;
    next_trees.reserve(current_trees.size());

    float max_delta_r = 0.0f;

    for (size_t i = 0; i < current_trees.size(); ++i)
    {
      const Tree &tree = current_trees[i];
      Species     sp_i = lookup.get(tree);

      // scale max radius locally from map
      float s = sampler.sample(tree.position.x, tree.position.y);
      float r_max_eff = lerp(sp_i.radius_min, sp_i.radius_max, s);

      // intrinsic potential growth increment
      float growth_potential = rate * (r_max_eff - sp_i.radius_min);

      // --- Evaluate Neighborhood Competition C_i

      float competition = 0.0f;

      for (size_t neighbor_idx : neighbors_idx[i])
      {
        const Tree &neighbor_tree = current_trees[neighbor_idx];
        float       dx = tree.position.x - neighbor_tree.position.x;
        float       dy = tree.position.y - neighbor_tree.position.y;
        float       dist = std::sqrt(dx * dx + dy * dy);

        float alpha = lookup.get_alpha(sp_i,
                                       neighbor_tree.class_id,
                                       competition_matrix);

        // crown overlap / encroachment interaction
        float crown_reach = tree.radius + neighbor_tree.radius;
        if (dist < crown_reach && tree.radius > 1e-6f)
        {
          float overlap = (crown_reach - dist) / tree.radius;
          competition += alpha * std::max(0.0f, overlap);
        }
      }

      float comp_clamped = std::clamp(competition, 0.0f, 1.0f);

      // --- Update Crown Radius

      float r_new = tree.radius + growth_potential * (1.0f - comp_clamped);
      r_new = std::clamp(r_new, sp_i.radius_min, r_max_eff);

      max_delta_r = std::max(max_delta_r, std::abs(r_new - tree.radius));

      // --- Check Suppression and Mortality

      if (prune_unviable && r_new < sp_i.radius_min)
      {
        continue; // suppressed plant dies
      }

      Tree grown_tree = tree;
      grown_tree.radius = r_new;
      next_trees.push_back(grown_tree);
    }

    current_trees = std::move(next_trees);

    // stop if system converged
    if (max_delta_r < 1e-5f)
    {
      break;
    }
  }

  Forest result(std::move(current_trees));

  if (prune_unviable)
  {
    result.prune_unviable(species, true);
  }

  return result;
}

} // namespace hmap
