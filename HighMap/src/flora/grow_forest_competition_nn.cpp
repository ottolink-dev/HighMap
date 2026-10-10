/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "point_sampler/metrics.hpp"
#include "point_sampler/point.hpp"

#include "highmap/array.hpp"
#include "highmap/flora/forest.hpp"
#include "highmap/flora/forest_growth.hpp"
#include "highmap/flora/species.hpp"
#include "highmap/flora/tree.hpp"
#include "highmap/internal/flora_utils.hpp"
#include "highmap/math/core.hpp"
#include "highmap/scatter/scatter_field.hpp"
#include "highmap/scatter/scatter_item.hpp"

namespace hmap
{

Forest grow_forest_competition_nn(const Forest               &forest,
                                  const std::vector<Species> &species,
                                  const InteractionMatrix &competition_matrix,
                                  const Array             &max_radius_scale,
                                  float            max_radius_scale_strength,
                                  bool             prune_unviable,
                                  const glm::vec4 &bbox)
{
  if (forest.empty()) return Forest();

  ForestScaleSampler sampler(max_radius_scale, max_radius_scale_strength, bbox);
  SpeciesLookup      lookup(species);

  // if only 1 tree, clamp according to its species if available
  if (forest.size() == 1)
  {
    Tree    t = forest[0];
    Species sp_i = lookup.get(t);
    float   s = sampler.sample(t.position.x, t.position.y);
    float   r_max_eff = lerp(sp_i.radius_min, sp_i.radius_max, s);
    t.radius = std::clamp(t.radius, sp_i.radius_min, r_max_eff);
    return Forest(std::vector<Tree>{t});
  }

  // --- Build 2D Point List for Nearest Neighbor Query

  std::vector<ps::Point<float, 2>> points;
  points.reserve(forest.size());
  for (const auto &tree : forest)
  {
    points.push_back({tree.position.x, tree.position.y});
  }

  auto neighbors_idx = ps::nearest_neighbors_indices(points, 1);

  std::vector<Tree> result;
  result.reserve(forest.size());

  for (size_t i = 0; i < forest.size(); ++i)
  {
    const Tree &tree = forest[i];
    Species     sp_i = lookup.get(tree);

    // find distance to closest neighbor
    size_t neighbor_idx = neighbors_idx[i].empty() ? i : neighbors_idx[i][0];
    const Tree &neighbor_tree = forest[neighbor_idx];

    float dx = tree.position.x - neighbor_tree.position.x;
    float dy = tree.position.y - neighbor_tree.position.y;
    float d_nn = std::sqrt(dx * dx + dy * dy);

    // determine pairwise competition factor alpha
    float alpha = lookup.get_alpha(sp_i,
                                   neighbor_tree.class_id,
                                   competition_matrix);

    // scale max_radius locally if array is provided
    float s = sampler.sample(tree.position.x, tree.position.y);
    float r_max_eff = lerp(sp_i.radius_min, sp_i.radius_max, s);

    // estimated crown radius
    float r_est = alpha * d_nn;

    // check viability
    if (prune_unviable && r_est < sp_i.radius_min)
    {
      continue; // cull choked plant
    }

    Tree grown_tree = tree;
    grown_tree.radius = std::clamp(r_est, sp_i.radius_min, r_max_eff);
    result.push_back(grown_tree);
  }

  return Forest(std::move(result));
}

} // namespace hmap
