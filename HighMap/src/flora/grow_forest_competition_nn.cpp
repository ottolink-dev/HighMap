/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "point_sampler/metrics.hpp"

#include "highmap/flora/forest_growth.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math.hpp"

#include <unordered_map>

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

  bool has_scale_array = !max_radius_scale.vector.empty() &&
                         validate_non_empty(max_radius_scale);
  float bbox_dx = bbox.y - bbox.x;
  float bbox_dy = bbox.w - bbox.z;
  bool  valid_bbox = (std::abs(bbox_dx) > 1e-7f && std::abs(bbox_dy) > 1e-7f);
  float strength = std::clamp(max_radius_scale_strength, 0.0f, 1.0f);

  auto sample_scale = [&](float x, float y) -> float
  {
    if (!has_scale_array || !valid_bbox || strength <= 0.0f) return 1.0f;

    float u_norm = std::clamp((x - bbox.x) / bbox_dx, 0.0f, 1.0f);
    float v_norm = std::clamp((y - bbox.z) / bbox_dy, 0.0f, 1.0f);

    float xn = u_norm * static_cast<float>(max_radius_scale.shape.x - 1);
    float yn = v_norm * static_cast<float>(max_radius_scale.shape.y - 1);

    int i = static_cast<int>(xn);
    int j = static_cast<int>(yn);

    float u = xn - static_cast<float>(i);
    float v = yn - static_cast<float>(j);
    float s_raw = std::clamp(max_radius_scale.get_value_bilinear_at(i, j, u, v),
                             0.0f,
                             1.0f);
    // blend: strength = 0 -> 1.0 (unchanged rmax), strength = 1 -> s_raw
    return lerp(1.f, s_raw, strength);
  };

  // if only 1 tree, clamp according to its species if available
  if (forest.size() == 1)
  {
    Tree   t = forest[0];
    size_t sid = t.species_id;
    if (sid < species.size())
    {
      float s = sample_scale(t.position.x, t.position.y);
      float r_max_eff = lerp(species[sid].radius_min,
                             species[sid].radius_max,
                             s);
      t.radius = std::clamp(t.radius, species[sid].radius_min, r_max_eff);
    }
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

  // --- Map Species Definitions for Fast Lookup

  std::unordered_map<uint32_t, Species> species_map;
  for (const auto &sp : species)
  {
    species_map[sp.id] = sp;
  }

  std::vector<Tree> result;
  result.reserve(forest.size());

  for (size_t i = 0; i < forest.size(); ++i)
  {
    const Tree &tree = forest[i];
    size_t      s_i = tree.species_id;

    // lookup species i traits
    Species sp_i(tree.species_id, tree.radius);
    auto    it_i = species_map.find(tree.species_id);
    if (it_i != species_map.end())
    {
      sp_i = it_i->second;
    }
    else if (s_i < species.size())
    {
      sp_i = species[s_i];
    }

    // find distance to closest neighbor
    size_t neighbor_idx = neighbors_idx[i].empty() ? i : neighbors_idx[i][0];
    const Tree &neighbor_tree = forest[neighbor_idx];
    size_t      s_j = neighbor_tree.species_id;

    float dx = tree.position.x - neighbor_tree.position.x;
    float dy = tree.position.y - neighbor_tree.position.y;
    float d_nn = std::sqrt(dx * dx + dy * dy);

    // determine pairwise competition factor alpha
    float alpha = sp_i.competition_factor;
    if (competition_matrix.size > 0)
    {
      float mat_val = competition_matrix.get(s_i, s_j);
      if (mat_val > 0.0f) alpha = mat_val;
    }

    // scale max_radius locally if array is provided
    float s = sample_scale(tree.position.x, tree.position.y);
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
