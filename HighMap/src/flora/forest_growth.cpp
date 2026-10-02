/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

#include "point_sampler/metrics.hpp"

#include "highmap/flora/forest_growth.hpp"
#include "highmap/internal/validation.hpp"

#include <unordered_map>

namespace hmap
{

// ============================================================================
//  InteractionMatrix Implementation
// ============================================================================

InteractionMatrix::InteractionMatrix(size_t num_species, float default_val)
    : size(num_species), values(num_species * num_species, default_val)
{
}

void InteractionMatrix::fill(float value)
{
  std::fill(values.begin(), values.end(), value);
}

float InteractionMatrix::get(size_t s1, size_t s2) const
{
  if (s1 >= size || s2 >= size) return 0.0f;
  return values[s1 * size + s2];
}

void InteractionMatrix::set(size_t s1, size_t s2, float value)
{
  if (s1 < size && s2 < size) values[s1 * size + s2] = value;
}

void InteractionMatrix::set_symmetric(size_t s1, size_t s2, float value)
{
  set(s1, s2, value);
  set(s2, s1, value);
}

InteractionMatrix InteractionMatrix::diagonal(const std::vector<float> &diag,
                                              float off_diag_val)
{
  InteractionMatrix mat(diag.size(), off_diag_val);
  for (size_t i = 0; i < diag.size(); ++i)
    mat.set(i, i, diag[i]);
  return mat;
}

InteractionMatrix InteractionMatrix::from_radii(const std::vector<float> &radii,
                                                float multiplier)
{
  size_t            num_species = radii.size();
  InteractionMatrix mat(num_species);

  for (size_t s1 = 0; s1 < num_species; ++s1)
  {
    for (size_t s2 = 0; s2 < num_species; ++s2)
    {
      mat.set(s1, s2, multiplier * (radii[s1] + radii[s2]));
    }
  }

  return mat;
}

InteractionMatrix InteractionMatrix::from_species(
    const std::vector<Species> &species,
    float                       multiplier)
{
  size_t            num_species = species.size();
  InteractionMatrix mat(num_species);

  for (size_t s1 = 0; s1 < num_species; ++s1)
  {
    for (size_t s2 = 0; s2 < num_species; ++s2)
    {
      mat.set(s1, s2, multiplier * (species[s1].radius + species[s2].radius));
    }
  }

  return mat;
}

InteractionMatrix InteractionMatrix::random(size_t   num_species,
                                            uint32_t seed,
                                            float    random_offset,
                                            bool     symmetric)
{
  InteractionMatrix                     mat(num_species);
  std::mt19937                          rng(seed);
  std::uniform_real_distribution<float> dist(1.0f - random_offset,
                                             1.0f + random_offset);

  if (symmetric)
  {
    for (size_t i = 0; i < num_species; ++i)
    {
      mat.set(i, i, dist(rng));
      for (size_t j = i + 1; j < num_species; ++j)
      {
        float val = dist(rng);
        mat.set_symmetric(i, j, val);
      }
    }
  }
  else
  {
    for (size_t i = 0; i < num_species; ++i)
    {
      for (size_t j = 0; j < num_species; ++j)
      {
        mat.set(i, j, dist(rng));
      }
    }
  }

  return mat;
}

InteractionMatrix InteractionMatrix::uniform(size_t num_species, float val)
{
  return InteractionMatrix(num_species, val);
}

// ============================================================================
//  Growth & Thinning Functions (Alphabetically Sorted)
// ============================================================================

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

    float xn = (x - bbox.x) / bbox_dx *
               static_cast<float>(max_radius_scale.shape.x - 1);
    float yn = (y - bbox.z) / bbox_dy *
               static_cast<float>(max_radius_scale.shape.y - 1);

    int i = static_cast<int>(xn);
    int j = static_cast<int>(yn);

    if (i >= 0 && i < max_radius_scale.shape.x && j >= 0 &&
        j < max_radius_scale.shape.y)
    {
      float u = xn - static_cast<float>(i);
      float v = yn - static_cast<float>(j);
      float s_raw = std::clamp(
          max_radius_scale.get_value_bilinear_at(i, j, u, v),
          0.0f,
          1.0f);
      // blend: strength = 0 -> 1.0 (unchanged rmax), strength = 1 -> s_raw
      return 1.0f + strength * (s_raw - 1.0f);
    }
    return 1.0f;
  };

  // if only 1 tree, clamp according to its species if available
  if (forest.size() == 1)
  {
    Tree   t = forest[0];
    size_t sid = t.species_id;
    if (sid < species.size())
    {
      float s = sample_scale(t.position.x, t.position.y);
      float r_max_eff = species[sid].radius_min +
                        s * (species[sid].radius_max - species[sid].radius_min);
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
    float r_max_eff = sp_i.radius_min + s * (sp_i.radius_max - sp_i.radius_min);

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

Forest thin_forest_soft_core(const Forest            &forest,
                             const InteractionMatrix &repulsion_distances,
                             const InteractionMatrix &repulsion_strengths,
                             size_t                   target_count,
                             uint32_t                 seed,
                             const glm::vec4         &bbox)
{
  if (forest.empty()) return Forest();

  size_t max_target = (target_count == 0) ? forest.size() : target_count;

  // determine maximum interaction distance
  float max_r = 0.0f;
  for (float r : repulsion_distances.values)
    max_r = std::max(max_r, r);

  // prioritize candidates randomly
  std::mt19937                          rng(seed);
  std::uniform_real_distribution<float> dist_uni(0.0f, 1.0f);

  std::vector<size_t> indices(forest.size());
  std::iota(indices.begin(), indices.end(), 0);

  std::vector<float> priorities(forest.size());
  for (size_t i = 0; i < forest.size(); ++i)
    priorities[i] = dist_uni(rng);

  std::sort(indices.begin(),
            indices.end(),
            [&](size_t a, size_t b) { return priorities[a] < priorities[b]; });

  if (max_r <= 1e-6f)
  {
    std::vector<Tree> accepted;
    accepted.reserve(std::min(max_target, forest.size()));
    for (size_t i = 0; i < std::min(max_target, indices.size()); ++i)
      accepted.push_back(forest[indices[i]]);
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

          float r = repulsion_distances.get(p.species_id, q.species_id);
          float theta = repulsion_strengths.get(p.species_id, q.species_id);

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
      accepted.push_back(p);
      grid[make_key(gx, gy)].push_back(accepted.size() - 1);

      if (accepted.size() >= max_target) break;
    }
  }

  return Forest(std::move(accepted));
}

} // namespace hmap
