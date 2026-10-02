/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

#include "highmap/flora/forest_seeding.hpp"
#include "highmap/geometry/inverse_sampler_2d.hpp"
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

// --- Helper Functions

static float sample_array_bilinear(const Array &arr, float un, float vn)
{
  if (!validate_non_empty(arr)) return 0.0f;
  float xn = std::clamp(un, 0.0f, 1.0f) * static_cast<float>(arr.shape.x - 1);
  float yn = std::clamp(vn, 0.0f, 1.0f) * static_cast<float>(arr.shape.y - 1);
  int   i = std::min(static_cast<int>(xn), arr.shape.x - 1);
  int   j = std::min(static_cast<int>(yn), arr.shape.y - 1);
  float u = xn - static_cast<float>(i);
  float v = yn - static_cast<float>(j);
  return arr.get_value_bilinear_at(i, j, u, v);
}

// ============================================================================
//  Seeding & Thinning Functions (Alphabetically Sorted)
// ============================================================================

Forest seed_forest_clusters(size_t                      species_count,
                            size_t                      tree_count,
                            const Array                &density,
                            const Array                &exclusion,
                            float                       cluster_spread,
                            size_t                      points_per_cluster,
                            const ForestSeedingOptions &options)
{
  // validate inputs
  if (species_count == 0 || tree_count == 0 || points_per_cluster == 0 ||
      !validate_non_empty(density))
    return Forest();

  float width = options.bbox.y - options.bbox.x;
  float height = options.bbox.w - options.bbox.z;
  if (width <= 1e-7f || height <= 1e-7f) return Forest();

  // --- Modulate Density by Exclusion Map

  Array effective_density = density;
  bool  has_exclusion = validate_non_empty(exclusion);

  if (has_exclusion)
  {
    if (exclusion.shape == density.shape)
    {
      for (int j = 0; j < density.shape.y; ++j)
        for (int i = 0; i < density.shape.x; ++i)
        {
          if (exclusion(i, j) >= options.exclusion_threshold)
            effective_density(i, j) = 0.0f;
        }
    }
    else
    {
      for (int j = 0; j < density.shape.y; ++j)
      {
        float v_coord = static_cast<float>(j) /
                        static_cast<float>(std::max(1, density.shape.y - 1));
        for (int i = 0; i < density.shape.x; ++i)
        {
          float u_coord = static_cast<float>(i) /
                          static_cast<float>(std::max(1, density.shape.x - 1));
          float ex_val = sample_array_bilinear(exclusion, u_coord, v_coord);
          if (ex_val >= options.exclusion_threshold)
            effective_density(i, j) = 0.0f;
        }
      }
    }
  }

  InverseSampler2D sampler(effective_density, options.seed, options.bbox);
  if (sampler.get_total_weight() <= 1e-7f) return Forest();

  // --- Determine Species Quotas

  std::vector<size_t> species_tree_counts(species_count, 0);
  if (options.species_weights.size() == species_count)
  {
    float total_w = 0.0f;
    for (float w : options.species_weights)
      total_w += std::max(0.0f, w);

    if (total_w > 1e-6f)
    {
      size_t allocated = 0;
      for (size_t s = 0; s < species_count; ++s)
      {
        float norm_w = std::max(0.0f, options.species_weights[s]) / total_w;
        species_tree_counts[s] = static_cast<size_t>(
            std::round(static_cast<float>(tree_count) * norm_w));
        allocated += species_tree_counts[s];
      }

      // adjust difference due to rounding
      if (allocated < tree_count)
        species_tree_counts[0] += (tree_count - allocated);
      else if (allocated > tree_count &&
               species_tree_counts[0] >= (allocated - tree_count))
        species_tree_counts[0] -= (allocated - tree_count);
    }
  }
  else
  {
    size_t base_count = tree_count / species_count;
    size_t remainder = tree_count % species_count;
    for (size_t s = 0; s < species_count; ++s)
      species_tree_counts[s] = base_count + (s < remainder ? 1 : 0);
  }

  // --- Generate Clusters per Species

  std::vector<Tree> candidate_trees;
  candidate_trees.reserve(tree_count);

  for (size_t s = 0; s < species_count; ++s)
  {
    if (species_tree_counts[s] == 0) continue;

    size_t target_s = species_tree_counts[s];
    size_t cluster_count = std::max<size_t>(
        1,
        (target_s + points_per_cluster - 1) / points_per_cluster);

    auto parent_samples = sampler.sample(cluster_count);

    float default_radius = (s < options.species_radii.size())
                               ? options.species_radii[s]
                               : 1.0f;

    size_t count_for_s = 0;
    for (size_t c = 0; c < cluster_count && count_for_s < target_s; ++c)
    {
      float cx = parent_samples[0][c];
      float cy = parent_samples[1][c];

      size_t pts_to_sample = std::min(points_per_cluster,
                                      target_s - count_for_s);
      if (pts_to_sample == 0) break;

      glm::vec4 cluster_bbox = {
          cx - cluster_spread,
          cx + cluster_spread,
          cy - cluster_spread,
          cy + cluster_spread,
      };

      for (size_t k = 0; k < pts_to_sample; ++k)
      {
        glm::vec2 pt = sampler.sample(cluster_bbox);
        candidate_trees.emplace_back(pt.x,
                                     pt.y,
                                     0.0f,
                                     static_cast<uint32_t>(s),
                                     default_radius);
        count_for_s++;
      }
    }
  }

  Forest forest(std::move(candidate_trees));
  return forest;
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
