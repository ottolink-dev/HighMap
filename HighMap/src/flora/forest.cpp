/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>

#include "point_sampler/metrics.hpp"

#include "highmap/flora/forest.hpp"
#include "highmap/flora/species.hpp"

#include <unordered_map>

namespace hmap
{

// ==========================================================================
//  Constructors
// ==========================================================================

Forest::Forest(const std::vector<Tree> &trees) : ScatterField(trees)
{
}

Forest::Forest(std::vector<Tree> &&trees) noexcept
    : ScatterField(std::move(trees))
{
}

Forest::Forest(const ScatterField &field) : ScatterField(field)
{
}

Forest::Forest(ScatterField &&field) noexcept : ScatterField(std::move(field))
{
}

Forest::Forest(const Cloud &cloud, uint32_t species_id, float default_radius)
    : ScatterField(cloud, species_id, default_radius)
{
}

void Forest::prune_unviable(const std::vector<Species> &species,
                            bool                        prune_collisions)
{
  if (items.size() < 2) return;

  // --- Map Species Definitions for Fast Lookup

  std::unordered_map<uint32_t, Species> species_map;
  for (const auto &sp : species)
  {
    species_map[sp.id] = sp;
  }

  auto get_species_traits = [&](const Tree &tree) -> Species
  {
    size_t s_i = tree.class_id;
    auto   it = species_map.find(tree.class_id);
    if (it != species_map.end())
    {
      return it->second;
    }
    if (s_i < species.size())
    {
      return species[s_i];
    }
    return Species(tree.class_id, tree.radius);
  };

  // --- KD-Tree Nearest Neighbor Query

  std::vector<ps::Point<float, 2>> points;
  points.reserve(items.size());
  for (const auto &tree : items)
  {
    points.push_back({tree.position.x, tree.position.y});
  }

  auto neighbors_idx = ps::nearest_neighbors_indices(points, 1);

  // --- Check Minimum Radius / Viability

  std::vector<Tree> viable;
  viable.reserve(items.size());

  for (size_t i = 0; i < items.size(); ++i)
  {
    const Tree &tree = items[i];
    Species     sp_i = get_species_traits(tree);

    size_t neighbor_idx = neighbors_idx[i].empty() ? i : neighbors_idx[i][0];
    const Tree &neighbor_tree = items[neighbor_idx];

    float dx = tree.position.x - neighbor_tree.position.x;
    float dy = tree.position.y - neighbor_tree.position.y;
    float d_nn = std::sqrt(dx * dx + dy * dy);

    float r_est = sp_i.competition_factor * d_nn;

    // cull choked plants that cannot meet species radius_min
    if (r_est < sp_i.radius_min)
    {
      continue;
    }

    viable.push_back(tree);
  }

  items = std::move(viable);

  if (prune_collisions)
  {
    ScatterField::prune_collisions();
  }
}

void Forest::reinforce_species_clusters(size_t iterations,
                                        size_t k_neighbors,
                                        bool   include_self)
{
  if (items.size() < 2 || k_neighbors == 0 || iterations == 0) return;

  // --- Query Nearest Neighbors Graph Once (Static Positions)

  std::vector<ps::Point<float, 2>> points;
  points.reserve(items.size());
  for (const auto &tree : items)
  {
    points.push_back({tree.position.x, tree.position.y});
  }

  size_t safe_k = std::min(k_neighbors, items.size() - 1);
  auto   neighbors_idx = ps::nearest_neighbors_indices(points, safe_k);

  // --- Iterative Majority Species Assignment

  std::vector<uint32_t> next_species(items.size());

  for (size_t iter = 0; iter < iterations; ++iter)
  {
    for (size_t i = 0; i < items.size(); ++i)
    {
      std::unordered_map<uint32_t, size_t> counts;

      if (include_self)
      {
        counts[items[i].class_id]++;
      }

      for (size_t neighbor_idx : neighbors_idx[i])
      {
        counts[items[neighbor_idx].class_id]++;
      }

      // pick dominant species (majority vote)
      uint32_t best_species = items[i].class_id;
      size_t   max_count = 0;

      // check current tree species first to favor status quo in case of ties
      auto self_it = counts.find(items[i].class_id);
      if (self_it != counts.end())
      {
        best_species = self_it->first;
        max_count = self_it->second;
      }

      for (const auto &[sp_id, count] : counts)
      {
        if (count > max_count)
        {
          max_count = count;
          best_species = sp_id;
        }
      }

      next_species[i] = best_species;
    }

    for (size_t i = 0; i < items.size(); ++i)
    {
      items[i].class_id = next_species[i];
    }
  }
}

void Forest::shuffle_species(float ratio, size_t k_neighbors, uint32_t seed)
{
  if (items.size() < 2 || ratio <= 0.0f || k_neighbors == 0) return;

  // --- Build 2D Point List and Query Nearest Neighbors

  std::vector<ps::Point<float, 2>> points;
  points.reserve(items.size());
  for (const auto &tree : items)
  {
    points.push_back({tree.position.x, tree.position.y});
  }

  size_t safe_k = std::min(k_neighbors, items.size() - 1);
  auto   neighbors_idx = ps::nearest_neighbors_indices(points, safe_k);

  // --- Select Candidate Trees to Shuffle

  std::mt19937 gen(seed);

  std::vector<size_t> perm(items.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::shuffle(perm.begin(), perm.end(), gen);

  size_t target_count = std::min(
      items.size(),
      static_cast<size_t>(
          std::round(ratio * static_cast<float>(items.size()))));

  // --- Perform Neighbor Species Swaps with Failsafe

  for (size_t c = 0; c < target_count; ++c)
  {
    size_t i = perm[c];

    // find neighbors with a differing species
    std::vector<size_t> valid_neighbors;
    for (size_t neighbor_idx : neighbors_idx[i])
    {
      if (items[neighbor_idx].class_id != items[i].class_id)
      {
        valid_neighbors.push_back(neighbor_idx);
      }
    }

    // if all neighbors share the same species, skip candidate
    if (valid_neighbors.empty()) continue;

    // pick one differing neighbor at random and swap species & radius
    std::uniform_int_distribution<size_t> dis(0, valid_neighbors.size() - 1);
    size_t chosen_neighbor = valid_neighbors[dis(gen)];

    std::swap(items[i].class_id, items[chosen_neighbor].class_id);
    std::swap(items[i].radius, items[chosen_neighbor].radius);
  }
}

std::string Forest::to_string() const
{
  std::ostringstream ss;
  ss.imbue(std::locale("C"));

  glm::vec4 bbox = get_bbox();

  // count trees and collect average radius per species
  std::map<uint32_t, size_t> species_counts;
  std::map<uint32_t, float>  species_radii_sum;
  for (const auto &tree : items)
  {
    species_counts[tree.class_id]++;
    species_radii_sum[tree.class_id] += tree.radius;
  }

  ss << "Forest Infos\n";
  ss << "--------------------------------\n";
  ss << " total trees : " << items.size() << "\n";
  ss << " species count : " << species_counts.size() << "\n";
  ss << " bbox : {" << bbox.x << ", " << bbox.y << ", " << bbox.z << ", "
     << bbox.w << "}\n";

  if (!species_counts.empty())
  {
    ss << " species breakdown :\n";
    for (const auto &[sp_id, count] : species_counts)
    {
      float avg_r = (count > 0) ? (species_radii_sum[sp_id] / float(count))
                                : 0.f;
      float pct = (items.empty())
                      ? 0.f
                      : (100.f * float(count) / float(items.size()));
      ss << "   - species " << sp_id << ": " << count << " trees ("
         << std::fixed << std::setprecision(1) << pct << "%, avg radius "
         << std::setprecision(4) << avg_r << ")\n";
    }
  }
  ss << "--------------------------------";

  return ss.str();
}

} // namespace hmap
