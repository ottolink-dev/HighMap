/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <locale>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "point_sampler/metrics.hpp"
#include "point_sampler/point.hpp"

#include "highmap/flora/forest.hpp"
#include "highmap/flora/species.hpp"
#include "highmap/flora/tree.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/scatter/scatter_field.hpp"
#include "highmap/scatter/scatter_item.hpp"

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
