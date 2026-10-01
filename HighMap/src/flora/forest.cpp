/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

#include "highmap/flora/forest.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

// ==========================================================================
//  Constructors
// ==========================================================================

Forest::Forest(const std::vector<Tree> &trees) : trees(trees)
{
}

Forest::Forest(std::vector<Tree> &&trees) noexcept : trees(std::move(trees))
{
}

Forest::Forest(const Cloud &cloud, uint32_t species_id, float default_radius)
{
  trees.reserve(cloud.size());
  for (const auto &p : cloud)
  {
    float r = (std::abs(p.v) > 1e-6f) ? p.v : default_radius;
    trees.emplace_back(p.x, p.y, 0.0f, species_id, r);
  }
}

// ==========================================================================
//  Container Accessors
// ==========================================================================

Tree &Forest::at(size_t index)
{
  if (index >= trees.size())
    throw std::out_of_range("Forest::at index out of range");
  return trees.at(index);
}

const Tree &Forest::at(size_t index) const
{
  if (index >= trees.size())
    throw std::out_of_range("Forest::at index out of range");
  return trees.at(index);
}

Tree &Forest::back()
{
  return trees.back();
}

const Tree &Forest::back() const
{
  return trees.back();
}

Tree &Forest::front()
{
  return trees.front();
}

const Tree &Forest::front() const
{
  return trees.front();
}

// ==========================================================================
//  Operations
// ==========================================================================

Forest Forest::filter_by_species(uint32_t species_id) const
{
  Forest filtered;
  for (const auto &tree : trees)
  {
    if (tree.species_id == species_id) filtered.push_back(tree);
  }
  return filtered;
}

glm::vec4 Forest::get_bbox() const
{
  if (trees.empty()) return {0.f, 1.f, 0.f, 1.f};

  float xmin = trees[0].position.x;
  float xmax = trees[0].position.x;
  float ymin = trees[0].position.y;
  float ymax = trees[0].position.y;

  for (const auto &tree : trees)
  {
    xmin = std::min(xmin, tree.position.x);
    xmax = std::max(xmax, tree.position.x);
    ymin = std::min(ymin, tree.position.y);
    ymax = std::max(ymax, tree.position.y);
  }

  return {xmin, xmax, ymin, ymax};
}

std::vector<uint32_t> Forest::get_species_ids() const
{
  std::set<uint32_t> unique_species;
  for (const auto &tree : trees)
    unique_species.insert(tree.species_id);

  return std::vector<uint32_t>(unique_species.begin(), unique_species.end());
}

void Forest::set_elevation_from_terrain(const Array     &elevation,
                                        const glm::vec4 &bbox)
{
  if (!validate_non_empty(elevation)) return;

  float dx = bbox.y - bbox.x;
  float dy = bbox.w - bbox.z;
  if (std::abs(dx) < 1e-7f || std::abs(dy) < 1e-7f) return;

  for (auto &tree : trees)
  {
    // scale to unit interval
    float xn = (tree.position.x - bbox.x) / dx;
    float yn = (tree.position.y - bbox.z) / dy;

    // scale to array shape
    xn *= static_cast<float>(elevation.shape.x - 1);
    yn *= static_cast<float>(elevation.shape.y - 1);

    int i = static_cast<int>(xn);
    int j = static_cast<int>(yn);

    // sample only within bounds
    if (i >= 0 && i < elevation.shape.x && j >= 0 && j < elevation.shape.y)
    {
      float u = xn - static_cast<float>(i);
      float v = yn - static_cast<float>(j);
      tree.position.z = elevation.get_value_bilinear_at(i, j, u, v);
    }
  }
}

Cloud Forest::to_cloud() const
{
  std::vector<Point> pts;
  pts.reserve(trees.size());
  for (const auto &tree : trees)
    pts.push_back(tree.to_point());
  return Cloud(std::move(pts));
}

} // namespace hmap
