/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "highmap/flora/forest_growth.hpp"
#include "highmap/internal/flora_utils.hpp"
#include "highmap/math.hpp"
#include "highmap/terrain_tri_mesh.hpp"

namespace hmap
{

Forest grow_forest_competition_voronoi(
    const Forest               &forest,
    const std::vector<Species> &species,
    const InteractionMatrix    &competition_matrix,
    const Array                &max_radius_scale,
    float                       max_radius_scale_strength,
    bool                        prune_unviable,
    const glm::vec4            &bbox)
{
  if (forest.empty()) return Forest();

  ForestScaleSampler sampler(max_radius_scale, max_radius_scale_strength, bbox);
  SpeciesLookup      lookup(species);

  // if fewer than 3 trees, triangulate cannot form triangles -> fallback to
  // species clamping
  if (forest.size() < 3)
  {
    std::vector<Tree> result;
    result.reserve(forest.size());
    for (const auto &t : forest)
    {
      Species sp_i = lookup.get(t);
      float   s = sampler.sample(t.position.x, t.position.y);
      float   r_max_eff = lerp(sp_i.radius_min, sp_i.radius_max, s);
      Tree    grown_tree = t;
      grown_tree.radius = std::clamp(t.radius, sp_i.radius_min, r_max_eff);
      result.push_back(grown_tree);
    }
    return Forest(std::move(result));
  }

  // --- Build 3D Points for Delaunay Triangulation (TerrainTriMesh)

  std::vector<glm::vec3> points;
  points.reserve(forest.size());
  for (const auto &tree : forest)
  {
    points.push_back(glm::vec3(tree.position.x, tree.position.y, 0.0f));
  }

  TerrainTriMesh     mesh(points);
  std::vector<float> vertex_areas = mesh.get_vertex_areas(false);

  // If mesh failed to produce any triangles (e.g. collinear points)
  if (mesh.get_triangles().empty() || vertex_areas.size() != forest.size())
  {
    std::vector<Tree> result;
    result.reserve(forest.size());
    for (const auto &t : forest)
    {
      Species sp_i = lookup.get(t);
      float   s = sampler.sample(t.position.x, t.position.y);
      float   r_max_eff = lerp(sp_i.radius_min, sp_i.radius_max, s);
      Tree    grown_tree = t;
      grown_tree.radius = std::clamp(t.radius, sp_i.radius_min, r_max_eff);
      result.push_back(grown_tree);
    }
    return Forest(std::move(result));
  }

  // --- Compute Potential Growing Space & Radii

  std::vector<Tree> result;
  result.reserve(forest.size());

  const float inv_pi = 1.0f / static_cast<float>(M_PI);

  for (size_t i = 0; i < forest.size(); ++i)
  {
    const Tree &tree = forest[i];
    Species     sp_i = lookup.get(tree);

    // determine species competition factor alpha
    float alpha = lookup.get_alpha(sp_i, sp_i.id, competition_matrix);

    // territory area and unconstrained equivalent radius
    float area_i = std::max(0.0f, vertex_areas[i]);
    float r_raw = std::sqrt(area_i * inv_pi);

    // scale max_radius locally if array is provided
    float s = sampler.sample(tree.position.x, tree.position.y);
    float r_max_eff = lerp(sp_i.radius_min, sp_i.radius_max, s);

    // estimated crown radius
    float r_est = alpha * r_raw;

    // check viability
    if (prune_unviable && r_est < sp_i.radius_min)
    {
      continue; // cull choked plant
    }

    Tree grown_tree = tree;
    grown_tree.radius = std::clamp(r_est, sp_i.radius_min, r_max_eff);
    result.push_back(grown_tree);
  }

  Forest grown_forest(std::move(result));

  if (prune_unviable) grown_forest.prune_unviable(species, true);

  return grown_forest;
}

} // namespace hmap
