/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "highmap/flora/forest_growth.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math.hpp"
#include "highmap/terrain_tri_mesh.hpp"

#include <unordered_map>

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

  // --- Map Species Definitions for Fast Lookup

  std::unordered_map<uint32_t, Species> species_map;
  for (const auto &sp : species)
  {
    species_map[sp.id] = sp;
  }

  auto get_species_traits = [&](const Tree &tree) -> Species
  {
    size_t s_i = tree.species_id;
    auto   it = species_map.find(tree.species_id);
    if (it != species_map.end())
    {
      return it->second;
    }
    if (s_i < species.size())
    {
      return species[s_i];
    }
    return Species(tree.species_id, tree.radius);
  };

  // if fewer than 3 trees, triangulate cannot form triangles -> fallback to
  // species clamping
  if (forest.size() < 3)
  {
    std::vector<Tree> result;
    result.reserve(forest.size());
    for (const auto &t : forest)
    {
      Species sp_i = get_species_traits(t);
      float   s = sample_scale(t.position.x, t.position.y);
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
      Species sp_i = get_species_traits(t);
      float   s = sample_scale(t.position.x, t.position.y);
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
    Species     sp_i = get_species_traits(tree);
    size_t      s_i = tree.species_id;

    // determine species competition factor alpha
    float alpha = sp_i.competition_factor;
    if (competition_matrix.size > 0)
    {
      float mat_val = competition_matrix.get(s_i, s_i);
      if (mat_val > 0.0f) alpha = mat_val;
    }

    // territory area and unconstrained equivalent radius
    float area_i = std::max(0.0f, vertex_areas[i]);
    float r_raw = std::sqrt(area_i * inv_pi);

    // scale max_radius locally if array is provided
    float s = sample_scale(tree.position.x, tree.position.y);
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
