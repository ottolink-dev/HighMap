/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/flora/forest_growth.hpp"
#include "highmap/flora/species.hpp"
#include "highmap/flora/tree.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math.hpp"

#include <unordered_map>

namespace hmap
{

/**
 * @struct ForestScaleSampler
 * @brief Helper for sampling spatial scaling maps over a 2D bounding box with
 * blending.
 */
struct ForestScaleSampler
{
  const Array *scale_array = nullptr;
  glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f};
  float        strength = 1.0f;
  float        bbox_dx = 1.0f;
  float        bbox_dy = 1.0f;
  bool         active = false;

  ForestScaleSampler(const Array     &array,
                     float            strength_factor,
                     const glm::vec4 &bounds)
      : scale_array(&array),
        bbox(bounds),
        strength(std::clamp(strength_factor, 0.0f, 1.0f)),
        bbox_dx(bounds.y - bounds.x),
        bbox_dy(bounds.w - bounds.z)
  {
    bool has_data = !array.vector.empty() && validate_non_empty(array);
    bool valid_bbox = (std::abs(bbox_dx) > 1e-7f && std::abs(bbox_dy) > 1e-7f);
    active = (has_data && valid_bbox && strength > 0.0f);
  }

  [[nodiscard]] float sample(float x, float y) const
  {
    if (!active || !scale_array) return 1.0f;

    float u_norm = std::clamp((x - bbox.x) / bbox_dx, 0.0f, 1.0f);
    float v_norm = std::clamp((y - bbox.z) / bbox_dy, 0.0f, 1.0f);

    float xn = u_norm * static_cast<float>(scale_array->shape.x - 1);
    float yn = v_norm * static_cast<float>(scale_array->shape.y - 1);

    int i = static_cast<int>(xn);
    int j = static_cast<int>(yn);

    float u = xn - static_cast<float>(i);
    float v = yn - static_cast<float>(j);
    float s_raw = std::clamp(scale_array->get_value_bilinear_at(i, j, u, v),
                             0.0f,
                             1.0f);
    return lerp(1.f, s_raw, strength);
  }
};

/**
 * @struct SpeciesLookup
 * @brief Helper for querying species traits and pairwise competition factors.
 */
struct SpeciesLookup
{
  std::unordered_map<uint32_t, Species> map;
  const std::vector<Species>           *list = nullptr;

  explicit SpeciesLookup(const std::vector<Species> &species) : list(&species)
  {
    for (const auto &sp : species)
    {
      map[sp.id] = sp;
    }
  }

  [[nodiscard]] Species get(const Tree &tree) const
  {
    auto it = map.find(tree.class_id);
    if (it != map.end())
    {
      return it->second;
    }
    if (list && tree.class_id < list->size())
    {
      return (*list)[tree.class_id];
    }
    return Species(tree.class_id, tree.radius);
  }

  [[nodiscard]] float get_alpha(const Species           &sp_i,
                                size_t                   s_j,
                                const InteractionMatrix &matrix) const
  {
    if (matrix.size > 0)
    {
      float mat_val = matrix.get(sp_i.id, s_j);
      if (mat_val > 0.0f) return mat_val;
    }
    return sp_i.competition_factor;
  }
};

} // namespace hmap
