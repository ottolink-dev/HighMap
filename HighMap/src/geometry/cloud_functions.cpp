/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/kd_tree.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

Array cloud_sdf_to_array(const Cloud &cloud,
                         glm::ivec2   shape,
                         glm::vec4    bbox_array,
                         const Array *p_noise_x,
                         const Array *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  Array array(shape);

  if (!validate_min_size(cloud, 1, "Cloud points")) return array;

  // --- KD-tree

  std::vector<float> x = cloud.get_x();
  std::vector<float> y = cloud.get_y();
  KDTreeContext      tree(x, y);

  // --- SDF

  // array base grid
  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, bbox_array, /* endpoint */ false);

#pragma omp parallel for default(none)                                         \
    shared(shape, p_noise_x, p_noise_y, xg, yg, tree, array)
  for (int j = 0; j < shape.y; ++j)
  {
    std::vector<size_t> indices;
    std::vector<float>  distances;

    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;
      float xi = xg[i] + dx;
      float yi = yg[j] + dy;

      tree.neighbor_search(xi,
                           yi,
                           /* k_neighbors */ 1,
                           indices,
                           distances);

      array(i, j) = std::sqrt(distances[0]);
    }
  }

  return array;
}

bool has_duplicates(const Cloud &cloud, float eps, bool xy_only)
{
  std::vector<glm::vec3> pts = cloud.to_vec3();

  std::sort(pts.begin(),
            pts.end(),
            [](const auto &a, const auto &b)
            { return a.x < b.x || (a.x == b.x && a.y < b.y); });

  if (xy_only)
  {
    for (size_t i = 1; i < pts.size(); ++i)
    {
      glm::vec2 p0 = {pts[i].x, pts[i].y};
      glm::vec2 p1 = {pts[i - 1].x, pts[i - 1].y};

      if (glm::distance(p0, p1) < eps) return true;
    }
  }
  else
  {
    for (size_t i = 1; i < pts.size(); ++i)
      if (glm::distance(pts[i], pts[i - 1]) < eps) return true;
  }

  return false;
}

std::vector<float> interpolate_values_from_array(const Cloud &cloud,
                                                 const Array &array,
                                                 glm::vec4    bbox)
{
  if (!validate_non_empty(array)) return std::vector<float>(cloud.size(), 0.f);

  const float      inv_width = 1.0f / (bbox.y - bbox.x);
  const float      inv_height = 1.0f / (bbox.w - bbox.z);
  const glm::ivec2 shape = {array.shape.x - 1, array.shape.y - 1};

  std::vector<float> values;
  values.reserve(cloud.size());

  for (const auto &p : cloud)
  {
    const float xn = (p.x - bbox.x) * inv_width;
    const float yn = (p.y - bbox.z) * inv_height;

    if (xn < 0.0f || xn > 1.0f || yn < 0.0f || yn > 1.0f)
    {
      values.push_back(0.0f);
      continue;
    }

    const float x_scaled = xn * shape.x;
    const float y_scaled = yn * shape.y;
    const int   i = static_cast<int>(x_scaled);
    const int   j = static_cast<int>(y_scaled);

    if (i >= 0 && i < array.shape.x && j >= 0 && j < array.shape.y)
    {
      const float uu = x_scaled - i;
      const float vv = y_scaled - j;
      values.push_back(array.get_value_bilinear_at(i, j, uu, vv));
    }
    else
    {
      values.push_back(0.0f);
    }
  }
  return values;
}

void rejection_filter_density(Cloud           &cloud,
                              const Array     &density_mask,
                              std::uint32_t    seed,
                              const glm::vec4 &bbox)
{
  if (!validate_non_empty(density_mask)) return;

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(0.f, 1.f);

  auto density_fct = make_xy_function_from_array(density_mask, bbox);

  cloud.erase(std::remove_if(cloud.begin(),
                             cloud.end(),
                             [&](const Point &p)
                             {
                               float rnd = dis(gen);
                               return (rnd > density_fct(p.x, p.y));
                             }),
              cloud.end());
}

Cloud scale(const Cloud &cloud, glm::vec2 scale, glm::vec2 center)
{
  Cloud result = cloud;

  for (auto &p : result)
    p = hmap::scale(p, scale, center);

  return result;
}

Cloud scale(const Cloud &cloud, float scale_factor, glm::vec2 center)
{
  return scale(cloud, glm::vec2(scale_factor, scale_factor), center);
}

} // namespace hmap
