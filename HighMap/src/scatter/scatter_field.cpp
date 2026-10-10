/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <locale>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/core/types.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include "cl_wrapper/device_manager.hpp"
#include "cl_wrapper/run.hpp"
#include "delaunator-cpp.hpp"
#include "point_sampler/metrics.hpp"
#include "point_sampler/point.hpp"
#include "point_sampler/relaxation.hpp"

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/range.hpp"
#include "highmap/scatter/scatter_field.hpp"
#include "highmap/scatter/scatter_item.hpp"

#include <format>
#include <unordered_map>

namespace hmap
{

// --- Helper Functions

static cv::Scalar get_scatter_class_color(uint32_t class_id)
{
  static const std::vector<cv::Scalar> palette = {
      cv::Scalar(60, 200, 60),   // vibrant green
      cv::Scalar(40, 140, 245),  // warm orange
      cv::Scalar(235, 205, 50),  // cyan / teal
      cv::Scalar(210, 60, 200),  // magenta
      cv::Scalar(50, 230, 255),  // yellow
      cv::Scalar(240, 90, 70),   // blue
      cv::Scalar(90, 80, 240),   // coral / red
      cv::Scalar(180, 230, 100), // lime
      cv::Scalar(200, 130, 250), // pink / lilac
  };

  if (class_id < palette.size()) return palette[class_id];

  // golden ratio hue rotation for arbitrary class counts
  float   hue = std::fmod(static_cast<float>(class_id) * 137.508f, 180.0f);
  cv::Mat hsv(1, 1, CV_8UC3, cv::Scalar(static_cast<uint8_t>(hue), 220, 240));
  cv::Mat bgr;
  cv::cvtColor(hsv, bgr, cv::COLOR_HSV2BGR);
  cv::Vec3b vec = bgr.at<cv::Vec3b>(0, 0);
  return cv::Scalar(vec[0], vec[1], vec[2]);
}

// ==========================================================================
//  Constructors
// ==========================================================================

ScatterField::ScatterField(const std::vector<ScatterItem> &items) : items(items)
{
}

ScatterField::ScatterField(std::vector<ScatterItem> &&items) noexcept
    : items(std::move(items))
{
}

ScatterField::ScatterField(const Cloud &cloud,
                           uint32_t     class_id,
                           float        default_radius)
{
  items.reserve(cloud.size());
  for (const auto &p : cloud)
  {
    float r = (std::abs(p.v) > 1e-6f) ? p.v : default_radius;
    items.emplace_back(p.x, p.y, 0.0f, class_id, r);
  }
}

// ==========================================================================
//  Container Accessors
// ==========================================================================

ScatterItem &ScatterField::at(size_t index)
{
  if (index >= items.size())
    throw std::out_of_range("ScatterField::at index out of range");
  return items.at(index);
}

const ScatterItem &ScatterField::at(size_t index) const
{
  if (index >= items.size())
    throw std::out_of_range("ScatterField::at index out of range");
  return items.at(index);
}

ScatterItem &ScatterField::back()
{
  return items.back();
}

const ScatterItem &ScatterField::back() const
{
  return items.back();
}

ScatterItem &ScatterField::front()
{
  return items.front();
}

const ScatterItem &ScatterField::front() const
{
  return items.front();
}

// ==========================================================================
//  Operations
// ==========================================================================

void ScatterField::densify(float default_radius)
{
  if (items.size() < 3) return;

  // --- Prepare 2D Coordinates for Delaunay Triangulation

  std::vector<double> coords;
  coords.reserve(2 * items.size());
  for (const auto &item : items)
  {
    coords.push_back(static_cast<double>(item.position.x));
    coords.push_back(static_cast<double>(item.position.y));
  }

  delaunator::Delaunator d(coords);
  const auto            &tri = d.triangles;
  if (tri.empty()) return;

  // --- Compute Circumcenters & Assign Majority Class

  std::vector<ScatterItem> new_items;
  new_items.reserve(tri.size() / 3);

  for (size_t k = 0; k < tri.size(); k += 3)
  {
    size_t i0 = tri[k];
    size_t i1 = tri[k + 1];
    size_t i2 = tri[k + 2];

    double ax = coords[2 * i0];
    double ay = coords[2 * i0 + 1];
    double bx = coords[2 * i1];
    double by = coords[2 * i1 + 1];
    double cx = coords[2 * i2];
    double cy = coords[2 * i2 + 1];

    auto [cx_center,
          cy_center] = delaunator::circumcenter(ax, ay, bx, by, cx, cy);

    if (!std::isfinite(cx_center) || !std::isfinite(cy_center)) continue;

    // determine majority class id among the 3 triangle vertices
    uint32_t s0 = items[i0].class_id;
    uint32_t s1 = items[i1].class_id;
    uint32_t s2 = items[i2].class_id;

    uint32_t majority_class = s0;
    if (s0 == s1 || s0 == s2)
    {
      majority_class = s0;
    }
    else if (s1 == s2)
    {
      majority_class = s1;
    }
    else
    {
      majority_class = s0;
    }

    new_items.emplace_back(static_cast<float>(cx_center),
                           static_cast<float>(cy_center),
                           0.0f,
                           majority_class,
                           default_radius);
  }

  // --- Append New Items to ScatterField

  items.insert(items.end(),
               std::make_move_iterator(new_items.begin()),
               std::make_move_iterator(new_items.end()));
}

ScatterField ScatterField::filter_by_class(uint32_t class_id) const
{
  ScatterField filtered;
  for (const auto &item : items)
  {
    if (item.class_id == class_id) filtered.push_back(item);
  }
  return filtered;
}

glm::vec4 ScatterField::get_bbox() const
{
  if (items.empty()) return {0.f, 1.f, 0.f, 1.f};

  float xmin = items[0].position.x;
  float xmax = items[0].position.x;
  float ymin = items[0].position.y;
  float ymax = items[0].position.y;

  for (const auto &item : items)
  {
    xmin = std::min(xmin, item.position.x);
    xmax = std::max(xmax, item.position.x);
    ymin = std::min(ymin, item.position.y);
    ymax = std::max(ymax, item.position.y);
  }

  return {xmin, xmax, ymin, ymax};
}

std::vector<uint32_t> ScatterField::get_class_ids() const
{
  std::set<uint32_t> unique_classes;
  for (const auto &item : items)
    unique_classes.insert(item.class_id);

  return std::vector<uint32_t>(unique_classes.begin(), unique_classes.end());
}

std::vector<float> ScatterField::get_radius() const
{
  std::vector<float> r;
  r.reserve(items.size());
  for (const auto &item : items)
    r.push_back(item.radius);
  return r;
}

std::vector<float> ScatterField::get_x() const
{
  std::vector<float> x;
  x.reserve(items.size());
  for (const auto &item : items)
    x.push_back(item.position.x);
  return x;
}

std::vector<float> ScatterField::get_y() const
{
  std::vector<float> y;
  y.reserve(items.size());
  for (const auto &item : items)
    y.push_back(item.position.y);
  return y;
}

std::vector<float> ScatterField::get_z() const
{
  std::vector<float> z;
  z.reserve(items.size());
  for (const auto &item : items)
    z.push_back(item.position.z);
  return z;
}

void ScatterField::perturb_positions(float dx, float dy, uint32_t seed)
{
  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(-1.0f, 1.0f);

  for (auto &item : items)
  {
    item.position.x += dx * dis(gen);
    item.position.y += dy * dis(gen);
  }
}

void ScatterField::prune_collisions()
{
  if (items.size() < 2) return;

  // sort candidates by radius descending
  std::vector<size_t> order(items.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(),
                   order.end(),
                   [&](size_t a, size_t b)
                   { return items[a].radius > items[b].radius; });

  float max_r = 0.0f;
  for (const auto &t : items)
  {
    max_r = std::max(max_r, t.radius);
  }

  float                                            cell_size = 2.0f * max_r;
  std::unordered_map<int64_t, std::vector<size_t>> grid;

  auto compute_cell = [&](float x, float y) -> std::pair<int, int>
  {
    if (cell_size <= 1e-7f) return {0, 0};
    int gx = static_cast<int>(std::floor(x / cell_size));
    int gy = static_cast<int>(std::floor(y / cell_size));
    return {gx, gy};
  };

  auto make_key = [](int gx, int gy) -> int64_t
  {
    return (static_cast<int64_t>(gx) << 32) ^
           (static_cast<int64_t>(gy) & 0xFFFFFFFF);
  };

  std::vector<ScatterItem> retained;
  retained.reserve(items.size());

  for (size_t idx : order)
  {
    const ScatterItem &cand = items[idx];
    auto [gx, gy] = compute_cell(cand.position.x, cand.position.y);

    bool collides = false;

    for (int dy = -1; dy <= 1 && !collides; ++dy)
    {
      for (int dx = -1; dx <= 1 && !collides; ++dx)
      {
        int64_t key = make_key(gx + dx, gy + dy);
        auto    it = grid.find(key);
        if (it == grid.end()) continue;

        for (size_t ret_idx : it->second)
        {
          const ScatterItem &other = retained[ret_idx];
          float              dist_x = cand.position.x - other.position.x;
          float              dist_y = cand.position.y - other.position.y;
          float              dist_sq = dist_x * dist_x + dist_y * dist_y;

          float min_dist = cand.radius + other.radius;
          if (dist_sq < min_dist * min_dist)
          {
            collides = true;
            break;
          }
        }
      }
    }

    if (!collides)
    {
      retained.push_back(cand);
      grid[make_key(gx, gy)].push_back(retained.size() - 1);
    }
  }

  items = std::move(retained);
}

void ScatterField::prune_collisions(const ScatterField &other)
{
  if (items.empty() || other.empty()) return;

  // find max radius across both fields to dimension the spatial grid
  float max_r_other = 0.0f;
  for (const auto &item : other)
  {
    max_r_other = std::max(max_r_other, item.radius);
  }

  float max_r_self = 0.0f;
  for (const auto &item : items)
  {
    max_r_self = std::max(max_r_self, item.radius);
  }

  float cell_size = max_r_other + max_r_self;
  if (cell_size <= 1e-7f) return;

  auto compute_cell = [&](float x, float y) -> std::pair<int, int>
  {
    int gx = static_cast<int>(std::floor(x / cell_size));
    int gy = static_cast<int>(std::floor(y / cell_size));
    return {gx, gy};
  };

  auto make_key = [](int gx, int gy) -> int64_t
  {
    return (static_cast<int64_t>(gx) << 32) ^
           (static_cast<int64_t>(gy) & 0xFFFFFFFF);
  };

  // populate spatial hash grid with reference items from `other`
  std::unordered_map<int64_t, std::vector<size_t>> grid;
  grid.reserve(other.size());

  for (size_t i = 0; i < other.size(); ++i)
  {
    auto [gx, gy] = compute_cell(other[i].position.x, other[i].position.y);
    grid[make_key(gx, gy)].push_back(i);
  }

  std::vector<ScatterItem> retained;
  retained.reserve(items.size());

  for (const auto &item : items)
  {
    auto [gx, gy] = compute_cell(item.position.x, item.position.y);
    bool collides = false;

    for (int dy = -1; dy <= 1 && !collides; ++dy)
    {
      for (int dx = -1; dx <= 1 && !collides; ++dx)
      {
        int64_t key = make_key(gx + dx, gy + dy);
        auto    it = grid.find(key);
        if (it == grid.end()) continue;

        for (size_t other_idx : it->second)
        {
          const ScatterItem &ref_item = other[other_idx];
          float              dist_x = item.position.x - ref_item.position.x;
          float              dist_y = item.position.y - ref_item.position.y;
          float              dist_sq = dist_x * dist_x + dist_y * dist_y;

          float min_dist = item.radius + ref_item.radius;
          if (dist_sq < min_dist * min_dist)
          {
            collides = true;
            break;
          }
        }
      }
    }

    if (!collides)
    {
      retained.push_back(item);
    }
  }

  items = std::move(retained);
}

void ScatterField::prune_density(const Array     &density_mask,
                                 float            target_ratio,
                                 uint32_t         seed,
                                 const glm::vec4 &bbox)
{
  if (!validate_non_empty(density_mask) || items.empty()) return;

  float r_target = std::clamp(target_ratio, 0.0f, 1.0f);
  if (r_target <= 1e-6f)
  {
    items.clear();
    return;
  }
  if (r_target >= 1.0f - 1e-6f)
  {
    return;
  }

  // --- Sample Local Density at Item Positions

  auto density_fct = make_xy_function_from_array(density_mask, bbox);
  std::vector<float> sampled_densities(items.size());
  for (size_t i = 0; i < items.size(); ++i)
  {
    sampled_densities[i] = std::max(
        0.0f,
        density_fct(items[i].position.x, items[i].position.y));
  }

  // --- Remap Acceptance Probabilities to Match Target Keep Ratio

  float low = -1.0f;
  float high = 1.0f;
  float c_best = 0.0f;

  for (int iter = 0; iter < 40; ++iter)
  {
    float mid = 0.5f * (low + high);
    float sum_p = 0.0f;
    for (float d : sampled_densities)
    {
      sum_p += std::clamp(d + mid, 0.0f, 1.0f);
    }
    float expected_ratio = sum_p / static_cast<float>(items.size());

    if (expected_ratio < r_target)
    {
      low = mid;
    }
    else
    {
      high = mid;
    }
    c_best = mid;
  }

  // --- Perform Density-Modulated Pruning

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(0.0f, 1.0f);

  std::vector<ScatterItem> retained;
  retained.reserve(
      static_cast<size_t>(std::ceil(r_target * float(items.size()))));

  for (size_t i = 0; i < items.size(); ++i)
  {
    float p_accept = std::clamp(sampled_densities[i] + c_best, 0.0f, 1.0f);
    if (dis(gen) <= p_accept)
    {
      retained.push_back(items[i]);
    }
  }

  items = std::move(retained);
}

void ScatterField::regularize_positions(size_t k_neighbors,
                                        float  step_size,
                                        size_t iterations)
{
  if (items.size() < 2 || k_neighbors == 0 || iterations == 0) return;

  std::vector<ps::Point<float, 2>> points;
  points.reserve(items.size());
  for (const auto &item : items)
  {
    points.push_back({item.position.x, item.position.y});
  }

  ps::relaxation_ktree<float, 2>(points, k_neighbors, step_size, iterations);

  for (size_t i = 0; i < items.size(); ++i)
  {
    items[i].position.x = points[i][0];
    items[i].position.y = points[i][1];
  }
}

void ScatterField::reinforce_class_clusters(size_t iterations,
                                            size_t k_neighbors,
                                            bool   include_self)
{
  if (items.size() < 2 || k_neighbors == 0 || iterations == 0) return;

  // --- Query Nearest Neighbors Graph Once (Static Positions)

  std::vector<ps::Point<float, 2>> points;
  points.reserve(items.size());
  for (const auto &item : items)
  {
    points.push_back({item.position.x, item.position.y});
  }

  size_t safe_k = std::min(k_neighbors, items.size() - 1);
  auto   neighbors_idx = ps::nearest_neighbors_indices(points, safe_k);

  // --- Iterative Majority Class Assignment

  std::vector<uint32_t> next_classes(items.size());

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

      // pick dominant class (majority vote)
      uint32_t best_class = items[i].class_id;
      size_t   max_count = 0;

      // check current item class first to favor status quo in case of ties
      auto self_it = counts.find(items[i].class_id);
      if (self_it != counts.end())
      {
        best_class = self_it->first;
        max_count = self_it->second;
      }

      for (const auto &[cls_id, count] : counts)
      {
        if (count > max_count)
        {
          max_count = count;
          best_class = cls_id;
        }
      }

      next_classes[i] = best_class;
    }

    for (size_t i = 0; i < items.size(); ++i)
    {
      items[i].class_id = next_classes[i];
    }
  }
}

void ScatterField::resolve_collisions(size_t iterations,
                                      float  tolerance,
                                      float  step_size,
                                      size_t triangulation_substep)
{
  if (items.size() < 2 || iterations == 0) return;

  float  overlap_factor = 1.0f - std::clamp(tolerance, 0.0f, 1.0f);
  size_t step_freq = std::max<size_t>(1, triangulation_substep);

  std::vector<std::pair<size_t, size_t>> edges;

  for (size_t iter = 0; iter < iterations; ++iter)
  {
    // --- Update Delaunay Triangulation Neighborhood
    if (iter % step_freq == 0 || edges.empty())
    {
      edges.clear();
      if (items.size() >= 3)
      {
        std::vector<double> coords;
        coords.reserve(2 * items.size());
        for (const auto &item : items)
        {
          coords.push_back(static_cast<double>(item.position.x));
          coords.push_back(static_cast<double>(item.position.y));
        }

        delaunator::Delaunator d(coords);
        const auto            &tri = d.triangles;
        const auto            &halfedges = d.halfedges;

        for (size_t e = 0; e < tri.size(); ++e)
        {
          int opposite = static_cast<int>(halfedges[e]);
          if (static_cast<int>(e) > opposite || opposite == -1)
          {
            size_t next_e = (e % 3 == 2) ? e - 2 : e + 1;
            size_t u = tri[e];
            size_t v = tri[next_e];
            edges.emplace_back(u, v);
          }
        }
      }
      else
      {
        edges.emplace_back(0, 1);
      }
    }

    bool had_collision = false;

    // --- Displace In Place Edge by Edge
    for (const auto &[i, j] : edges)
    {
      float min_dist = (items[i].radius + items[j].radius) * overlap_factor;
      float dx = items[j].position.x - items[i].position.x;
      float dy = items[j].position.y - items[i].position.y;
      float dist_sq = dx * dx + dy * dy;

      if (dist_sq < min_dist * min_dist)
      {
        had_collision = true;
        float     dist = std::sqrt(dist_sq);
        glm::vec2 dir;

        if (dist > 1e-6f)
        {
          dir = glm::vec2(dx / dist, dy / dist);
        }
        else
        {
          // break symmetry deterministically if points are coincident
          float angle = static_cast<float>((i + j + iter) % 360) *
                        (3.14159265f / 180.0f);
          dir = glm::vec2(std::cos(angle), std::sin(angle));
          dist = 0.0f;
        }

        float     penetration = min_dist - dist;
        glm::vec2 delta = 0.5f * step_size * penetration * dir;

        items[i].position.x -= delta.x;
        items[i].position.y -= delta.y;
        items[j].position.x += delta.x;
        items[j].position.y += delta.y;
      }
    }

    if (!had_collision) break;
  }
}

void ScatterField::resolve_collisions(const ScatterField &other,
                                      size_t              iterations,
                                      float               tolerance,
                                      float               step_size,
                                      size_t              triangulation_substep)
{
  if (items.empty() || other.empty() || iterations == 0) return;

  float  overlap_factor = 1.0f - std::clamp(tolerance, 0.0f, 1.0f);
  size_t step_freq = std::max<size_t>(1, triangulation_substep);

  size_t n_self = items.size();
  size_t n_other = other.size();
  size_t total_n = n_self + n_other;

  std::vector<std::pair<size_t, size_t>> cross_edges;

  for (size_t iter = 0; iter < iterations; ++iter)
  {
    // --- Update Joint Delaunay Triangulation Neighborhood
    if (iter % step_freq == 0 || cross_edges.empty())
    {
      cross_edges.clear();
      if (total_n >= 3)
      {
        std::vector<double> coords;
        coords.reserve(2 * total_n);
        for (const auto &item : items)
        {
          coords.push_back(static_cast<double>(item.position.x));
          coords.push_back(static_cast<double>(item.position.y));
        }
        for (const auto &item : other)
        {
          coords.push_back(static_cast<double>(item.position.x));
          coords.push_back(static_cast<double>(item.position.y));
        }

        delaunator::Delaunator d(coords);
        const auto            &tri = d.triangles;
        const auto            &halfedges = d.halfedges;

        for (size_t e = 0; e < tri.size(); ++e)
        {
          int opposite = static_cast<int>(halfedges[e]);
          if (static_cast<int>(e) > opposite || opposite == -1)
          {
            size_t next_e = (e % 3 == 2) ? e - 2 : e + 1;
            size_t u = tri[e];
            size_t v = tri[next_e];

            // identify edges connecting self (< n_self) and other (>= n_self)
            if (u < n_self && v >= n_self)
            {
              cross_edges.emplace_back(u, v - n_self);
            }
            else if (v < n_self && u >= n_self)
            {
              cross_edges.emplace_back(v, u - n_self);
            }
          }
        }
      }
      else
      {
        for (size_t i = 0; i < n_self; ++i)
        {
          for (size_t j = 0; j < n_other; ++j)
          {
            cross_edges.emplace_back(i, j);
          }
        }
      }
    }

    bool had_collision = false;

    // --- Displace In Place Edge by Edge from Static Other Items
    for (const auto &[self_idx, other_idx] : cross_edges)
    {
      const auto &self_item = items[self_idx];
      const auto &other_item = other[other_idx];

      float min_dist = (self_item.radius + other_item.radius) * overlap_factor;
      float dx = self_item.position.x - other_item.position.x;
      float dy = self_item.position.y - other_item.position.y;
      float dist_sq = dx * dx + dy * dy;

      if (dist_sq < min_dist * min_dist)
      {
        had_collision = true;
        float     dist = std::sqrt(dist_sq);
        glm::vec2 dir;

        if (dist > 1e-6f)
        {
          dir = glm::vec2(dx / dist, dy / dist);
        }
        else
        {
          float angle = static_cast<float>((self_idx + other_idx + iter) %
                                           360) *
                        (3.14159265f / 180.0f);
          dir = glm::vec2(std::cos(angle), std::sin(angle));
          dist = 0.0f;
        }

        float     penetration = min_dist - dist;
        glm::vec2 delta = step_size * penetration * dir;

        items[self_idx].position.x += delta.x;
        items[self_idx].position.y += delta.y;
      }
    }

    if (!had_collision) break;
  }
}

void ScatterField::set_elevation_from_terrain(const Array     &elevation,
                                              const glm::vec4 &bbox)
{
  if (!validate_non_empty(elevation)) return;

  float dx = bbox.y - bbox.x;
  float dy = bbox.w - bbox.z;
  if (std::abs(dx) < 1e-7f || std::abs(dy) < 1e-7f) return;

  for (auto &item : items)
  {
    // scale to unit interval
    float xn = (item.position.x - bbox.x) / dx;
    float yn = (item.position.y - bbox.z) / dy;

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
      item.position.z = elevation.get_value_bilinear_at(i, j, u, v);
    }
  }
}

void ScatterField::shuffle_classes(float    ratio,
                                   size_t   k_neighbors,
                                   uint32_t seed)
{
  if (items.size() < 2 || ratio <= 0.0f || k_neighbors == 0) return;

  // --- Build 2D Point List and Query Nearest Neighbors

  std::vector<ps::Point<float, 2>> points;
  points.reserve(items.size());
  for (const auto &item : items)
  {
    points.push_back({item.position.x, item.position.y});
  }

  size_t safe_k = std::min(k_neighbors, items.size() - 1);
  auto   neighbors_idx = ps::nearest_neighbors_indices(points, safe_k);

  // --- Select Candidate Items to Shuffle

  std::mt19937 gen(seed);

  std::vector<size_t> perm(items.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::shuffle(perm.begin(), perm.end(), gen);

  size_t target_count = std::min(
      items.size(),
      static_cast<size_t>(
          std::round(ratio * static_cast<float>(items.size()))));

  // --- Perform Neighbor Class Swaps with Failsafe

  for (size_t c = 0; c < target_count; ++c)
  {
    size_t i = perm[c];

    // find neighbors with a differing class
    std::vector<size_t> valid_neighbors;
    for (size_t neighbor_idx : neighbors_idx[i])
    {
      if (items[neighbor_idx].class_id != items[i].class_id)
      {
        valid_neighbors.push_back(neighbor_idx);
      }
    }

    // if all neighbors share the same class, skip candidate
    if (valid_neighbors.empty()) continue;

    // pick one differing neighbor at random and swap class & radius
    std::uniform_int_distribution<size_t> dis(0, valid_neighbors.size() - 1);
    size_t chosen_neighbor = valid_neighbors[dis(gen)];

    std::swap(items[i].class_id, items[chosen_neighbor].class_id);
    std::swap(items[i].radius, items[chosen_neighbor].radius);
  }
}

Cloud ScatterField::to_cloud() const
{
  std::vector<Point> pts;
  pts.reserve(items.size());
  for (const auto &item : items)
    pts.push_back(item.to_point());
  return Cloud(std::move(pts));
}

void ScatterField::to_csv(const std::string &fname) const
{
  std::ofstream f(fname, std::ios::out);
  if (!f.is_open()) throw std::runtime_error("Failed to open file: " + fname);

  f.imbue(std::locale("C"));
  f << std::fixed << std::setprecision(9);

  for (const auto &item : items)
  {
    f << item.position.x << ',' << item.position.y << ',' << item.position.z
      << ',' << item.class_id << ',' << item.radius << '\n';
  }
}

Array ScatterField::to_density_map(glm::ivec2              shape,
                                   float                   sigma,
                                   std::optional<uint32_t> class_id,
                                   bool                    weighted_by_area,
                                   glm::vec4               bbox) const
{
  if (!validate_shape(shape)) return Array();

  float dx = bbox.y - bbox.x;
  float dy = bbox.w - bbox.z;
  if (std::abs(dx) < 1e-9f || std::abs(dy) < 1e-9f) return Array(shape, 0.0f);

  // --- Determine effective Gaussian sigma (in world coordinates)

  float pixel_size_x = dx / static_cast<float>(shape.x);
  float pixel_size_y = dy / static_cast<float>(shape.y);
  float mean_pixel_size = 0.5f * (pixel_size_x + pixel_size_y);

  float effective_sigma = sigma;
  if (effective_sigma <= 0.0f)
  {
    effective_sigma = std::max(2.5f * mean_pixel_size,
                               0.03f * std::min(dx, dy));
  }

  float sigma_sq = effective_sigma * effective_sigma;
  float inv_2sigma_sq = 1.0f / (2.0f * sigma_sq);
  float norm_factor = 1.0f / (2.0f * static_cast<float>(M_PI) * sigma_sq);

  // Splat footprint radius (3 sigma in pixels)
  int radius_px_x = static_cast<int>(
      std::ceil(3.0f * effective_sigma / pixel_size_x));
  int radius_px_y = static_cast<int>(
      std::ceil(3.0f * effective_sigma / pixel_size_y));

  Array density(shape, 0.0f);

  for (const auto &item : items)
  {
    if (class_id.has_value() && item.class_id != class_id.value())
    {
      continue;
    }

    float weight = 1.0f;
    if (weighted_by_area)
    {
      weight = static_cast<float>(M_PI) * item.radius * item.radius;
    }

    // Item center in continuous pixel coordinates
    float u = (item.position.x - bbox.x) / dx;
    float v = (item.position.y - bbox.z) / dy;

    float center_px_x = u * static_cast<float>(shape.x - 1);
    float center_px_y = v * static_cast<float>(shape.y - 1);

    int min_i = std::max(
        0,
        static_cast<int>(std::floor(center_px_x - float(radius_px_x))));
    int max_i = std::min(
        shape.x - 1,
        static_cast<int>(std::ceil(center_px_x + float(radius_px_x))));
    int min_j = std::max(
        0,
        static_cast<int>(std::floor(center_px_y - float(radius_px_y))));
    int max_j = std::min(
        shape.y - 1,
        static_cast<int>(std::ceil(center_px_y + float(radius_px_y))));

    for (int j = min_j; j <= max_j; ++j)
    {
      float y_world = bbox.z + (static_cast<float>(j) /
                                static_cast<float>(shape.y - 1)) *
                                   dy;
      float diff_y = y_world - item.position.y;
      float diff_y_sq = diff_y * diff_y;

      for (int i = min_i; i <= max_i; ++i)
      {
        float x_world = bbox.x + (static_cast<float>(i) /
                                  static_cast<float>(shape.x - 1)) *
                                     dx;
        float diff_x = x_world - item.position.x;
        float dist_sq = diff_x * diff_x + diff_y_sq;

        float g = norm_factor * std::exp(-dist_sq * inv_2sigma_sq);
        density(i, j) += weight * g;
      }
    }
  }

  return density;
}

Array ScatterField::to_heightmap(glm::ivec2              shape,
                                 ScatterShape            shape_type,
                                 float                   height_radius_ratio,
                                 std::optional<uint32_t> class_id,
                                 Array                  *p_rock_map,
                                 uint32_t                seed,
                                 glm::vec4               bbox) const
{
  if (!validate_shape(shape))
  {
    if (p_rock_map) *p_rock_map = Array();
    return Array();
  }
  if (p_rock_map) *p_rock_map = Array(shape, 0.0f);
  if (items.empty()) return Array(shape, 0.0f);

  int   num_items = static_cast<int>(items.size());
  Array heightmap(shape, 0.0f);

  std::vector<float>    pos_rad(num_items * 4);
  std::vector<uint32_t> class_ids(num_items);

  for (int i = 0; i < num_items; ++i)
  {
    pos_rad[i * 4 + 0] = items[i].position.x;
    pos_rad[i * 4 + 1] = items[i].position.y;
    pos_rad[i * 4 + 2] = items[i].position.z;
    pos_rad[i * 4 + 3] = items[i].radius;
    class_ids[i] = items[i].class_id;
  }

  int filter_class = class_id.has_value() ? static_cast<int>(*class_id) : -1;

  if (clwrapper::DeviceManager::get_instance().is_ready())
  {
    auto run = clwrapper::Run("scatter_to_heightmap");
    run.bind_buffer<float>("heightmap", heightmap.vector);
    gpu::helper_bind_optional_buffer(run, "rock_map", p_rock_map);
    run.bind_buffer<float>("items_pos_rad", pos_rad);
    run.bind_buffer<uint32_t>("items_class", class_ids);

    run.bind_arguments(num_items,
                       shape.x,
                       shape.y,
                       bbox,
                       static_cast<int>(shape_type),
                       height_radius_ratio,
                       seed,
                       filter_class,
                       p_rock_map ? 1 : 0);

    run.write_buffer("heightmap");
    if (p_rock_map)
    {
      run.write_buffer("rock_map");
    }
    run.write_buffer("items_pos_rad");
    run.write_buffer("items_class");

    run.execute_async(num_items);
    run.finish();

    run.read_buffer("heightmap");
    if (p_rock_map)
    {
      run.read_buffer("rock_map");
      clamp(*p_rock_map, 0.0f, 1.0f);
    }
  }
  else
  {
    // fallback CPU implementation
    float xmin = bbox.x;
    float xmax = bbox.y;
    float ymin = bbox.z;
    float ymax = bbox.w;
    float dx_dom = xmax - xmin;
    float dy_dom = ymax - ymin;
    if (dx_dom <= 0.0f || dy_dom <= 0.0f) return heightmap;

    float px = dx_dom / static_cast<float>(shape.x);
    float py = dy_dom / static_cast<float>(shape.y);

    for (int i = 0; i < num_items; ++i)
    {
      if (filter_class >= 0 &&
          class_ids[i] != static_cast<uint32_t>(filter_class))
      {
        continue;
      }

      float cx = pos_rad[i * 4 + 0];
      float cy = pos_rad[i * 4 + 1];
      float radius = pos_rad[i * 4 + 3];
      if (radius <= 0.0f) continue;

      int ix_min = std::clamp(
          static_cast<int>(std::floor((cx - radius - xmin) / px)),
          0,
          shape.x - 1);
      int ix_max = std::clamp(
          static_cast<int>(std::ceil((cx + radius - xmin) / px)),
          0,
          shape.x - 1);
      int iy_min = std::clamp(
          static_cast<int>(std::floor((cy - radius - ymin) / py)),
          0,
          shape.y - 1);
      int iy_max = std::clamp(
          static_cast<int>(std::ceil((cy + radius - ymin) / py)),
          0,
          shape.y - 1);

      float h_max = radius * height_radius_ratio;

      constexpr int max_verts = 8;
      int           n_verts = 6;
      float         thetas[8];
      float         r_verts[8];

      if (shape_type == SCATTER_SHAPE_POLYGON ||
          shape_type == SCATTER_SHAPE_PYRAMID)
      {
        std::mt19937 rng(seed ^ (static_cast<uint32_t>(i) * 1999u + 17u));
        std::uniform_real_distribution<float> dist_uniform(0.0f, 1.0f);
        n_verts = 5 + static_cast<int>(dist_uniform(rng) * 3.99f);
        if (n_verts > max_verts) n_verts = max_verts;
        float base_rot = dist_uniform(rng) * 2.0f * static_cast<float>(M_PI);

        for (int k = 0; k < n_verts; ++k)
        {
          thetas[k] = base_rot + (2.0f * static_cast<float>(M_PI) *
                                  static_cast<float>(k)) /
                                     static_cast<float>(n_verts);
          r_verts[k] = radius * (0.65f + 0.35f * dist_uniform(rng));
        }
      }

      for (int iy = iy_min; iy <= iy_max; ++iy)
      {
        float wy = ymin + (static_cast<float>(iy) + 0.5f) * py;
        float dy = wy - cy;

        for (int ix = ix_min; ix <= ix_max; ++ix)
        {
          float wx = xmin + (static_cast<float>(ix) + 0.5f) * px;
          float dx = wx - cx;
          float dist = std::sqrt(dx * dx + dy * dy);
          float shape_val = 0.0f;

          if (shape_type == SCATTER_SHAPE_DISK)
          {
            if (dist < radius)
            {
              float u = dist / radius;
              shape_val = std::sqrt(std::max(0.0f, 1.0f - u * u));
            }
          }
          else if (shape_type == SCATTER_SHAPE_CONE)
          {
            if (dist < radius)
            {
              float u = dist / radius;
              shape_val = (1.0f - u);
            }
          }
          else if (shape_type == SCATTER_SHAPE_SMOOTH_DOME)
          {
            if (dist < radius)
            {
              float u = dist / radius;
              shape_val = 0.5f *
                          (1.0f + std::cos(u * static_cast<float>(M_PI)));
            }
          }
          else if (shape_type == SCATTER_SHAPE_POLYGON ||
                   shape_type == SCATTER_SHAPE_PYRAMID)
          {
            if (dist < 1e-7f)
            {
              shape_val = 1.0f;
            }
            else
            {
              float phi = std::atan2(dy, dx);
              if (phi < 0.0f) phi += 2.0f * static_cast<float>(M_PI);

              float r_poly = radius;
              for (int k = 0; k < n_verts; ++k)
              {
                int   next_k = (k + 1) % n_verts;
                float t1 = thetas[k];
                float t2 = thetas[next_k];
                while (t1 < 0.0f)
                  t1 += 2.0f * static_cast<float>(M_PI);
                while (t1 >= 2.0f * static_cast<float>(M_PI))
                  t1 -= 2.0f * static_cast<float>(M_PI);
                while (t2 < 0.0f)
                  t2 += 2.0f * static_cast<float>(M_PI);
                while (t2 >= 2.0f * static_cast<float>(M_PI))
                  t2 -= 2.0f * static_cast<float>(M_PI);

                bool in_sector = (t1 < t2) ? (phi >= t1 && phi <= t2)
                                           : (phi >= t1 || phi <= t2);
                if (in_sector)
                {
                  float d_theta = t2 - t1;
                  if (d_theta < 0.0f)
                    d_theta += 2.0f * static_cast<float>(M_PI);
                  float d_phi1 = phi - t1;
                  if (d_phi1 < 0.0f) d_phi1 += 2.0f * static_cast<float>(M_PI);
                  float d_phi2 = t2 - phi;
                  if (d_phi2 < 0.0f) d_phi2 += 2.0f * static_cast<float>(M_PI);

                  float r1 = r_verts[k];
                  float r2 = r_verts[next_k];
                  float denom = r2 * std::sin(d_phi2) + r1 * std::sin(d_phi1);
                  if (denom > 1e-6f)
                  {
                    r_poly = (r1 * r2 * std::sin(d_theta)) / denom;
                  }
                  break;
                }
              }

              if (dist < r_poly)
              {
                float u = dist / r_poly;
                if (shape_type == SCATTER_SHAPE_POLYGON)
                {
                  shape_val = std::sqrt(std::max(0.0f, 1.0f - u * u));
                }
                else
                {
                  shape_val = (1.0f - u);
                }
              }
            }
          }

          if (shape_val > 0.0f)
          {
            float dz = h_max * shape_val;
            if (dz > 0.0f)
            {
              heightmap(ix, iy) += dz;
            }
            if (p_rock_map)
            {
              (*p_rock_map)(ix, iy) += shape_val;
            }
          }
        }
      }
    }

    if (p_rock_map)
    {
      clamp(*p_rock_map, 0.0f, 1.0f);
    }
  }

  return heightmap;
}

// --- Internal Rendering Helper

static cv::Mat render_scatter_field_mat(const std::vector<ScatterItem> &items,
                                        glm::ivec2                      shape,
                                        const Array &background,
                                        glm::vec4    bbox)
{
  if (!validate_shape(shape)) return {};

  cv::Mat img;

  if (validate_non_empty(background))
  {
    float vmin = background.min();
    float vmax = background.max();
    float range = (std::abs(vmax - vmin) > 1e-6f) ? (vmax - vmin) : 1.0f;

    cv::Mat bg_float(background.shape.y,
                     background.shape.x,
                     CV_32FC1,
                     const_cast<float *>(background.vector.data()));

    cv::Mat bg_resized;
    if (background.shape.x != shape.x || background.shape.y != shape.y)
    {
      cv::resize(bg_float,
                 bg_resized,
                 cv::Size(shape.x, shape.y),
                 0,
                 0,
                 cv::INTER_LINEAR);
    }
    else
    {
      bg_resized = bg_float.clone();
    }

    cv::Mat bg_8u;
    bg_resized.convertTo(bg_8u, CV_8UC1, 255.0 / range, -vmin * 255.0 / range);

    // flip vertically so row 0 is top
    cv::flip(bg_8u, bg_8u, 0);

    cv::cvtColor(bg_8u, img, cv::COLOR_GRAY2BGR);
  }
  else
  {
    img = cv::Mat(shape.y, shape.x, CV_8UC3, cv::Scalar(35, 35, 35));
  }

  // --- Draw Scatter Items

  float width = bbox.y - bbox.x;
  float height = bbox.w - bbox.z;
  if (width <= 1e-7f) width = 1.0f;
  if (height <= 1e-7f) height = 1.0f;

  float scale_x = static_cast<float>(shape.x) / width;
  float scale_y = static_cast<float>(shape.y) / height;
  float scale = 0.5f * (scale_x + scale_y);

  for (const auto &item : items)
  {
    float u = (item.position.x - bbox.x) / width;
    float v = (item.position.y - bbox.z) / height;

    int px = static_cast<int>(std::round(u * static_cast<float>(shape.x - 1)));
    int py = static_cast<int>(
        std::round((1.0f - v) * static_cast<float>(shape.y - 1)));

    int r_px = std::max(1, static_cast<int>(std::round(item.radius * scale)));

    cv::Scalar fill_color = get_scatter_class_color(item.class_id);
    cv::Scalar edge_color = cv::Scalar(fill_color[0] * 0.5,
                                       fill_color[1] * 0.5,
                                       fill_color[2] * 0.5);

    cv::circle(img,
               cv::Point(px, py),
               r_px,
               fill_color,
               cv::FILLED,
               cv::LINE_AA);
    if (r_px > 1)
    {
      cv::circle(img, cv::Point(px, py), r_px, edge_color, 1, cv::LINE_AA);
    }
  }

  return img;
}

std::vector<uint8_t> ScatterField::to_img_8bit(glm::ivec2   shape,
                                               const Array &background,
                                               glm::vec4    bbox,
                                               bool         flip_y) const
{
  cv::Mat img = render_scatter_field_mat(items, shape, background, bbox);
  if (img.empty()) return {};

  if (!flip_y)
  {
    cv::flip(img, img, 0);
  }

  cv::Mat rgb;
  cv::cvtColor(img, rgb, cv::COLOR_BGR2RGB);

  if (!rgb.isContinuous()) rgb = rgb.clone();

  std::vector<uint8_t> buffer(rgb.total() * rgb.elemSize());
  std::memcpy(buffer.data(), rgb.data, buffer.size());
  return buffer;
}

void ScatterField::to_png(const std::string &fname,
                          glm::ivec2         shape,
                          const Array       &background,
                          glm::vec4          bbox) const
{
  cv::Mat img = render_scatter_field_mat(items, shape, background, bbox);
  if (img.empty()) return;
  cv::imwrite(fname, img);
}

std::string ScatterField::to_string() const
{
  std::ostringstream ss;
  ss.imbue(std::locale("C"));

  glm::vec4 bbox = get_bbox();

  std::map<uint32_t, size_t> class_counts;
  std::map<uint32_t, float>  class_radii_sum;
  for (const auto &item : items)
  {
    class_counts[item.class_id]++;
    class_radii_sum[item.class_id] += item.radius;
  }

  ss << "ScatterField Infos\n";
  ss << "--------------------------------\n";
  ss << " total items : " << items.size() << "\n";
  ss << " class count : " << class_counts.size() << "\n";
  ss << " bbox : {" << bbox.x << ", " << bbox.y << ", " << bbox.z << ", "
     << bbox.w << "}\n";

  if (!class_counts.empty())
  {
    ss << " class breakdown :\n";
    for (const auto &[c_id, count] : class_counts)
    {
      float avg_r = (count > 0) ? (class_radii_sum[c_id] / float(count)) : 0.f;
      float pct = (items.empty())
                      ? 0.f
                      : (100.f * float(count) / float(items.size()));
      ss << "   - class " << c_id << ": " << count << " items (" << std::fixed
         << std::setprecision(1) << pct << "%, avg radius "
         << std::setprecision(4) << avg_r << ")\n";
    }
  }
  ss << "--------------------------------";

  return ss.str();
}

// ============================================================================
//  Functions (Alphabetically Sorted)
// ============================================================================

ScatterField merge_scatter_field(const ScatterField &field1,
                                 const ScatterField &field2,
                                 bool                merge_by_class)
{
  return merge_scatter_fields({field1, field2}, merge_by_class);
}

ScatterField merge_scatter_fields(const std::vector<ScatterField> &fields,
                                  bool merge_by_class)
{
  ScatterField result;
  size_t       total_size = 0;
  for (const auto &f : fields)
    total_size += f.size();

  result.reserve(total_size);

  for (size_t field_idx = 0; field_idx < fields.size(); ++field_idx)
  {
    const auto &f = fields[field_idx];
    for (const auto &item : f)
    {
      if (merge_by_class)
      {
        result.push_back(item);
      }
      else
      {
        ScatterItem remapped_item = item;
        remapped_item.class_id = static_cast<uint32_t>(field_idx);
        result.push_back(std::move(remapped_item));
      }
    }
  }

  return result;
}

Array scatter_field_to_heightmap(const ScatterField     &field,
                                 glm::ivec2              shape,
                                 ScatterShape            shape_type,
                                 float                   height_radius_ratio,
                                 std::optional<uint32_t> class_id,
                                 Array                  *p_rock_map,
                                 uint32_t                seed,
                                 glm::vec4               bbox)
{
  return field.to_heightmap(shape,
                            shape_type,
                            height_radius_ratio,
                            class_id,
                            p_rock_map,
                            seed,
                            bbox);
}

} // namespace hmap
