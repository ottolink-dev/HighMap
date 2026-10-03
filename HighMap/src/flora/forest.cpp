/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include "delaunator-cpp.hpp"
#include "point_sampler/metrics.hpp"

#include "highmap/flora/forest.hpp"
#include "highmap/flora/species.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

// --- Helper Functions

static cv::Scalar get_species_color(uint32_t species_id)
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

  if (species_id < palette.size()) return palette[species_id];

  // golden ratio hue rotation for arbitrary species counts
  float   hue = std::fmod(static_cast<float>(species_id) * 137.508f, 180.0f);
  cv::Mat hsv(1, 1, CV_8UC3, cv::Scalar(static_cast<uint8_t>(hue), 220, 240));
  cv::Mat bgr;
  cv::cvtColor(hsv, bgr, cv::COLOR_HSV2BGR);
  cv::Vec3b vec = bgr.at<cv::Vec3b>(0, 0);
  return cv::Scalar(vec[0], vec[1], vec[2]);
}

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

void Forest::densify(float default_radius)
{
  if (trees.size() < 3) return;

  // --- Prepare 2D Coordinates for Delaunay Triangulation

  std::vector<double> coords;
  coords.reserve(2 * trees.size());
  for (const auto &tree : trees)
  {
    coords.push_back(static_cast<double>(tree.position.x));
    coords.push_back(static_cast<double>(tree.position.y));
  }

  delaunator::Delaunator d(coords);
  const auto            &tri = d.triangles;
  if (tri.empty()) return;

  // --- Compute Circumcenters & Assign Majority Species

  std::vector<Tree> new_trees;
  new_trees.reserve(tri.size() / 3);

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

    // Determine majority species id among the 3 triangle vertices
    uint32_t s0 = trees[i0].species_id;
    uint32_t s1 = trees[i1].species_id;
    uint32_t s2 = trees[i2].species_id;

    uint32_t majority_species = s0;
    if (s0 == s1 || s0 == s2)
    {
      majority_species = s0;
    }
    else if (s1 == s2)
    {
      majority_species = s1;
    }
    else
    {
      majority_species = s0;
    }

    new_trees.emplace_back(static_cast<float>(cx_center),
                           static_cast<float>(cy_center),
                           0.0f,
                           majority_species,
                           default_radius);
  }

  // --- Append New Trees to Forest

  trees.insert(trees.end(),
               std::make_move_iterator(new_trees.begin()),
               std::make_move_iterator(new_trees.end()));
}

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

void Forest::prune_density(const Array     &density_mask,
                           float            target_ratio,
                           uint32_t         seed,
                           const glm::vec4 &bbox)
{
  if (!validate_non_empty(density_mask) || trees.empty()) return;

  float r_target = std::clamp(target_ratio, 0.0f, 1.0f);
  if (r_target <= 1e-6f)
  {
    trees.clear();
    return;
  }
  if (r_target >= 1.0f - 1e-6f)
  {
    return;
  }

  // --- Sample Local Density at Tree Positions

  auto density_fct = make_xy_function_from_array(density_mask, bbox);
  std::vector<float> sampled_densities(trees.size());
  for (size_t i = 0; i < trees.size(); ++i)
  {
    sampled_densities[i] = std::max(
        0.0f,
        density_fct(trees[i].position.x, trees[i].position.y));
  }

  // --- Remap Acceptance Probabilities to Match Target Keep Ratio

  // find shift c in [-1, 1] such that expected ratio
  // mean(clamp(d + c, 0, 1)) approximately equals r_target
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
    float expected_ratio = sum_p / static_cast<float>(trees.size());

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

  std::vector<Tree> retained;
  retained.reserve(
      static_cast<size_t>(std::ceil(r_target * float(trees.size()))));

  for (size_t i = 0; i < trees.size(); ++i)
  {
    float p_accept = std::clamp(sampled_densities[i] + c_best, 0.0f, 1.0f);
    if (dis(gen) <= p_accept)
    {
      retained.push_back(trees[i]);
    }
  }

  trees = std::move(retained);
}

void Forest::prune_unviable(const std::vector<Species> &species,
                            bool                        prune_collisions)
{
  if (trees.size() < 2) return;

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

  // --- KD-Tree Nearest Neighbor Query

  std::vector<ps::Point<float, 2>> points;
  points.reserve(trees.size());
  for (const auto &tree : trees)
  {
    points.push_back({tree.position.x, tree.position.y});
  }

  auto neighbors_idx = ps::nearest_neighbors_indices(points, 1);

  // --- Check Minimum Radius / Viability

  std::vector<Tree> viable;
  viable.reserve(trees.size());

  for (size_t i = 0; i < trees.size(); ++i)
  {
    const Tree &tree = trees[i];
    Species     sp_i = get_species_traits(tree);

    size_t neighbor_idx = neighbors_idx[i].empty() ? i : neighbors_idx[i][0];
    const Tree &neighbor_tree = trees[neighbor_idx];

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

  // --- Collision Pruning (Optional)

  if (!prune_collisions || viable.size() < 2)
  {
    trees = std::move(viable);
    return;
  }

  // Sort candidates by radius descending
  std::vector<size_t> order(viable.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(),
                   order.end(),
                   [&](size_t a, size_t b)
                   { return viable[a].radius > viable[b].radius; });

  float max_r = 0.0f;
  for (const auto &t : viable)
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

  std::vector<Tree> retained;
  retained.reserve(viable.size());

  for (size_t idx : order)
  {
    const Tree &cand = viable[idx];
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
          const Tree &other = retained[ret_idx];
          float       dist_x = cand.position.x - other.position.x;
          float       dist_y = cand.position.y - other.position.y;
          float       dist_sq = dist_x * dist_x + dist_y * dist_y;

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

  trees = std::move(retained);
}

void Forest::reinforce_species_clusters(size_t iterations,
                                        size_t k_neighbors,
                                        bool   include_self)
{
  if (trees.size() < 2 || k_neighbors == 0 || iterations == 0) return;

  // --- Query Nearest Neighbors Graph Once (Static Positions)

  std::vector<ps::Point<float, 2>> points;
  points.reserve(trees.size());
  for (const auto &tree : trees)
  {
    points.push_back({tree.position.x, tree.position.y});
  }

  size_t safe_k = std::min(k_neighbors, trees.size() - 1);
  auto   neighbors_idx = ps::nearest_neighbors_indices(points, safe_k);

  // --- Iterative Majority Species Assignment

  std::vector<uint32_t> next_species(trees.size());

  for (size_t iter = 0; iter < iterations; ++iter)
  {
    for (size_t i = 0; i < trees.size(); ++i)
    {
      std::unordered_map<uint32_t, size_t> counts;

      if (include_self)
      {
        counts[trees[i].species_id]++;
      }

      for (size_t neighbor_idx : neighbors_idx[i])
      {
        counts[trees[neighbor_idx].species_id]++;
      }

      // Pick dominant species (majority vote)
      uint32_t best_species = trees[i].species_id;
      size_t   max_count = 0;

      // Check current tree species first to favor status quo in case of ties
      auto self_it = counts.find(trees[i].species_id);
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

    for (size_t i = 0; i < trees.size(); ++i)
    {
      trees[i].species_id = next_species[i];
    }
  }
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

void Forest::shuffle_species(float ratio, size_t k_neighbors, uint32_t seed)
{
  if (trees.size() < 2 || ratio <= 0.0f || k_neighbors == 0) return;

  // --- Build 2D Point List and Query Nearest Neighbors

  std::vector<ps::Point<float, 2>> points;
  points.reserve(trees.size());
  for (const auto &tree : trees)
  {
    points.push_back({tree.position.x, tree.position.y});
  }

  size_t safe_k = std::min(k_neighbors, trees.size() - 1);
  auto   neighbors_idx = ps::nearest_neighbors_indices(points, safe_k);

  // --- Select Candidate Trees to Shuffle

  std::mt19937 gen(seed);

  std::vector<size_t> perm(trees.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::shuffle(perm.begin(), perm.end(), gen);

  size_t target_count = std::min(
      trees.size(),
      static_cast<size_t>(
          std::round(ratio * static_cast<float>(trees.size()))));

  // --- Perform Neighbor Species Swaps with Failsafe

  for (size_t c = 0; c < target_count; ++c)
  {
    size_t i = perm[c];

    // find neighbors with a differing species
    std::vector<size_t> valid_neighbors;
    for (size_t neighbor_idx : neighbors_idx[i])
    {
      if (trees[neighbor_idx].species_id != trees[i].species_id)
      {
        valid_neighbors.push_back(neighbor_idx);
      }
    }

    // if all neighbors share the same species, skip candidate
    if (valid_neighbors.empty()) continue;

    // pick one differing neighbor at random and swap species & radius
    std::uniform_int_distribution<size_t> dis(0, valid_neighbors.size() - 1);
    size_t chosen_neighbor = valid_neighbors[dis(gen)];

    std::swap(trees[i].species_id, trees[chosen_neighbor].species_id);
    std::swap(trees[i].radius, trees[chosen_neighbor].radius);
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

void Forest::to_csv(const std::string &fname) const
{
  std::ofstream f(fname, std::ios::out);
  if (!f.is_open()) throw std::runtime_error("Failed to open file: " + fname);

  // use C locale for consistent number formatting
  f.imbue(std::locale("C"));
  f << std::fixed << std::setprecision(9);

  for (const auto &tree : trees)
  {
    f << tree.position.x << ',' << tree.position.y << ',' << tree.position.z
      << ',' << tree.species_id << ',' << tree.radius << '\n';
  }
}

Array Forest::to_density_map(glm::ivec2              shape,
                             float                   sigma,
                             std::optional<uint32_t> species_id,
                             bool                    weighted_by_crown,
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
    // Auto bandwidth: Silverman's rule of thumb heuristic or 2.5 grid cells
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

  for (const auto &tree : trees)
  {
    if (species_id.has_value() && tree.species_id != species_id.value())
    {
      continue;
    }

    float weight = 1.0f;
    if (weighted_by_crown)
    {
      weight = static_cast<float>(M_PI) * tree.radius * tree.radius;
    }

    // Tree center in continuous pixel coordinates
    float u = (tree.position.x - bbox.x) / dx;
    float v = (tree.position.y - bbox.z) / dy;

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
      float diff_y = y_world - tree.position.y;
      float diff_y_sq = diff_y * diff_y;

      for (int i = min_i; i <= max_i; ++i)
      {
        float x_world = bbox.x + (static_cast<float>(i) /
                                  static_cast<float>(shape.x - 1)) *
                                     dx;
        float diff_x = x_world - tree.position.x;
        float dist_sq = diff_x * diff_x + diff_y_sq;

        float g = norm_factor * std::exp(-dist_sq * inv_2sigma_sq);
        density(i, j) += weight * g;
      }
    }
  }

  return density;
}

void Forest::to_png(const std::string &fname,
                    glm::ivec2         shape,
                    const Array       &background,
                    glm::vec4          bbox) const
{
  if (!validate_shape(shape)) return;

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

  // --- Draw Trees

  float width = bbox.y - bbox.x;
  float height = bbox.w - bbox.z;
  if (width <= 1e-7f) width = 1.0f;
  if (height <= 1e-7f) height = 1.0f;

  float scale_x = static_cast<float>(shape.x) / width;
  float scale_y = static_cast<float>(shape.y) / height;
  float scale = 0.5f * (scale_x + scale_y);

  for (const auto &tree : trees)
  {
    float u = (tree.position.x - bbox.x) / width;
    float v = (tree.position.y - bbox.z) / height;

    int px = static_cast<int>(std::round(u * static_cast<float>(shape.x - 1)));
    int py = static_cast<int>(
        std::round((1.0f - v) * static_cast<float>(shape.y - 1)));

    int r_px = std::max(1, static_cast<int>(std::round(tree.radius * scale)));

    cv::Scalar fill_color = get_species_color(tree.species_id);
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

  cv::imwrite(fname, img);
}

std::string Forest::to_string() const
{
  std::ostringstream ss;
  ss.imbue(std::locale("C"));

  glm::vec4 bbox = get_bbox();

  // count trees and collect average radius per species
  std::map<uint32_t, size_t> species_counts;
  std::map<uint32_t, float>  species_radii_sum;
  for (const auto &tree : trees)
  {
    species_counts[tree.species_id]++;
    species_radii_sum[tree.species_id] += tree.radius;
  }

  ss << "Forest Infos\n";
  ss << "--------------------------------\n";
  ss << " total trees : " << trees.size() << "\n";
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
      float pct = (trees.empty())
                      ? 0.f
                      : (100.f * float(count) / float(trees.size()));
      ss << "   - species " << sp_id << ": " << count << " trees ("
         << std::fixed << std::setprecision(1) << pct << "%, avg radius "
         << std::setprecision(4) << avg_r << ")\n";
    }
  }
  ss << "--------------------------------";

  return ss.str();
}

} // namespace hmap
