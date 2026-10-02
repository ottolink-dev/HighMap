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

#include "highmap/flora/forest.hpp"
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

void Forest::rejection_filter_density(const Array     &density_mask,
                                      uint32_t         seed,
                                      const glm::vec4 &bbox)
{
  if (!validate_non_empty(density_mask)) return;

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(0.f, 1.f);

  auto density_fct = make_xy_function_from_array(density_mask, bbox);

  trees.erase(std::remove_if(trees.begin(),
                             trees.end(),
                             [&](const Tree &tree)
                             {
                               float rnd = dis(gen);
                               return (rnd > density_fct(tree.position.x,
                                                         tree.position.y));
                             }),
              trees.end());
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
      float avg_r = (count > 0) ? (species_radii_sum[sp_id] / float(count)) : 0.f;
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
