/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include "point_sampler/halton.hpp"
#include "point_sampler/hammersley.hpp"
#include "point_sampler/jittered_grid.hpp"
#include "point_sampler/latin_hypercube_sampling.hpp"
#include "point_sampler/point.hpp"
#include "point_sampler/poisson_disk_sampling.hpp"
#include "point_sampler/random.hpp"
#include "point_sampler/rejection_sampling.hpp"
#include "point_sampler/utils.hpp"

#include "highmap/array.hpp"
#include "highmap/geometry/point_sampling.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

size_t helper_estimate_count(const glm::vec4 &bbox, float distance)
{
  size_t count = static_cast<size_t>(2.f * (bbox.y - bbox.x) / distance *
                                     (bbox.w - bbox.z) / distance);
  return count;
}

std::array<std::pair<float, float>, 2> bbox_to_ranges2d(const glm::vec4 &bbox)
{
  std::array<std::pair<float, float>, 2> ranges = {
      std::make_pair(bbox.x, bbox.y),
      std::make_pair(bbox.z, bbox.w)};
  return ranges;
}

std::function<float(const ps::Point<float, 2> &)>
make_pointwise_function_from_array(const Array &array, const glm::vec4 &bbox)
{
  if (!validate_non_empty(array))
    return [](const ps::Point<float, 2> &) -> float { return 0.f; };

  return [&array, &bbox](const ps::Point<float, 2> &p) -> float
  {
    float x = (p[0] - bbox.x) / (bbox.y - bbox.x);
    float y = (p[1] - bbox.z) / (bbox.w - bbox.z);

    x = std::clamp(x, 0.f, 1.f);
    y = std::clamp(y, 0.f, 1.f);

    float xn = x * (array.shape.x - 1);
    float yn = y * (array.shape.y - 1);

    int   i = static_cast<int>(xn);
    int   j = static_cast<int>(yn);
    float u = xn - i;
    float v = yn - j;

    return array.get_value_bilinear_at(i, j, u, v);
  };
}

std::array<std::vector<float>, 2> random_points(
    size_t                     count,
    std::uint32_t              seed,
    const PointSamplingMethod &method,
    const glm::vec4           &bbox)
{
  std::vector<ps::Point<float, 2>> points;
  auto                             ranges = bbox_to_ranges2d(bbox);

  switch (method)
  {
  case PointSamplingMethod::RND_RANDOM:
  {
    points = ps::random<float, 2>(count, ranges, seed);
  }
  break;
  //
  case PointSamplingMethod::RND_HALTON:
  {
    points = ps::halton<float, 2>(count, ranges, seed);
  }
  break;
  //
  case PointSamplingMethod::RND_HAMMERSLEY:
  {
    points = ps::hammersley<float, 2>(count, ranges, seed);
  }
  break;
  //
  case PointSamplingMethod::RND_LHS:
  {
    points = ps::latin_hypercube_sampling<float, 2>(count, ranges, seed);
  }
  break;
    //
  }

  return ps::split_by_dimension(points);
}

std::array<std::vector<float>, 2> random_points_density(size_t        count,
                                                        const Array  &density,
                                                        std::uint32_t seed,
                                                        const glm::vec4 &bbox)
{
  if (!validate_non_empty(density))
    return {std::vector<float>{}, std::vector<float>{}};

  auto ranges = bbox_to_ranges2d(bbox);
  auto density_fct = make_pointwise_function_from_array(density, bbox);

  auto points = ps::rejection_sampling<float, 2>(count,
                                                 ranges,
                                                 density_fct,
                                                 seed);
  return ps::split_by_dimension(points);
}

std::array<std::vector<float>, 2> random_points_inverse_sampling(
    size_t           count,
    const Array     &array,
    std::uint32_t    seed,
    const glm::vec4 &bbox)
{
  if (!validate_non_empty(array) || count == 0)
    return {std::vector<float>{}, std::vector<float>{}};

  const int nx = array.shape.x;
  const int ny = array.shape.y;

  if (nx < 1 || ny < 1) return {std::vector<float>{}, std::vector<float>{}};

  // 1D case handling
  if (nx == 1 && ny == 1)
  {
    float cx = 0.5f * (bbox.x + bbox.y);
    float cy = 0.5f * (bbox.z + bbox.w);
    return {std::vector<float>(count, cx), std::vector<float>(count, cy)};
  }

  // Row marginal CDF (R: prefix sums of row sums)
  std::vector<float> r_cdf(ny);
  float              running_row_sum = 0.f;
  for (int y = 0; y < ny; ++y)
  {
    float row_sum = 0.f;
    for (int x = 0; x < nx; ++x)
      row_sum += std::max(0.f, array(x, y));
    running_row_sum += row_sum;
    r_cdf[y] = running_row_sum;
  }

  const float total_weight = r_cdf.back();

  // If array is all zero, fallback to uniform random points within bbox
  if (total_weight <= 0.f)
    return random_points(count, seed, PointSamplingMethod::RND_RANDOM, bbox);

  std::mt19937                          rng(seed);
  std::uniform_real_distribution<float> dist(0.f, 1.f);

  std::vector<float> x_coords(count);
  std::vector<float> y_coords(count);
  std::vector<float> c_cdf(nx);

  const float bbox_dx = bbox.y - bbox.x;
  const float bbox_dy = bbox.w - bbox.z;
  const float denom_y = (ny > 1) ? static_cast<float>(ny - 1) : 1.f;
  const float denom_x = (nx > 1) ? static_cast<float>(nx - 1) : 1.f;

  for (size_t k = 0; k < count; ++k)
  {
    // Row sampling
    float u = dist(rng) * total_weight;
    auto  it_y = std::upper_bound(r_cdf.begin(), r_cdf.end(), u);
    int   i = static_cast<int>(std::distance(r_cdf.begin(), it_y));
    if (i >= ny) i = ny - 1;

    float lo = (i > 0) ? r_cdf[i - 1] : 0.f;
    float hi = r_cdf[i];
    float delta_r = hi - lo;
    float y0 = (delta_r > 0.f) ? (static_cast<float>(i) + (u - lo) / delta_r)
                               : static_cast<float>(i);

    // Conditional CDF over x, blended between rows
    int   r0 = std::clamp(i, 0, ny - 1);
    int   r1 = std::min(i + 1, ny - 1);
    float t = std::clamp(y0 - static_cast<float>(i), 0.f, 1.f);

    float acc = 0.f;
    for (int j = 0; j < nx; ++j)
    {
      float v0 = std::max(0.f, array(j, r0));
      float v1 = std::max(0.f, array(j, r1));
      acc += (1.f - t) * v0 + t * v1;
      c_cdf[j] = acc;
    }

    float row_weight = c_cdf.back();
    float x0 = 0.f;
    if (row_weight > 0.f)
    {
      float u2 = dist(rng) * row_weight;
      auto  it_x = std::upper_bound(c_cdf.begin(), c_cdf.end(), u2);
      int   j = static_cast<int>(std::distance(c_cdf.begin(), it_x));
      if (j >= nx) j = nx - 1;

      float lo2 = (j > 0) ? c_cdf[j - 1] : 0.f;
      float hi2 = c_cdf[j];
      float delta_c = hi2 - lo2;
      x0 = (delta_c > 0.f) ? (static_cast<float>(j) + (u2 - lo2) / delta_c)
                           : static_cast<float>(j);
    }
    else
    {
      x0 = dist(rng) * static_cast<float>(nx - 1);
    }

    // Map from grid coordinate [0, nx-1] and [0, ny-1] to bbox
    float norm_x = (nx > 1) ? std::clamp(x0 / denom_x, 0.f, 1.f) : 0.5f;
    float norm_y = (ny > 1) ? std::clamp(y0 / denom_y, 0.f, 1.f) : 0.5f;

    x_coords[k] = bbox.x + norm_x * bbox_dx;
    y_coords[k] = bbox.z + norm_y * bbox_dy;
  }

  return {std::move(x_coords), std::move(y_coords)};
}

std::array<std::vector<float>, 2> random_points_distance(float         min_dist,
                                                         std::uint32_t seed,
                                                         const glm::vec4 &bbox)
{
  auto   ranges = bbox_to_ranges2d(bbox);
  size_t count = helper_estimate_count(bbox, min_dist);

  auto points = ps::poisson_disk_sampling_uniform(count,
                                                  ranges,
                                                  min_dist,
                                                  seed);
  return ps::split_by_dimension(points);
}

std::array<std::vector<float>, 2> random_points_distance(float         min_dist,
                                                         float         max_dist,
                                                         const Array  &density,
                                                         std::uint32_t seed,
                                                         const glm::vec4 &bbox)
{
  if (!validate_non_empty(density))
    return {std::vector<float>{}, std::vector<float>{}};

  auto ranges = bbox_to_ranges2d(bbox);

  // convert density [0, 1] to scale (when density = 0, enforce
  // max_dist, when density = 1, enforce min_dist)
  Array scale = (1.f - density) * (max_dist / min_dist - 1.f) + 1.f;
  auto  scale_fct = make_pointwise_function_from_array(scale, bbox);

  // estimate a maximum count using the minimum distance
  size_t count = helper_estimate_count(bbox, min_dist);

  auto points = ps::poisson_disk_sampling<float, 2>(count,
                                                    ranges,
                                                    min_dist,
                                                    scale_fct,
                                                    seed);
  return ps::split_by_dimension(points);
}

std::array<std::vector<float>, 2> random_points_distance_power_law(
    float            dist_min,
    float            dist_max,
    float            alpha,
    std::uint32_t    seed,
    const glm::vec4 &bbox)
{
  auto   ranges = bbox_to_ranges2d(bbox);
  size_t count = helper_estimate_count(bbox, dist_min);

  auto points = ps::poisson_disk_sampling_power_law<float, 2>(count,
                                                              dist_min,
                                                              dist_max,
                                                              alpha,
                                                              ranges,
                                                              seed);
  return ps::split_by_dimension(points);
}

std::array<std::vector<float>, 2> random_points_distance_weibull(
    float            dist_min,
    float            lambda,
    float            k,
    std::uint32_t    seed,
    const glm::vec4 &bbox)
{
  auto   ranges = bbox_to_ranges2d(bbox);
  size_t count = helper_estimate_count(bbox, dist_min);

  auto points = ps::poisson_disk_sampling_weibull<float, 2>(count,
                                                            lambda,
                                                            k,
                                                            dist_min,
                                                            ranges,
                                                            seed);
  return ps::split_by_dimension(points);
}

std::array<std::vector<float>, 2> random_points_jittered(
    size_t           count,
    const glm::vec2 &jitter_amount,
    const glm::vec2 &stagger_ratio,
    std::uint32_t    seed,
    const glm::vec4 &bbox)
{
  auto                 ranges = bbox_to_ranges2d(bbox);
  std::array<float, 2> jt = {jitter_amount.x, jitter_amount.y};
  std::array<float, 2> sr = {stagger_ratio.x, stagger_ratio.y};

  auto points = ps::jittered_grid(count, ranges, jt, sr, seed);
  return ps::split_by_dimension(points);
}

} // namespace hmap
