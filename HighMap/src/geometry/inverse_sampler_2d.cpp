/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>

#include "highmap/geometry/point_sampling.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

InverseSampler2D::InverseSampler2D(const Array     &density,
                                   std::uint32_t    seed,
                                   const glm::vec4 &bbox)
    : density_(density), bbox_(bbox), rng_(seed)
{
  this->rebuild_cdf();
}

const glm::vec4 &InverseSampler2D::get_bbox() const
{
  return this->bbox_;
}

const Array &InverseSampler2D::get_density() const
{
  return this->density_;
}

float InverseSampler2D::get_total_weight() const
{
  if (this->r_cdf_.empty()) return 0.f;
  return this->r_cdf_.back();
}

void InverseSampler2D::rebuild_cdf()
{
  const int ny = this->density_.shape.y;
  const int nx = this->density_.shape.x;

  if (nx <= 0 || ny <= 0)
  {
    this->row_sums_.clear();
    this->r_cdf_.clear();
    return;
  }

  this->row_sums_.resize(ny);
  this->r_cdf_.resize(ny);

  float running_sum = 0.f;
  for (int y = 0; y < ny; ++y)
  {
    float row_sum = 0.f;
    for (int x = 0; x < nx; ++x)
      row_sum += std::max(0.f, this->density_(x, y));
    this->row_sums_[y] = row_sum;
    running_sum += row_sum;
    this->r_cdf_[y] = running_sum;
  }
}

void InverseSampler2D::reset(const Array     &density,
                             std::uint32_t    seed,
                             const glm::vec4 &bbox)
{
  this->density_ = density;
  this->bbox_ = bbox;
  this->rng_.seed(seed);
  this->rebuild_cdf();
}

glm::vec2 InverseSampler2D::sample()
{
  return this->sample(this->dist_(this->rng_), this->dist_(this->rng_));
}

glm::vec2 InverseSampler2D::sample(float u_in, float u2_in)
{
  const int nx = this->density_.shape.x;
  const int ny = this->density_.shape.y;

  if (nx <= 0 || ny <= 0)
    return glm::vec2(0.5f * (this->bbox_.x + this->bbox_.y),
                     0.5f * (this->bbox_.z + this->bbox_.w));

  if (nx == 1 && ny == 1)
    return glm::vec2(0.5f * (this->bbox_.x + this->bbox_.y),
                     0.5f * (this->bbox_.z + this->bbox_.w));

  const float total_weight = this->get_total_weight();
  const float bbox_dx = this->bbox_.y - this->bbox_.x;
  const float bbox_dy = this->bbox_.w - this->bbox_.z;
  const float denom_y = (ny > 1) ? static_cast<float>(ny - 1) : 1.f;
  const float denom_x = (nx > 1) ? static_cast<float>(nx - 1) : 1.f;

  if (total_weight <= 0.f)
  {
    float norm_x = std::clamp(u_in, 0.f, 1.f);
    float norm_y = std::clamp(u2_in, 0.f, 1.f);
    return glm::vec2(this->bbox_.x + norm_x * bbox_dx,
                     this->bbox_.z + norm_y * bbox_dy);
  }

  // Row sampling
  float u = std::clamp(u_in, 0.f, 1.f) * total_weight;
  auto  it_y = std::upper_bound(this->r_cdf_.begin(), this->r_cdf_.end(), u);
  int   i = static_cast<int>(std::distance(this->r_cdf_.begin(), it_y));
  if (i >= ny) i = ny - 1;

  float lo = (i > 0) ? this->r_cdf_[i - 1] : 0.f;
  float hi = this->r_cdf_[i];
  float delta_r = hi - lo;
  float y0 = (delta_r > 0.f) ? (static_cast<float>(i) + (u - lo) / delta_r)
                             : static_cast<float>(i);

  // Conditional CDF over x, blended between rows
  int   r0 = std::clamp(i, 0, ny - 1);
  int   r1 = std::min(i + 1, ny - 1);
  float t = std::clamp(y0 - static_cast<float>(i), 0.f, 1.f);

  std::vector<float> c_cdf(nx);
  float              acc = 0.f;
  for (int j = 0; j < nx; ++j)
  {
    float v0 = std::max(0.f, this->density_(j, r0));
    float v1 = std::max(0.f, this->density_(j, r1));
    acc += (1.f - t) * v0 + t * v1;
    c_cdf[j] = acc;
  }

  float row_weight = c_cdf.back();
  float x0 = 0.f;
  if (row_weight > 0.f)
  {
    float u2 = std::clamp(u2_in, 0.f, 1.f) * row_weight;
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
    x0 = u2_in * static_cast<float>(nx - 1);
  }

  float norm_x = (nx > 1) ? std::clamp(x0 / denom_x, 0.f, 1.f) : 0.5f;
  float norm_y = (ny > 1) ? std::clamp(y0 / denom_y, 0.f, 1.f) : 0.5f;

  return glm::vec2(this->bbox_.x + norm_x * bbox_dx,
                   this->bbox_.z + norm_y * bbox_dy);
}

std::array<std::vector<float>, 2> InverseSampler2D::sample(size_t count)
{
  std::vector<float> x_coords(count);
  std::vector<float> y_coords(count);

  for (size_t k = 0; k < count; ++k)
  {
    glm::vec2 p = this->sample();
    x_coords[k] = p.x;
    y_coords[k] = p.y;
  }

  return {std::move(x_coords), std::move(y_coords)};
}

glm::vec2 InverseSampler2D::sample(const glm::vec4 &restricted_bbox)
{
  return this->sample(restricted_bbox,
                      this->dist_(this->rng_),
                      this->dist_(this->rng_));
}

glm::vec2 InverseSampler2D::sample(const glm::vec4 &restricted_bbox,
                                   float            u_in,
                                   float            u2_in)
{
  const int nx = this->density_.shape.x;
  const int ny = this->density_.shape.y;

  if (nx <= 0 || ny <= 0)
    return glm::vec2(0.5f * (restricted_bbox.x + restricted_bbox.y),
                     0.5f * (restricted_bbox.z + restricted_bbox.w));

  const float bbox_dx = this->bbox_.y - this->bbox_.x;
  const float bbox_dy = this->bbox_.w - this->bbox_.z;

  if (bbox_dx <= 1e-7f || bbox_dy <= 1e-7f)
    return glm::vec2(this->bbox_.x, this->bbox_.z);

  // clamp restricted bbox to the sampler domain
  float xmin = std::clamp(std::min(restricted_bbox.x, restricted_bbox.y),
                          this->bbox_.x,
                          this->bbox_.y);
  float xmax = std::clamp(std::max(restricted_bbox.x, restricted_bbox.y),
                          this->bbox_.x,
                          this->bbox_.y);
  float ymin = std::clamp(std::min(restricted_bbox.z, restricted_bbox.w),
                          this->bbox_.z,
                          this->bbox_.w);
  float ymax = std::clamp(std::max(restricted_bbox.z, restricted_bbox.w),
                          this->bbox_.z,
                          this->bbox_.w);

  if (xmax <= xmin || ymax <= ymin) return glm::vec2(xmin, ymin);

  const float denom_x = (nx > 1) ? static_cast<float>(nx - 1) : 1.f;
  const float denom_y = (ny > 1) ? static_cast<float>(ny - 1) : 1.f;

  float norm_xmin = (xmin - this->bbox_.x) / bbox_dx;
  float norm_xmax = (xmax - this->bbox_.x) / bbox_dx;
  float norm_ymin = (ymin - this->bbox_.z) / bbox_dy;
  float norm_ymax = (ymax - this->bbox_.z) / bbox_dy;

  float x0_min = norm_xmin * denom_x;
  float x0_max = norm_xmax * denom_x;
  float y0_min = norm_ymin * denom_y;
  float y0_max = norm_ymax * denom_y;

  int jmin = std::clamp(static_cast<int>(std::floor(x0_min)), 0, nx - 1);
  int jmax = std::clamp(static_cast<int>(std::ceil(x0_max)), 0, nx - 1);
  int imin = std::clamp(static_cast<int>(std::floor(y0_min)), 0, ny - 1);
  int imax = std::clamp(static_cast<int>(std::ceil(y0_max)), 0, ny - 1);

  int num_rows = imax - imin + 1;
  int num_cols = jmax - jmin + 1;

  std::vector<float> local_r_cdf(num_rows, 0.f);

  float running_r = 0.f;
  for (int r = imin; r <= imax; ++r)
  {
    float row_sum = 0.f;
    for (int c = jmin; c <= jmax; ++c)
    {
      row_sum += std::max(0.f, this->density_(c, r));
    }
    int idx = r - imin;
    running_r += row_sum;
    local_r_cdf[idx] = running_r;
  }

  float total_sub_weight = local_r_cdf.empty() ? 0.f : local_r_cdf.back();

  if (total_sub_weight <= 0.f)
  {
    float u1_c = std::clamp(u_in, 0.f, 1.f);
    float u2_c = std::clamp(u2_in, 0.f, 1.f);
    return glm::vec2(xmin + u1_c * (xmax - xmin), ymin + u2_c * (ymax - ymin));
  }

  // Row sampling within restricted region
  float u = std::clamp(u_in, 0.f, 1.f) * total_sub_weight;
  auto  it_y = std::upper_bound(local_r_cdf.begin(), local_r_cdf.end(), u);
  int   local_i = static_cast<int>(std::distance(local_r_cdf.begin(), it_y));
  if (local_i >= num_rows) local_i = num_rows - 1;

  int global_i = imin + local_i;

  float lo = (local_i > 0) ? local_r_cdf[local_i - 1] : 0.f;
  float hi = local_r_cdf[local_i];
  float delta_r = hi - lo;
  float y0 = (delta_r > 0.f)
                 ? (static_cast<float>(global_i) + (u - lo) / delta_r)
                 : static_cast<float>(global_i);

  y0 = std::clamp(y0, y0_min, y0_max);

  // Column sampling within restricted region
  int   r0 = std::clamp(static_cast<int>(std::floor(y0)), 0, ny - 1);
  int   r1 = std::min(r0 + 1, ny - 1);
  float t = std::clamp(y0 - static_cast<float>(r0), 0.f, 1.f);

  std::vector<float> c_cdf(num_cols, 0.f);
  float              acc = 0.f;
  for (int c = jmin; c <= jmax; ++c)
  {
    float v0 = std::max(0.f, this->density_(c, r0));
    float v1 = std::max(0.f, this->density_(c, r1));
    acc += (1.f - t) * v0 + t * v1;
    c_cdf[c - jmin] = acc;
  }

  float row_weight = c_cdf.empty() ? 0.f : c_cdf.back();
  float x0 = 0.f;
  if (row_weight > 0.f)
  {
    float u2 = std::clamp(u2_in, 0.f, 1.f) * row_weight;
    auto  it_x = std::upper_bound(c_cdf.begin(), c_cdf.end(), u2);
    int   local_j = static_cast<int>(std::distance(c_cdf.begin(), it_x));
    if (local_j >= num_cols) local_j = num_cols - 1;

    int global_j = jmin + local_j;

    float lo2 = (local_j > 0) ? c_cdf[local_j - 1] : 0.f;
    float hi2 = c_cdf[local_j];
    float delta_c = hi2 - lo2;
    x0 = (delta_c > 0.f) ? (static_cast<float>(global_j) + (u2 - lo2) / delta_c)
                         : static_cast<float>(global_j);
  }
  else
  {
    x0 = x0_min + u2_in * (x0_max - x0_min);
  }

  x0 = std::clamp(x0, x0_min, x0_max);

  float norm_x = (nx > 1) ? std::clamp(x0 / denom_x, 0.f, 1.f) : 0.5f;
  float norm_y = (ny > 1) ? std::clamp(y0 / denom_y, 0.f, 1.f) : 0.5f;

  return glm::vec2(this->bbox_.x + norm_x * bbox_dx,
                   this->bbox_.z + norm_y * bbox_dy);
}

void InverseSampler2D::set_bbox(const glm::vec4 &bbox)
{
  this->bbox_ = bbox;
}

void InverseSampler2D::set_seed(std::uint32_t seed)
{
  this->rng_.seed(seed);
}

void InverseSampler2D::update_region(const Array      &updated_array,
                                     const glm::ivec2 &range_i,
                                     const glm::ivec2 &range_j)
{
  if (!validate_non_empty(this->density_)) return;

  const int nx = this->density_.shape.x;
  const int ny = this->density_.shape.y;

  int imin = std::clamp(std::min(range_i.x, range_i.y), 0, nx - 1);
  int imax = std::clamp(std::max(range_i.x, range_i.y), 0, nx - 1);
  int jmin = std::clamp(std::min(range_j.x, range_j.y), 0, ny - 1);
  int jmax = std::clamp(std::max(range_j.x, range_j.y), 0, ny - 1);

  if (imin > imax || jmin > jmax) return;

  for (int y = jmin; y <= jmax; ++y)
  {
    for (int x = imin; x <= imax; ++x)
    {
      if (x < updated_array.shape.x && y < updated_array.shape.y)
        this->density_(x, y) = updated_array(x, y);
    }
    // Recompute row sum for row y
    float row_sum = 0.f;
    for (int x = 0; x < nx; ++x)
      row_sum += std::max(0.f, this->density_(x, y));
    this->row_sums_[y] = row_sum;
  }

  this->update_row_cdf_from(jmin);
}

void InverseSampler2D::update_region(
    const std::function<float(int i, int j, float current)> &updater,
    const glm::ivec2                                        &range_i,
    const glm::ivec2                                        &range_j)
{
  if (!validate_non_empty(this->density_) || !updater) return;

  const int nx = this->density_.shape.x;
  const int ny = this->density_.shape.y;

  int imin = std::clamp(std::min(range_i.x, range_i.y), 0, nx - 1);
  int imax = std::clamp(std::max(range_i.x, range_i.y), 0, nx - 1);
  int jmin = std::clamp(std::min(range_j.x, range_j.y), 0, ny - 1);
  int jmax = std::clamp(std::max(range_j.x, range_j.y), 0, ny - 1);

  if (imin > imax || jmin > jmax) return;

  for (int y = jmin; y <= jmax; ++y)
  {
    for (int x = imin; x <= imax; ++x)
    {
      float cur = this->density_(x, y);
      this->density_(x, y) = updater(x, y, cur);
    }
    // Recompute row sum for row y
    float row_sum = 0.f;
    for (int x = 0; x < nx; ++x)
      row_sum += std::max(0.f, this->density_(x, y));
    this->row_sums_[y] = row_sum;
  }

  this->update_row_cdf_from(jmin);
}

void InverseSampler2D::update_row_cdf_from(int jmin)
{
  const int ny = static_cast<int>(this->r_cdf_.size());
  if (ny == 0) return;

  float running_sum = (jmin > 0) ? this->r_cdf_[jmin - 1] : 0.f;
  for (int y = jmin; y < ny; ++y)
  {
    running_sum += this->row_sums_[y];
    this->r_cdf_[y] = running_sum;
  }
}

} // namespace hmap
