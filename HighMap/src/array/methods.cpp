/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include <opencv2/core.hpp>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/interpolate/interpolate_array.hpp"
#include "highmap/range.hpp"
#include "highmap/vectors.hpp"

namespace hmap
{

void Array::argmax(float &max, int &im, int &jm) const
{
  max = -std::numeric_limits<float>::infinity();
  im = -1;
  jm = -1;

  if (!validate_non_empty(*this)) return;

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      float v = (*this)(i, j);
      if (v > max)
      {
        max = v;
        im = i;
        jm = j;
      }
    }
}

std::vector<float> Array::col_to_vector(int j)
{
  std::vector<float> vec(this->shape.x);
  for (int i = 0; i < this->shape.x; i++)
    vec[i] = (*this)(i, j);
  return vec;
}

void Array::depose_amount_bilinear_at(int   i,
                                      int   j,
                                      float u,
                                      float v,
                                      float amount)
{
  (*this)(i, j) += amount * (1 - u) * (1 - v);
  (*this)(i + 1, j) += amount * u * (1 - v);
  (*this)(i, j + 1) += amount * (1 - u) * v;
  (*this)(i + 1, j + 1) += amount * u * v;
}

void Array::depose_amount_kernel_bilinear_at(int   i,
                                             int   j,
                                             float u,
                                             float v,
                                             int   ir,
                                             float amount)
{
  Array kernel = Array(glm::ivec2(2 * ir + 1, 2 * ir + 1));

  // compute kernel first
  for (int p = -ir; p < ir + 1; p++)
  {
    for (int q = -ir; q < ir + 1; q++)
    {
      float x = (float)p - u;
      float y = (float)q - v;
      float r = std::max(0.f, 1.f - std::hypot(x, y));
      kernel(p + ir, q + ir) = r;
    }
  }
  kernel.normalize();

  // perform deposition
  this->depose_amount_kernel_at(i, j, kernel, amount);
}

void Array::depose_amount_kernel_at(int          i,
                                    int          j,
                                    const Array &kernel,
                                    float        amount)
{
  const int ir = (kernel.shape.x - 1) / 2;
  const int jr = (kernel.shape.y - 1) / 2;

  for (int p = 0; p < kernel.shape.x; p++)
  {
    for (int q = 0; q < kernel.shape.y; q++)
    {
      (*this)(i + p - ir, j + q - jr) += amount * kernel(p, q);
    }
  }
}

void Array::dump(const std::string &fname) const
{
  this->infos(fname);
  this->to_png_grayscale(fname, CV_16U);
}

void Array::dump_histogram(int bins, const std::string &msg) const
{
  std::cout << "Array: " << msg << "\n";
  std::cout << make_histogram(this->vector, bins, 6);
}

Array Array::extract_slice(glm::ivec4 idx) const
{
  if (!validate_slice(*this, idx)) return Array();

  Array array_out = Array(glm::ivec2(idx.y - idx.x, idx.w - idx.z));

  for (int j = idx.z; j < idx.w; j++)
    for (int i = idx.x; i < idx.y; i++)
      array_out(i - idx.x, j - idx.z) = (*this)(i, j);

  return array_out;
}

Array Array::extract_slice(int i1, int i2, int j1, int j2) const
{
  glm::ivec4 idx(i1, i2, j1, j2);
  return this->extract_slice(idx);
}

float Array::get_gradient_x_at(int i, int j) const
{
  return 0.5f * ((*this)(i + 1, j) - (*this)(i - 1, j));
}

float Array::get_gradient_y_at(int i, int j) const
{
  return 0.5f * ((*this)(i, j + 1) - (*this)(i, j - 1));
}

float Array::get_gradient_x_bilinear_at(int i, int j, float u, float v) const
{
  float f00 = (*this)(i, j) - (*this)(i - 1, j);
  float f10 = (*this)(i + 1, j) - (*this)(i, j);
  float f01 = (*this)(i, j + 1) - (*this)(i - 1, j + 1);
  float f11 = (*this)(i + 1, j + 1) - (*this)(i, j + 1);

  float a10 = f10 - f00;
  float a01 = f01 - f00;
  float a11 = f11 - f10 - f01 + f00;

  return f00 + a10 * u + a01 * v + a11 * u * v;
}

float Array::get_gradient_y_bilinear_at(int i, int j, float u, float v) const
{
  float f00 = (*this)(i, j) - (*this)(i, j - 1);
  float f10 = (*this)(i + 1, j) - (*this)(i + 1, j - 1);
  float f01 = (*this)(i, j + 1) - (*this)(i, j);
  float f11 = (*this)(i + 1, j + 1) - (*this)(i + 1, j);

  float a10 = f10 - f00;
  float a01 = f01 - f00;
  float a11 = f11 - f10 - f01 + f00;

  return f00 + a10 * u + a01 * v + a11 * u * v;
}

glm::vec3 Array::get_normal_at(int i, int j) const
{
  glm::vec3 normal;

  normal.x = -this->get_gradient_x_at(i, j);
  normal.y = -this->get_gradient_y_at(i, j);
  normal.z = 1.f;

  float norm = std::hypot(normal.x, normal.y, normal.z);
  normal /= norm;

  return normal;
}

size_t Array::get_sizeof() const
{
  return sizeof(float) * this->vector.size();
}

glm::vec2 Array::normalization_coeff(float vmin, float vmax) const
{
  float a = 0.f;
  float b = 0.f;
  if (vmin != vmax)
  {
    a = 1.f / (vmax - vmin);
    b = -vmin / (vmax - vmin);
  }
  return glm::vec2(a, b);
}

float Array::get_value_bicubic_at(int i, int j, float u, float v) const
{
  float arr[4][4];

  // Get the 4x4 surrounding grid points
  for (int n = -1; n <= 2; ++n)
    for (int m = -1; m <= 2; ++m)
    {
      int ip = std::clamp(i + m, 0, this->shape.x - 1);
      int jp = std::clamp(j + n, 0, this->shape.y - 1);
      arr[m + 1][n + 1] = (*this)(ip, jp);
    }

  // interpolate in the x direction
  float col_results[4];
  for (int i = 0; i < 4; ++i)
    col_results[i] = cubic_interpolate(arr[i], v);

  // interpolate in the y direction
  return cubic_interpolate(col_results, u);
}

float Array::get_value_bilinear_at(int i, int j, float u, float v) const
{
  if (i == this->shape.x - 1)
  {
    i = this->shape.x - 2;
    u = 1.f;
  }

  if (j == this->shape.y - 1)
  {
    j = this->shape.y - 2;
    v = 1.f;
  }

  float a10 = (*this)(i + 1, j) - (*this)(i, j);
  float a01 = (*this)(i, j + 1) - (*this)(i, j);
  float a11 = (*this)(i + 1, j + 1) - (*this)(i + 1, j) - (*this)(i, j + 1) +
              (*this)(i, j);

  return (*this)(i, j) + a10 * u + a01 * v + a11 * u * v;
}

float Array::get_value_bilinear(float x, float y, glm::vec4 bbox) const
{
  float dx = bbox.y - bbox.x;
  float dy = bbox.w - bbox.z;

  float u_norm = (dx != 0.f) ? (x - bbox.x) / dx : 0.f;
  float v_norm = (dy != 0.f) ? (y - bbox.z) / dy : 0.f;

  float x_scaled = std::clamp(u_norm, 0.f, 1.f) * (this->shape.x - 1);
  float y_scaled = std::clamp(v_norm, 0.f, 1.f) * (this->shape.y - 1);

  int i = static_cast<int>(x_scaled);
  int j = static_cast<int>(y_scaled);

  i = std::clamp(i, 0, this->shape.x - 1);
  j = std::clamp(j, 0, this->shape.y - 1);

  float u = x_scaled - static_cast<float>(i);
  float v = y_scaled - static_cast<float>(j);

  return this->get_value_bilinear_at(i, j, u, v);
}

float Array::get_value_nearest(float x, float y, glm::vec4 bbox) const
{
  int i = (int)(std::clamp((x - bbox.x) / (bbox.y - bbox.x), 0.f, 1.f) *
                (this->shape.x - 1));
  int j = (int)(std::clamp((y - bbox.z) / (bbox.w - bbox.z), 0.f, 1.f) *
                (this->shape.y - 1));
  return (*this)(i, j);
}

int Array::linear_index(int i, int j) const
{
  return j * this->shape.x + i;
}

glm::ivec2 Array::linear_index_reverse(int k) const
{
  glm::ivec2 ij;
  ij.y = std::floor(k / shape.x);
  ij.x = k - ij.y * shape.x;
  return ij;
}

float Array::max(ICommunicator &comm) const
{
  if (!validate_non_empty(*this)) return 0.f;
  float local_max = *std::max_element(this->vector.begin(), this->vector.end());
  return comm.all_reduce(local_max, ReduceOp::Max);
}

float Array::mean(ICommunicator &comm) const
{
  if (!validate_non_empty(*this)) return 0.f;

  float total_sum = this->sum(comm);
  float total_size = static_cast<float>(this->size());

  if (comm.size() > 1)
  {
    total_size = comm.all_reduce(total_size, ReduceOp::Sum);
  }

  return (total_size > 0.f) ? (total_sum / total_size) : 0.f;
}

float Array::median() const
{
  if (!validate_non_empty(*this)) return 0.f;
  return compute_median(this->vector);
}

float Array::min(ICommunicator &comm) const
{
  if (!validate_non_empty(*this)) return 0.f;
  float local_min = *std::min_element(this->vector.begin(), this->vector.end());
  return comm.all_reduce(local_min, ReduceOp::Min);
}

void Array::normalize(ICommunicator &comm)
{
  if (!validate_non_empty(*this)) return;

  float total_sum = this->sum(comm);

  if (total_sum != 0.f)
  {
    for (float &v : this->vector)
      v /= total_sum;
  }
}

float Array::ptp(ICommunicator &comm) const
{
  if (!validate_non_empty(*this)) return 0.f;
  return this->max(comm) - this->min(comm);
}

glm::vec2 Array::range(ICommunicator &comm) const
{
  if (!validate_non_empty(*this)) return glm::vec2(0.f, 0.f);

  auto [min_it, max_it] = std::minmax_element(this->vector.begin(),
                                              this->vector.end());
  const float local_min = *min_it;
  const float local_max = *max_it;

  float global_min = comm.all_reduce(local_min, ReduceOp::Min);
  float global_max = comm.all_reduce(local_max, ReduceOp::Max);

  return glm::vec2(global_min, global_max);
}

glm::vec2 Array::range_percentile(float  p_low,
                                  float  p_high,
                                  size_t bins,
                                  bool   ignore_zeros) const
{
  if (!validate_non_empty(*this) ||
      !validate_parameter_range(glm::vec2(p_low, p_high),
                                0.f,
                                1.f,
                                "Percentile range"))
    return glm::vec2(0.f, 0.f);

  const std::vector<float> &data = this->vector;

  float min_val = 0.f;
  float max_val = 0.f;

  if (ignore_zeros)
  {
    bool found = false;
    for (float v : data)
    {
      if (v != 0.f)
      {
        if (!found)
        {
          min_val = v;
          max_val = v;
          found = true;
        }
        else
        {
          min_val = std::min(min_val, v);
          max_val = std::max(max_val, v);
        }
      }
    }

    if (!found) return glm::vec2(0.f, 0.f);
  }
  else
  {
    const glm::vec2 range = this->range();
    min_val = range.x;
    max_val = range.y;
  }

  if (min_val == max_val) return glm::vec2(min_val, max_val);

  std::vector<size_t> hist(bins, 0);
  size_t              total = 0;

  // fill histogram
  for (float v : data)
  {
    if (ignore_zeros && v == 0.f) continue;

    size_t idx = std::min<size_t>(
        bins - 1,
        size_t((v - min_val) / (max_val - min_val) * bins));
    hist[idx]++;
    ++total;
  }

  if (total == 0) return glm::vec2(0.f, 0.f);

  // accumulate
  size_t low_target = size_t(p_low * total);
  size_t high_target = size_t(p_high * total);

  size_t acc = 0;
  float  low = min_val;
  float  high = max_val;

  for (size_t i = 0; i < bins; ++i)
  {
    acc += hist[i];

    if (acc >= low_target && low == min_val)
      low = min_val + (max_val - min_val) * (float(i) / bins);

    if (acc >= high_target)
    {
      high = min_val + (max_val - min_val) * (float(i) / bins);
      break;
    }
  }

  return {low, high};
}

Array Array::remapped() const
{
  Array out = *this;
  remap(out);
  return out;
}

Array Array::resample_to_shape(glm::ivec2 new_shape) const
{
  if (!validate_shape(new_shape)) return Array();
  if (this->shape == new_shape) return *this;

  return this->resample_to_shape_bilinear(new_shape);
}

Array Array::resample_to_shape_bicubic(glm::ivec2 new_shape) const
{
  if (!validate_shape(new_shape)) return Array();
  if (this->shape == new_shape) return *this;

  Array array_out = Array(new_shape);
  interpolate_array_bicubic(*this, array_out);

  return array_out;
}

Array Array::resample_to_shape_bilinear(glm::ivec2 new_shape) const
{
  if (!validate_shape(new_shape)) return Array();
  if (this->shape == new_shape) return *this;

  Array array_out = Array(new_shape);
  interpolate_array_bilinear(*this, array_out);

  return array_out;
}

Array Array::resample_to_shape_nearest(glm::ivec2 new_shape) const
{
  if (!validate_shape(new_shape)) return Array();
  if (this->shape == new_shape) return *this;

  Array array_out = Array(new_shape);
  interpolate_array_nearest(*this, array_out);

  return array_out;
}

std::vector<float> Array::row_to_vector(int i)
{
  std::vector<float> vec(this->shape.y);
  for (int j = 0; j < this->shape.y; j++)
    vec[j] = (*this)(i, j);
  return vec;
}

void Array::set_slice(glm::ivec4 idx, float value)
{
  if (!validate_slice(*this, idx)) return;

  for (int i = idx.x; i < idx.y; i++)
    for (int j = idx.z; j < idx.w; j++)
      (*this)(i, j) = value;
}

void Array::set_slice(glm::ivec4 idx, const Array &array)
{
  if (!validate_slice_for_array(*this, idx, array)) return;

  for (int i = idx.x; i < idx.y; i++)
    for (int j = idx.z; j < idx.w; j++)
      (*this)(i, j) = array(i - idx.x, j - idx.z);
}

int Array::size() const
{
  return this->shape.x * this->shape.y;
}

float Array::std(ICommunicator &comm) const
{
  if (!validate_non_empty(*this)) return 0.f;

  float global_mean = this->mean(comm);
  float local_sq_diff_sum = 0.f;

  for (float v : this->vector)
  {
    float diff = v - global_mean;
    local_sq_diff_sum += diff * diff;
  }

  float total_sq_diff_sum = comm.all_reduce(local_sq_diff_sum, ReduceOp::Sum);
  float total_size = static_cast<float>(this->size());

  if (comm.size() > 1)
  {
    total_size = comm.all_reduce(total_size, ReduceOp::Sum);
  }

  return (total_size > 0.f) ? std::sqrt(total_sq_diff_sum / total_size) : 0.f;
}

float Array::sum(ICommunicator &comm) const
{
  if (!validate_non_empty(*this)) return 0.f;

  float local_sum = std::accumulate(this->vector.begin(),
                                    this->vector.end(),
                                    0.f);
  return comm.all_reduce(local_sum, ReduceOp::Sum);
}

std::vector<float> Array::unique_values(ICommunicator &comm) const
{
  if (!validate_non_empty(*this)) return {};
  std::vector<float> v = this->vector;
  vector_unique_values(v);
  return comm.all_union(v);
}

} // namespace hmap
