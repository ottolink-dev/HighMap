/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <filesystem>
#include <initializer_list>
#include <string>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/export/image.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

Array::Array()
{
}

Array::Array(glm::ivec2 shape) : shape(shape)
{
  if (!validate_shape(shape))
  {
    this->shape = glm::ivec2(0, 0);
    return;
  }
  this->vector.resize(this->shape.x * this->shape.y);
}

Array::Array(glm::ivec2 shape, float value) : shape(shape)
{
  if (!validate_shape(shape))
  {
    this->shape = glm::ivec2(0, 0);
    return;
  }
  this->vector.resize(this->shape.x * this->shape.y);
  std::fill(this->vector.begin(), this->vector.end(), value);
}

Array::Array(const std::string &filename, bool flip_j)
{
  bool remap = true;

  std::filesystem::path file_path(filename);
  std::filesystem::path ext = file_path.extension();
  if (ext.string() == ".exr") remap = false;

  *this = read_to_array(filename, flip_j, remap);
}

Array::Array(const std::vector<std::vector<float>> &data)
{
  if (data.empty() || data[0].empty())
  {
    (void)validate_shape(glm::ivec2(0, 0));
    return;
  }

  this->shape = glm::ivec2(data.size(), data[0].size());
  *this = Array(shape);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      (*this)(i, j) = data[i][j];
}

Array::Array(const std::initializer_list<std::initializer_list<float>> &data)
{
  std::vector<std::vector<float>> tmp;
  tmp.reserve(data.size());

  for (const auto &row : data)
    tmp.emplace_back(row); // vector<float>(initializer_list<float>)

  *this = Array(tmp);
}

glm::ivec2 Array::get_shape()
{
  return shape;
}

std::vector<float> Array::get_vector() const
{
  return this->vector;
}

void Array::set_shape(glm::ivec2 new_shape)
{
  if (!validate_shape(new_shape))
  {
    this->shape = glm::ivec2(0, 0);
    this->vector.clear();
    return;
  }
  this->shape = new_shape;
  this->vector.resize(this->shape.x * this->shape.y);
}

Array &Array::operator=(const float value)
{
  std::fill(this->vector.begin(), this->vector.end(), value);
  return *this;
}

Array &Array::operator*=(const float value)
{
  std::transform(this->vector.begin(),
                 this->vector.end(),
                 this->vector.begin(),
                 [&value](float v) { return v * value; });
  return *this;
}

Array &Array::operator*=(const Array &array)
{
  if (!validate_same_shape(*this, array)) return *this;

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array.vector.begin(),
                 this->vector.begin(),
                 [](float v, float a) { return v * a; });
  return *this;
}

Array &Array::operator/=(const float value)
{
  (void)validate_not_zero(value);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 this->vector.begin(),
                 [&value](float v) { return v / value; });
  return *this;
}

Array &Array::operator/=(const Array &array)
{
  if (!validate_same_shape(*this, array)) return *this;

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array.vector.begin(),
                 this->vector.begin(),
                 [](float v, float a) { return v / a; });
  return *this;
}

Array &Array::operator+=(const float value)
{
  std::transform(this->vector.begin(),
                 this->vector.end(),
                 this->vector.begin(),
                 [&value](float v) { return v + value; });
  return *this;
}

Array &Array::operator+=(const Array &array)
{
  if (!validate_same_shape(*this, array)) return *this;

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array.vector.begin(),
                 this->vector.begin(),
                 [](float v, float a) { return v + a; });
  return *this;
}

Array &Array::operator-=(const float value)
{
  std::transform(this->vector.begin(),
                 this->vector.end(),
                 this->vector.begin(),
                 [&value](float v) { return v - value; });
  return *this;
}

Array &Array::operator-=(const Array &array)
{
  if (!validate_same_shape(*this, array)) return *this;

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array.vector.begin(),
                 this->vector.begin(),
                 [](float v, float a) { return v - a; });
  return *this;
}

Array Array::operator*(const float value) const
{
  Array array_out = Array(this->shape);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array_out.vector.begin(),
                 [&value](float v) { return v * value; });
  return array_out;
}

Array Array::operator*(const Array &array) const
{
  if (!validate_same_shape(*this, array)) return Array();

  Array array_out = Array(array.shape);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array.vector.begin(),
                 array_out.vector.begin(),
                 [](float a, float b) { return a * b; });
  return array_out;
}

Array operator*(const float value, const Array &array) // friend function
{
  Array array_out = Array(array.shape);

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [&value](float v) { return v * value; });
  return array_out;
}

Array Array::operator/(const float value) const
{
  (void)validate_not_zero(value);

  Array array_out = Array(this->shape);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array_out.vector.begin(),
                 [&value](float v) { return v / value; });
  return array_out;
}

Array Array::operator/(const Array &array) const
{
  if (!validate_same_shape(*this, array)) return Array();

  Array array_out = Array(array.shape);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array.vector.begin(),
                 array_out.vector.begin(),
                 [](float a, float b) { return a / b; });
  return array_out;
}

Array operator/(const float value, const Array &array) // friend function
{
  Array array_out = Array(array.shape);

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [&value](float v) { return value / v; });
  return array_out;
}

Array Array::operator+(const float value) const
{
  Array array_out = Array(this->shape);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array_out.vector.begin(),
                 [&value](float v) { return v + value; });
  return array_out;
}

Array Array::operator+(const Array &array) const
{
  if (!validate_same_shape(*this, array)) return Array();

  Array array_out = Array(array.shape);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array.vector.begin(),
                 array_out.vector.begin(),
                 [](float a, float b) { return a + b; });
  return array_out;
}

Array operator+(const float value, const Array &array) // friend function
{
  Array array_out = Array(array.shape);

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [&value](float v) { return value + v; });
  return array_out;
}

Array Array::operator-() const
{
  Array array_out = Array(this->shape);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return -v; });
  return array_out;
}

Array Array::operator-(float value) const
{
  Array array_out = Array(this->shape);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array_out.vector.begin(),
                 [&value](float v) { return v - value; });
  return array_out;
}

Array Array::operator-(const Array &array) const
{
  if (!validate_same_shape(*this, array)) return Array();

  Array array_out = Array(array.shape);

  std::transform(this->vector.begin(),
                 this->vector.end(),
                 array.vector.begin(),
                 array_out.vector.begin(),
                 [](float a, float b) { return a - b; });
  return array_out;
}

const Array operator-(const float value, const Array &array) // friend function
{
  Array array_out = Array(array.shape);

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [&value](float v) { return value - v; });
  return array_out;
}

} // namespace hmap
