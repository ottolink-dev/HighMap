/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <utility>
#include <vector>

#include "highmap/geometry/cell_path.hpp"

namespace hmap
{

CellPath::CellPath(const std::vector<glm::ivec2> &indices) : indices(indices)
{
}

CellPath::CellPath(std::vector<glm::ivec2> &&indices) noexcept
    : indices(std::move(indices))
{
}

glm::ivec2 &CellPath::back()
{
  return this->indices.back();
}

const glm::ivec2 &CellPath::back() const
{
  return this->indices.back();
}

std::vector<glm::ivec2>::iterator CellPath::begin() noexcept
{
  return this->indices.begin();
}

std::vector<glm::ivec2>::const_iterator CellPath::begin() const noexcept
{
  return this->indices.begin();
}

void CellPath::clear() noexcept
{
  this->indices.clear();
}

bool CellPath::empty() const noexcept
{
  return this->indices.empty();
}

std::vector<glm::ivec2>::iterator CellPath::end() noexcept
{
  return this->indices.end();
}

std::vector<glm::ivec2>::const_iterator CellPath::end() const noexcept
{
  return this->indices.end();
}

glm::ivec2 &CellPath::front()
{
  return this->indices.front();
}

const glm::ivec2 &CellPath::front() const
{
  return this->indices.front();
}

std::vector<glm::ivec2> &CellPath::get_indices()
{
  return this->indices;
}

const std::vector<glm::ivec2> &CellPath::get_indices() const
{
  return this->indices;
}

glm::ivec2 &CellPath::operator[](size_t index)
{
  return this->indices[index];
}

const glm::ivec2 &CellPath::operator[](size_t index) const
{
  return this->indices[index];
}

void CellPath::push_back(const glm::ivec2 &idx)
{
  this->indices.push_back(idx);
}

size_t CellPath::size() const noexcept
{
  return this->indices.size();
}

} // namespace hmap
