/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

void find_flow_apex(const Array &z, std::vector<int> &is, std::vector<int> &js)
{
  is.clear();
  js.clear();
  if (!validate_non_empty(z)) return;

  Array d8 = flow_direction_d8(z);
  Array nidp = d8_compute_ndip(d8);

  const int nx = z.shape.x;
  const int ny = z.shape.y;

  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      if (nidp(i, j) == 0)
      {
        is.push_back(i);
        js.push_back(j);
      }
    }
}

void find_flow_sinks(const Array &z, std::vector<int> &is, std::vector<int> &js)
{
  is.clear();
  js.clear();
  if (!validate_non_empty(z)) return;

  const auto         &neighbors = neighborhood::MOORE_8;
  const std::uint32_t nb = neighbors.size();

  for (int j = 1; j < z.shape.y - 1; j++)
    for (int i = 1; i < z.shape.x - 1; i++)
    {
      // count neighbor cells with a higher elevation than the
      // current cell
      int n_higher_cells = 0;
      for (size_t k = 0; k < nb; k++)
      {
        int ik = i + neighbors[k].offset.x;
        int jk = j + neighbors[k].offset.y;
        if (z(i, j) < z(ik, jk)) n_higher_cells++;
      }

      if (n_higher_cells == 8)
      {
        is.push_back(i);
        js.push_back(j);
      }
    }
}

std::vector<glm::ivec2> find_flow_sinks(const Array &z)
{
  if (!validate_non_empty(z)) return {};

  std::vector<glm::ivec2> indices;

  const auto         &neighbors = neighborhood::MOORE_8;
  const std::uint32_t nb = neighbors.size();

  for (int j = 1; j < z.shape.y - 1; j++)
    for (int i = 1; i < z.shape.x - 1; i++)
    {
      // count neighbor cells with a higher elevation than the
      // current cell
      int n_higher_cells = 0;
      for (size_t k = 0; k < nb; k++)
      {
        int ik = i + neighbors[k].offset.x;
        int jk = j + neighbors[k].offset.y;
        if (z(i, j) < z(ik, jk)) n_higher_cells++;
      }

      if (n_higher_cells == 8) indices.push_back({i, j});
    }

  return indices;
}

std::vector<glm::ivec2> find_flow_sinks_border(const Array &z)
{
  if (!validate_non_empty(z)) return {};

  std::vector<glm::ivec2> indices;

  const int rows = z.shape.x;
  const int cols = z.shape.y;

  const auto &neighbors = neighborhood::MOORE_8;

  auto is_inside = [&](int i, int j)
  { return i >= 0 && j >= 0 && i < rows && j < cols; };

  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
    {
      // only border cells
      if (!(i == 0 || j == 0 || i == rows - 1 || j == cols - 1)) continue;

      int valid_neighbors = 0;
      int higher_neighbors = 0;

      for (const auto &nbr : neighbors)
      {
        int ni = i + nbr.offset.x;
        int nj = j + nbr.offset.y;

        if (!is_inside(ni, nj)) continue;

        valid_neighbors++;

        if (z(i, j) < z(ni, nj)) higher_neighbors++;
      }

      // sink: all valid neighbors are higher
      if (valid_neighbors > 0 && higher_neighbors == valid_neighbors)
        indices.push_back(glm::ivec2(i, j));
    }

  return indices;
}

} // namespace hmap
