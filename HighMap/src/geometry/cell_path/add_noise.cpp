/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "highmap/functions.hpp"
#include "highmap/geometry/cell_path.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

// --- Functions

void add_line_bresenham(CellPath &out, glm::ivec2 a, glm::ivec2 b)
{
  int x0 = a.x, y0 = a.y;
  int x1 = b.x, y1 = b.y;

  int dx = std::abs(x1 - x0);
  int dy = std::abs(y1 - y0);

  int sx = (x0 < x1) ? 1 : -1;
  int sy = (y0 < y1) ? 1 : -1;

  int err = dx - dy;

  while (true)
  {
    out.push_back({x0, y0});

    if (x0 == x1 && y0 == y1) break;

    int e2 = 2 * err;

    if (e2 > -dy)
    {
      err -= dy;
      x0 += sx;
    }

    if (e2 < dx)
    {
      err += dx;
      y0 += sy;
    }
  }
}

void add_noise(CellPath         &path,
               std::uint32_t     seed,
               float             kw,
               float             amp,
               const glm::ivec2 &shape,
               NoiseType         noise_type,
               int               octaves,
               float             weight,
               float             persistence,
               float             lacunarity)
{
  if (!validate_shape(shape)) return;
  if (!validate_min_size(path, 2, "CellPath indices")) return;

  std::unique_ptr<NoiseFunction> p = create_noise_function_from_type(noise_type,
                                                                     {kw, kw},
                                                                     seed);

  hmap::FbmFunction f = hmap::FbmFunction(std::move(p),
                                          octaves,
                                          weight,
                                          persistence,
                                          lacunarity);

  auto noise_fct = f.get_delegate();

  // direction from first to last point
  glm::vec2 p0 = glm::vec2(path.front());
  glm::vec2 p1 = glm::vec2(path.back());

  glm::vec2 dir = p1 - p0;
  float     len = glm::length(dir);
  if (len < 1e-6f) return;

  dir /= len;

  // perpendicular direction (2D rotation)
  glm::vec2 perp(-dir.y, dir.x);

  // prevent zero division / stable sampling scale
  float inv_len = 1.0f / len;

  for (size_t i = 0; i < path.size(); ++i)
  {
    glm::vec2 p = glm::vec2(path[i]);

    // projection along main direction (0 → 1)
    float t = glm::dot(p - p0, dir) * inv_len;

    // noise sampled along path
    float n = noise_fct(t, 0.f, 0.f);

    // apply transverse displacement
    p += perp * (n * amp);

    p.x = std::clamp(int(p.x), 0, shape.x - 1);
    p.y = std::clamp(int(p.y), 0, shape.y - 1);

    path[i] = glm::ivec2(p);
  }

  enforce_path_adjacency(path);
}

void enforce_path_adjacency(CellPath &path)
{
  if (path.size() < 2) return;

  CellPath corrected;
  CellPath segment;

  for (size_t i = 1; i < path.size(); ++i)
  {
    glm::ivec2 prev = corrected.empty() ? path.front() : corrected.back();
    glm::ivec2 curr = path[i];

    if (prev == curr) continue;

    segment.clear();
    add_line_bresenham(segment, prev, curr);

    // skip the first point of segment if corrected already has it
    size_t start_k = corrected.empty() ? 0 : 1;
    for (size_t k = start_k; k < segment.size(); ++k)
    {
      corrected.push_back(segment[k]);
    }
  }

  if (corrected.empty() && !path.empty())
  {
    corrected.push_back(path.front());
  }

  path = std::move(corrected);
}

bool is_path_adjacent(const CellPath &path)
{
  if (path.size() < 2) return true;

  for (size_t i = 1; i < path.size(); ++i)
  {
    const glm::ivec2 delta = path[i] - path[i - 1];
    if (std::abs(delta.x) > 1 || std::abs(delta.y) > 1) return false;
  }

  return true;
}

} // namespace hmap
