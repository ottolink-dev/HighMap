/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/primitives/geo.hpp"

namespace hmap
{

Array island_chain_land_mask(glm::ivec2    shape,
                             const Path   &path,
                             std::uint32_t seed,
                             int           island_count,
                             float         island_radius,
                             float         size_falloff,
                             float         size_jitter,
                             float         scatter,
                             float         displacement,
                             NoiseType     noise_type,
                             float         kw,
                             int           octaves,
                             float         weight,
                             float         persistence,
                             float         lacunarity,
                             glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_non_empty(path, "Path")) return Array(shape);

  Array mask(shape);

  if (island_count < 1) return mask;

  const std::vector<float> arc = path.get_arc_length();

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(-1.f, 1.f);

  for (int k = 0; k < island_count; k++)
  {
    // deterministic draw order: always 3 draws per island
    float jt = dis(gen);
    float jn = dis(gen);
    float jr = dis(gen);

    float t = island_count == 1 ? 0.5f : float(k) / float(island_count - 1);
    t += scatter * jt;

    glm::vec2 pos, tg;
    glm::vec3 point = path.sample_at(t, &arc, &tg);
    pos = {point.x, point.y};

    glm::vec2 normal = {-tg.y, tg.x};
    pos += scatter * jn * normal;

    // signed size falloff: +1 shrinks toward the path end, -1 toward start
    float fall = size_falloff >= 0.f ? 1.f - size_falloff * t
                                     : 1.f + size_falloff * (1.f - t);
    float radius_factor = std::max(0.f, fall) * (1.f + size_jitter * jr);
    float rk = island_radius * radius_factor;

    if (rk <= 0.f) continue;

    with_noise_function(
        noise_type,
        glm::vec2(0.5f * kw, 0.5f * kw),
        seed + std::uint32_t(k),
        [&](const auto &base_fn)
        {
          FbmFunction f(base_fn, octaves, weight, persistence, lacunarity);

          // pixel window around the island (same x/y mapping as
          // island_land_mask)
          float lx = bbox.y - bbox.x;
          float ly = bbox.w - bbox.z;
          float margin = rk + displacement;

          int i1 = std::clamp(
              int((pos.x - margin - bbox.x) / lx * float(shape.x - 1)),
              0,
              shape.x - 1);
          int i2 = std::clamp(
              int((pos.x + margin - bbox.x) / lx * float(shape.x - 1)) + 1,
              0,
              shape.x - 1);
          int j1 = std::clamp(
              int((pos.y - margin - bbox.z) / ly * float(shape.y - 1)),
              0,
              shape.y - 1);
          int j2 = std::clamp(
              int((pos.y + margin - bbox.z) / ly * float(shape.y - 1)) + 1,
              0,
              shape.y - 1);

          for (int j = j1; j <= j2; j++)
            for (int i = i1; i <= i2; i++)
            {
              float x = lx * float(i) / float(shape.x - 1) + bbox.x - pos.x;
              float y = ly * float(j) / float(shape.y - 1) + bbox.z - pos.y;

              float r = std::hypot(x, y);
              float theta = std::atan2(y, x);

              float dr = radius_factor * displacement *
                         f(std::cos(theta), std::sin(theta), 0.f);

              if (r < rk + dr) mask(i, j) = 1.f;
            }
        });
  }

  return mask;
}

} // namespace hmap
