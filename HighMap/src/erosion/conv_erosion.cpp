/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/convolve.hpp"
#include "highmap/erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/particles.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

std::vector<glm::ivec2> compute_particle_path(const Array  &z,
                                              glm::ivec2    start,
                                              bool          allow_uphill,
                                              float         randomness,
                                              std::uint32_t seed,
                                              float         inertia,
                                              float         exit_weight,
                                              int           border_margin,
                                              Mat<uint8_t> &visited)
{
  const glm::ivec2 &shape = z.shape;

  std::vector<glm::ivec2> path;
  path.reserve(1024);
  path.push_back(start);

  visited(start) = 1;

  glm::ivec2 current = start;
  glm::vec2  prev_dir(0.f);

  std::mt19937                          rng(seed);
  std::uniform_real_distribution<float> rand01(0.f, 1.f);

  static const glm::ivec2 dirs[8] =
      {{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
  static const float dists[8] =
      {1.f, 1.f, 1.f, 1.f, 1.41421356f, 1.41421356f, 1.41421356f, 1.41421356f};
  static const glm::vec2 norm_dirs[8] = {{1.f, 0.f},
                                         {-1.f, 0.f},
                                         {0.f, 1.f},
                                         {0.f, -1.f},
                                         {0.70710678f, 0.70710678f},
                                         {0.70710678f, -0.70710678f},
                                         {-0.70710678f, 0.70710678f},
                                         {-0.70710678f, -0.70710678f}};

  const glm::vec2 half_shape(0.5f * static_cast<float>(shape.x - 1),
                             0.5f * static_cast<float>(shape.y - 1));
  const glm::vec2 inv_half_shape(half_shape.x > 0.f ? 1.f / half_shape.x : 0.f,
                                 half_shape.y > 0.f ? 1.f / half_shape.y : 0.f);

  while (true)
  {
    float current_h = z(current);

    glm::ivec2 best_next = current;
    float      best_score = -std::numeric_limits<float>::infinity();
    bool       found_move = false;

    // direction from the center toward the current particle position
    glm::vec2 center(0.5f * static_cast<float>(shape.x - 1),
                     0.5f * static_cast<float>(shape.y - 1));

    glm::vec2 exit_dir = glm::vec2(current) - center;
    exit_dir /= glm::vec2(shape.x - 1, shape.y - 1);

    float len2 = std::clamp(glm::dot(exit_dir, exit_dir), 0.f, 1.f);

    if (len2 > 0.f)
    {
      exit_dir /= std::sqrt(len2);
    }
    else
    {
      // particle exactly at the center: no outward preference
      exit_dir = glm::vec2(0.f);
    }

    for (int d_idx = 0; d_idx < 8; ++d_idx)
    {
      const auto &d = dirs[d_idx];
      glm::ivec2  n = current + d;

      if (n.x < 0 || n.y < 0 || n.x >= shape.x || n.y >= shape.y) continue;
      if (visited(n)) continue;

      float nh = z(n);
      float dist = dists[d_idx];
      float slope = (current_h - nh) / dist;
      bool  downhill = slope > 0.f;

      if (!downhill && !allow_uphill) continue;

      const glm::vec2 &dir = norm_dirs[d_idx];

      float dir_align = 0.f;
      if (prev_dir != glm::vec2(0.f)) dir_align = glm::dot(prev_dir, dir);

      float score;

      if (downhill)
      {
        // downhill: gravity + inertia
        score = slope + inertia * dir_align;
      }
      else
      {
        // uphill: smallest climb + inertia + preference toward the nearest exit
        float exit_align = glm::dot(dir, exit_dir);

        // negative uphill penalty (scale exit_weight with distance
        // from the center, the closer the particle is from the
        // boundary, the more the particle is attracted to it)
        score = slope + inertia * dir_align + len2 * exit_weight * exit_align;
      }

      if (randomness > 0.f)
      {
        score += randomness * rand01(rng);
      }

      if (!found_move || score > best_score)
      {
        best_score = score;
        best_next = n;
        found_move = true;
      }
    }

    if (!found_move) break;

    glm::ivec2 step = best_next - current;

    prev_dir = glm::normalize(glm::vec2(step));

    current = best_next;

    visited(current) = 1;
    path.push_back(current);

    if (current.x < border_margin || current.y < border_margin ||
        current.x >= shape.x - border_margin ||
        current.y >= shape.y - border_margin)
      break;

    // safe guard, keep the number of iterations reasonable
    if (path.size() >= static_cast<size_t>(shape.x + shape.y)) break;
  }

  for (const auto &p : path)
  {
    visited(p) = 0;
  }

  return path;
}

void conv_erosion(Array        &z,
                  std::uint32_t seed,
                  int           iterations,
                  int           particle_count,
                  int           ir_min,
                  int           ir_max,
                  float         size_distrib_exp,
                  float         erosion_strength,
                  float         randomness,
                  float         exit_forcing,
                  int           gradient_ir,
                  float         gradient_exp,
                  float         gradient_strength_min)
{
  if (!validate_non_empty(z)) return;

  float bulk_amp = 2.f;
  float filling_strength = 0.f;

  // --- Setup

  const glm::ivec2 &shape = z.shape;
  const int         border_margin = 1;

  std::mt19937                          rng(seed);
  std::uniform_real_distribution<float> rand01(0.f, 1.f);

  const float size_scale = static_cast<float>(ir_min) /
                           static_cast<float>(ir_max);

  // erosion kernel
  const glm::ivec2 kernel_shape = {2 * ir_max + 1, 2 * ir_max + 1};
  const Array      kernel = cone(kernel_shape, 2.f);

  // --- Working array (add bulk and depression filling)

  Array z_wrk = z;

  auto prim = hmap::PrimitiveType::PRIM_CONE_SMOOTH;

  if (filling_strength > 0.f)
  {
    depression_filling_priority_flood(z_wrk, /* apply_post_filter */ true);
    z_wrk = lerp(z, z_wrk, filling_strength);
  }

  if (bulk_amp > 0.f) z_wrk = hmap::bulkify(z_wrk, prim, bulk_amp);

  // --- Main loop

  Mat<uint8_t> visited(shape);
  Array        mask(shape);
  Array        size(shape);

  for (int it = 0; it < iterations; ++it)
  {
    mask = 0.f;
    size = 0.f;

    for (int k = 0; k < particle_count; ++k)
    {
      glm::ivec2 spawn = spawn_uniform(rng, shape);

      const std::uint32_t particle_seed = seed ^
                                          (static_cast<uint32_t>(it) * 6271u +
                                           static_cast<uint32_t>(k) * 2053u);

      std::vector<glm::ivec2> path = compute_particle_path(z_wrk,
                                                           spawn,
                                                           true,
                                                           randomness,
                                                           particle_seed,
                                                           /* inertia */ 0.f,
                                                           exit_forcing,
                                                           border_margin,
                                                           visited);

      if (path.empty()) continue;

      const float size_factor = std::pow(rand01(rng), size_distrib_exp);
      const float local_size = lerp(size_scale, 1.f, size_factor);
      const float inv_sum = 1.f / static_cast<float>(path.size());

      for (size_t i = 0; i < path.size(); ++i)
      {
        const glm::ivec2 &p = path[i];

        float t = (i + 1) * inv_sum;

        float v0 = mask(p);
        float v1 = lerp(0.1f, 1.f, t);
        float vsize = local_size * v1;

        mask(p) = std::max(v0, v1);
        size(p) = std::max(size(p), vsize);
      }
    }

    // --- Apply erosion

    Array erosion_field = gpu::sparse_max_convolution(mask, size, kernel, 0.f);

    // gradient scaling
    if (gradient_ir > 0 && gradient_exp != 0.f)
    {
      int   effective_ir = std::min(gradient_ir,
                                  std::max(1, std::min(shape.x, shape.y) / 4));
      Array gn = gradient_norm_filtered(z_wrk, effective_ir);
      if (gn.max() > gn.min())
      {
        remap(gn, gradient_strength_min, 1.f);
        erosion_field *= pow(gn, gradient_exp);
      }
    }

    hmap::laplace(erosion_field, 0.125f, 1);

    z_wrk -= erosion_strength * erosion_field;
  }

  // --- Output

  if (bulk_amp > 0.f) z_wrk = hmap::bulkify(z_wrk, prim, -bulk_amp);

  z = z_wrk;
}

} // namespace hmap::gpu
