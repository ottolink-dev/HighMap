/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>

#include <glm/geometric.hpp>

#include "highmap/math/core.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/primitives/geo.hpp"

namespace hmap
{

namespace
{

struct NearestPlates
{
  int   a = 0;
  int   b = 0;
  int   c = 0;
  float da = 0.0f;
  float db = 0.0f;
  float dc = 0.0f;
};

float cross2d(glm::vec2 a, glm::vec2 b)
{
  return a.x * b.y - a.y * b.x;
}

NearestPlates find_nearest_plates(const std::vector<TectonicPlate> &plates,
                                  glm::vec2                         p)
{
  NearestPlates res;
  res.da = std::numeric_limits<float>::max();
  res.db = std::numeric_limits<float>::max();
  res.dc = std::numeric_limits<float>::max();

  for (int i = 0; i < static_cast<int>(plates.size()); ++i)
  {
    glm::vec2 diff = p - plates[i].center;
    float     dist = glm::length(diff);

    if (dist < res.da)
    {
      res.dc = res.db;
      res.c = res.b;
      res.db = res.da;
      res.b = res.a;
      res.da = dist;
      res.a = i;
    }
    else if (dist < res.db)
    {
      res.dc = res.db;
      res.c = res.b;
      res.db = dist;
      res.b = i;
    }
    else if (dist < res.dc)
    {
      res.dc = dist;
      res.c = i;
    }
  }

  return res;
}

void evaluate_pair_features(int                               idx0,
                            int                               idx1,
                            float                             d0,
                            float                             d1,
                            const std::vector<TectonicPlate> &plates,
                            float                             width,
                            float                             relief_amp,
                            float                             subduction_offset,
                            float                             fault_amp,
                            float                             craggy,
                            float                            &out_feature_h,
                            float                            &out_range,
                            float                            &out_fault)
{
  bool cont0 = plates[idx0].is_continental;
  bool cont1 = plates[idx1].is_continental;

  if (!cont0 && !cont1) return;

  glm::vec2 c0 = plates[idx0].center;
  glm::vec2 c1 = plates[idx1].center;
  glm::vec2 diff = c1 - c0;
  float     len = glm::length(diff);
  glm::vec2 n = (len > 1e-6f) ? (diff / len) : glm::vec2(1.0f, 0.0f);

  float d = std::abs(d1 - d0) * 0.5f;
  float d_norm = d / width;
  float near_falloff = std::exp(-d_norm * d_norm * 4.0f);
  float wide_falloff = std::exp(-d_norm * d_norm * 0.8f);

  glm::vec2 rel_v = plates[idx0].velocity - plates[idx1].velocity;
  float     closing = glm::dot(rel_v, n);
  float     sliding = std::abs(cross2d(rel_v, n));

  if (closing > 0.0f)
  {
    float conv_strength = std::min(1.0f, closing);

    if (cont0 && cont1)
    {
      float r = near_falloff * craggy * conv_strength;
      out_range += r;
      out_feature_h += r * relief_amp;
    }
    else
    {
      float d_cont = cont0 ? d : -d;
      float range_d = (d_cont - subduction_offset * width) / width;
      float trench_d = (d_cont + 0.35f * width) / width;

      float inland_falloff = std::exp(-range_d * range_d * 4.0f);
      float trench_falloff = std::exp(-trench_d * trench_d * 4.0f);

      float r = inland_falloff * craggy * conv_strength * 0.85f;
      out_range += r;
      out_feature_h += r * relief_amp;
      out_feature_h -= trench_falloff * 0.3f * relief_amp * conv_strength;
    }
  }
  else
  {
    float div_strength = std::min(1.0f, -closing);

    if (cont0 && cont1)
    {
      out_feature_h += (0.12f * wide_falloff - 0.32f * near_falloff) *
                       relief_amp * div_strength;
    }
    else
    {
      out_feature_h += near_falloff * 0.2f * relief_amp * div_strength;
    }
  }

  float fault_val = near_falloff * std::cos(d_norm * 3.0f) *
                    std::min(1.0f, sliding);
  out_feature_h += fault_val * fault_amp;
  out_fault += std::abs(fault_val);
}

// --- Plate scattering

std::vector<TectonicPlate> scatter_tectonic_plates(int   plate_count,
                                                   uint  seed,
                                                   float continental_ratio)
{
  int count = std::max(1, plate_count);
  int g = static_cast<int>(std::ceil(std::sqrt(static_cast<float>(count))));
  int total_cells = g * g;

  std::vector<int> cell_indices(total_cells);
  std::iota(cell_indices.begin(), cell_indices.end(), 0);

  std::mt19937 rng(seed);
  std::shuffle(cell_indices.begin(), cell_indices.end(), rng);

  std::uniform_real_distribution<float> dist_01(0.0f, 1.0f);
  std::uniform_real_distribution<float> dist_jitter(0.15f, 0.85f);
  std::uniform_real_distribution<float> dist_speed(0.4f, 1.0f);
  std::uniform_real_distribution<float> dist_angle(0.0f, 6.28318530718f);

  std::vector<TectonicPlate> plates;
  plates.reserve(count);

  float cell_size = 1.0f / static_cast<float>(g);

  for (int i = 0; i < count; ++i)
  {
    int cell = cell_indices[i];
    int cx = cell % g;
    int cy = cell / g;

    TectonicPlate plate;
    plate.center.x = (static_cast<float>(cx) + dist_jitter(rng)) * cell_size;
    plate.center.y = (static_cast<float>(cy) + dist_jitter(rng)) * cell_size;

    float angle = dist_angle(rng);
    float speed = dist_speed(rng);
    plate.velocity = {std::cos(angle) * speed, std::sin(angle) * speed};

    plate.is_continental = (i == 0) ? true : (dist_01(rng) < continental_ratio);

    plates.push_back(plate);
  }

  return plates;
}

} // namespace

// --- Per-sample evaluation

Array tectonic_relief(glm::ivec2                        shape,
                      const std::vector<TectonicPlate> &plates,
                      float                             roughness,
                      float                             boundary_width,
                      float                             relief_amp,
                      float                             subduction_offset,
                      float                             fault_amp,
                      const Array                      *p_dx,
                      const Array                      *p_dy,
                      const Array                      *p_coast_noise,
                      const Array                      *p_ridge_noise,
                      Array                            *p_ranges,
                      Array                            *p_faults,
                      const UvRect                     &uv_rect)
{
  Array height(shape);

  if (p_ranges) *p_ranges = Array(shape);
  if (p_faults) *p_faults = Array(shape);

  if (plates.empty()) return height;

  float width = std::max(1e-4f, boundary_width);
  float du = (uv_rect.u_max - uv_rect.u_min) / static_cast<float>(shape.x);
  float dv = (uv_rect.v_max - uv_rect.v_min) / static_cast<float>(shape.y);

  for (int j = 0; j < shape.y; ++j)
  {
    float v = uv_rect.v_min + (static_cast<float>(j) + 0.5f) * dv;

    for (int i = 0; i < shape.x; ++i)
    {
      float u = uv_rect.u_min + (static_cast<float>(i) + 0.5f) * du;

      float     px = u + (p_dx ? (*p_dx)(i, j) : 0.0f);
      float     pz = v + (p_dy ? (*p_dy)(i, j) : 0.0f);
      glm::vec2 p(px, pz);

      // 1. nearest plates
      NearestPlates np = find_nearest_plates(plates, p);

      if (plates.size() == 1)
      {
        float cont = plates[0].is_continental ? 1.0f : 0.0f;
        height(i, j) = cont * 0.5f - (1.0f - cont) * 0.4f;
        continue;
      }

      // 2. continuous continental mask from nearest continental vs oceanic
      // plate
      float c_val = p_coast_noise ? (*p_coast_noise)(i, j) : 0.0f;
      float min_cont_dist = std::numeric_limits<float>::max();
      float min_ocean_dist = std::numeric_limits<float>::max();
      bool  has_cont = false;
      bool  has_ocean = false;

      for (int k = 0; k < static_cast<int>(plates.size()); ++k)
      {
        float dk = glm::length(p - plates[k].center);
        if (plates[k].is_continental)
        {
          min_cont_dist = std::min(min_cont_dist, dk);
          has_cont = true;
        }
        else
        {
          min_ocean_dist = std::min(min_ocean_dist, dk);
          has_ocean = true;
        }
      }

      float cont = 0.0f;
      if (has_cont && !has_ocean)
      {
        cont = 1.0f;
      }
      else if (has_cont && has_ocean)
      {
        float margin = (min_ocean_dist - min_cont_dist) * 0.5f;
        float perturbed = margin + c_val * roughness * width * 0.6f;
        cont = smoothstep3(perturbed, -width * 0.5f, width * 0.5f);
      }

      // base elevation: continents above 0, oceans below 0
      float base_h = cont * 0.35f - (1.0f - cont) * 0.45f;
      base_h += c_val * roughness * 0.08f;

      // 3. boundary features for nearest pair (a, b), faded at triple junctions
      float r_val = p_ridge_noise ? (*p_ridge_noise)(i, j) : 0.0f;
      float craggy = 0.45f + 0.55f * (1.0f - std::abs(r_val));
      float range_val = 0.0f;
      float fault_val = 0.0f;
      float feature_h = 0.0f;

      evaluate_pair_features(np.a,
                             np.b,
                             np.da,
                             np.db,
                             plates,
                             width,
                             relief_amp,
                             subduction_offset,
                             fault_amp,
                             craggy,
                             feature_h,
                             range_val,
                             fault_val);

      if (plates.size() > 2)
      {
        float tj_fade = smoothstep3((np.dc - np.db) * 0.5f, 0.0f, width * 0.5f);
        feature_h *= tj_fade;
        range_val *= tj_fade;
        fault_val *= tj_fade;
      }

      height(i, j) = base_h + feature_h;

      if (p_ranges) (*p_ranges)(i, j) = range_val;
      if (p_faults) (*p_faults)(i, j) = fault_val;
    }
  }

  return height;
}

Array tectonic_relief(glm::ivec2    shape,
                      int           plate_count,
                      uint          seed,
                      float         continental_ratio,
                      float         roughness,
                      float         boundary_width,
                      float         relief_amp,
                      float         subduction_offset,
                      float         fault_amp,
                      const Array  *p_dx,
                      const Array  *p_dy,
                      const Array  *p_coast_noise,
                      const Array  *p_ridge_noise,
                      Array        *p_ranges,
                      Array        *p_faults,
                      const UvRect &uv_rect)
{
  std::vector<TectonicPlate> plates = scatter_tectonic_plates(
      plate_count,
      seed,
      continental_ratio);

  return tectonic_relief(shape,
                         plates,
                         roughness,
                         boundary_width,
                         relief_amp,
                         subduction_offset,
                         fault_amp,
                         p_dx,
                         p_dy,
                         p_coast_noise,
                         p_ridge_noise,
                         p_ranges,
                         p_faults,
                         uv_rect);
}

} // namespace hmap
