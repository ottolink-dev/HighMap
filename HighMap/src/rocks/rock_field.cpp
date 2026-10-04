/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>

#include "highmap/functions.hpp"
#include "highmap/geometry/inverse_sampler_2d.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/rocks/rock_field.hpp"
#include "highmap/rocks/rock_simulation.hpp"

#include <unordered_map>

namespace hmap
{

// ==========================================================================
//  Constructors
// ==========================================================================

RockField::RockField(const std::vector<Rock> &rocks) : ScatterField(rocks)
{
}

RockField::RockField(std::vector<Rock> &&rocks) noexcept
    : ScatterField(std::move(rocks))
{
}

RockField::RockField(const ScatterField &field) : ScatterField(field)
{
}

RockField::RockField(ScatterField &&field) noexcept
    : ScatterField(std::move(field))
{
}

RockField::RockField(const Cloud &cloud,
                     uint32_t     class_id,
                     float        default_radius)
    : ScatterField(cloud, class_id, default_radius)
{
}

// ==========================================================================
//  Rock-Specific Operations
// ==========================================================================

void RockField::apply_power_law_distribution(float    alpha,
                                             float    r_min,
                                             float    r_max,
                                             uint32_t seed)
{
  if (items.empty()) return;

  float r0 = std::min(r_min, r_max);
  float r1 = std::max(r_min, r_max);
  if (r0 <= 0.0f) r0 = 1e-4f;
  if (r1 <= r0) r1 = r0 * 1.01f;

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(0.0f, 1.0f);

  float beta = alpha - 1.0f;

  for (auto &rock : items)
  {
    float u = dis(gen);
    float radius = r0;

    if (std::abs(beta) < 1e-4f)
    {
      // log-uniform fallback when alpha == 1
      radius = r0 * std::exp(u * std::log(r1 / r0));
    }
    else
    {
      float inv_r0 = std::pow(r0, -beta);
      float inv_r1 = std::pow(r1, -beta);
      float val = inv_r0 - u * (inv_r0 - inv_r1);
      radius = std::pow(val, -1.0f / beta);
    }

    rock.radius = std::clamp(radius, r0, r1);
  }
}

void RockField::apply_power_law_distribution(const RockDistribution &dist,
                                             uint32_t                seed)
{
  apply_power_law_distribution(dist.power_law_alpha,
                               dist.radius_min,
                               dist.radius_max,
                               seed);
}

void RockField::apply_slope_sorting(const Array     &slope,
                                    float            sorting_strength,
                                    const glm::vec4 &bbox)
{
  if (items.size() < 2 || !validate_non_empty(slope)) return;

  float strength = std::clamp(sorting_strength, 0.0f, 1.0f);
  if (strength <= 1e-6f) return;

  // --- Sample Local Slopes at Rock Positions

  auto               slope_fct = make_xy_function_from_array(slope, bbox);
  std::vector<float> sampled_slopes(items.size());
  for (size_t i = 0; i < items.size(); ++i)
  {
    sampled_slopes[i] = slope_fct(items[i].position.x, items[i].position.y);
  }

  // --- Sort Rocks by Local Slope Ascending (Gentle to Steep)

  std::vector<size_t> slope_order(items.size());
  std::iota(slope_order.begin(), slope_order.end(), 0);
  std::stable_sort(slope_order.begin(),
                   slope_order.end(),
                   [&](size_t a, size_t b)
                   { return sampled_slopes[a] < sampled_slopes[b]; });

  // --- Collect and Sort Existing Radii Descending (Largest to Smallest)

  std::vector<float> sorted_radii;
  sorted_radii.reserve(items.size());
  for (const auto &rock : items)
  {
    sorted_radii.push_back(rock.radius);
  }
  std::sort(sorted_radii.begin(), sorted_radii.end(), std::greater<float>());

  // --- Blend Radii Towards Gravitational Slope Sorting

  for (size_t rank = 0; rank < items.size(); ++rank)
  {
    size_t rock_idx = slope_order[rank];
    float  target_radius = sorted_radii[rank];
    items[rock_idx].radius = (1.0f - strength) * items[rock_idx].radius +
                             strength * target_radius;
  }
}

void RockField::pack_interstitial_rocks(size_t           count,
                                        const Array     &density,
                                        float            r_min,
                                        float            r_max,
                                        uint32_t         class_id,
                                        uint32_t         seed,
                                        const glm::vec4 &bbox)
{
  if (count == 0) return;

  glm::vec4 effective_bbox = bbox;
  if (std::abs(effective_bbox.y - effective_bbox.x) < 1e-6f)
  {
    effective_bbox = get_bbox();
  }

  float width = effective_bbox.y - effective_bbox.x;
  float height = effective_bbox.w - effective_bbox.z;
  if (width <= 1e-6f || height <= 1e-6f) return;

  // --- Build Spatial Grid from Existing Rocks

  float max_r = r_max;
  for (const auto &item : items)
  {
    max_r = std::max(max_r, item.radius);
  }

  float                                            cell_size = 2.0f * max_r;
  std::unordered_map<int64_t, std::vector<size_t>> grid;

  auto compute_cell = [&](float x, float y) -> std::pair<int, int>
  {
    if (cell_size <= 1e-7f) return {0, 0};
    int gx = static_cast<int>(std::floor(x / cell_size));
    int gy = static_cast<int>(std::floor(y / cell_size));
    return {gx, gy};
  };

  auto make_key = [](int gx, int gy) -> int64_t
  {
    return (static_cast<int64_t>(gx) << 32) ^
           (static_cast<int64_t>(gy) & 0xFFFFFFFF);
  };

  for (size_t i = 0; i < items.size(); ++i)
  {
    auto [gx, gy] = compute_cell(items[i].position.x, items[i].position.y);
    grid[make_key(gx, gy)].push_back(i);
  }

  // --- Sample Interstitial Rock Candidates

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis_x(effective_bbox.x,
                                              effective_bbox.y);
  std::uniform_real_distribution<float> dis_y(effective_bbox.z,
                                              effective_bbox.w);
  std::uniform_real_distribution<float> dis_r(r_min, r_max);

  bool has_density = !density.vector.empty() && density.shape.x > 0 &&
                     density.shape.y > 0;
  std::unique_ptr<InverseSampler2D> sampler;
  if (has_density)
  {
    sampler = std::make_unique<InverseSampler2D>(density, seed, effective_bbox);
    if (sampler->get_total_weight() <= 1e-7f)
    {
      has_density = false;
      sampler.reset();
    }
  }

  size_t max_attempts = count * 20;
  size_t accepted = 0;

  for (size_t attempt = 0; attempt < max_attempts && accepted < count;
       ++attempt)
  {
    float cx = 0.0f;
    float cy = 0.0f;
    if (has_density)
    {
      auto pt = sampler->sample();
      cx = pt.x;
      cy = pt.y;
    }
    else
    {
      cx = dis_x(gen);
      cy = dis_y(gen);
    }
    float cr = dis_r(gen);

    auto [gx, gy] = compute_cell(cx, cy);
    bool collides = false;

    for (int dy = -1; dy <= 1 && !collides; ++dy)
    {
      for (int dx = -1; dx <= 1 && !collides; ++dx)
      {
        int64_t key = make_key(gx + dx, gy + dy);
        auto    it = grid.find(key);
        if (it == grid.end()) continue;

        for (size_t other_idx : it->second)
        {
          const auto &other = items[other_idx];
          float       dist_x = cx - other.position.x;
          float       dist_y = cy - other.position.y;
          float       dist_sq = dist_x * dist_x + dist_y * dist_y;

          float min_dist = cr + other.radius;
          if (dist_sq < min_dist * min_dist)
          {
            collides = true;
            break;
          }
        }
      }
    }

    if (!collides)
    {
      items.emplace_back(cx, cy, 0.0f, class_id, cr);
      grid[make_key(gx, gy)].push_back(items.size() - 1);
      accepted++;
    }
  }
}

void RockField::simulate_physics(const Array                 &elevation,
                                 const Array                 *p_friction_map,
                                 const RockSimulationOptions &options)
{
  *this = simulate_rock_trajectories(*this, elevation, p_friction_map, options);
}

std::string RockField::to_string() const
{
  std::ostringstream ss;
  ss.imbue(std::locale("C"));

  glm::vec4 bbox = get_bbox();

  std::map<uint32_t, size_t> rock_counts;
  std::map<uint32_t, float>  rock_radii_sum;
  for (const auto &rock : items)
  {
    rock_counts[rock.class_id]++;
    rock_radii_sum[rock.class_id] += rock.radius;
  }

  ss << "RockField Infos\n";
  ss << "--------------------------------\n";
  ss << " total rocks : " << items.size() << "\n";
  ss << " rock classes count : " << rock_counts.size() << "\n";
  ss << " bbox : {" << bbox.x << ", " << bbox.y << ", " << bbox.z << ", "
     << bbox.w << "}\n";

  if (!rock_counts.empty())
  {
    ss << " rock classes breakdown :\n";
    for (const auto &[r_id, count] : rock_counts)
    {
      float avg_r = (count > 0) ? (rock_radii_sum[r_id] / float(count)) : 0.f;
      float pct = (items.empty())
                      ? 0.f
                      : (100.f * float(count) / float(items.size()));
      ss << "   - rock class " << r_id << ": " << count << " rocks ("
         << std::fixed << std::setprecision(1) << pct << "%, avg radius "
         << std::setprecision(4) << avg_r << ")\n";
    }
  }
  ss << "--------------------------------";

  return ss.str();
}

} // namespace hmap
