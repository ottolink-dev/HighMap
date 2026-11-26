/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "highmap/features/poi.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/range.hpp"

namespace hmap
{

namespace
{

struct Candidate
{
  glm::ivec2 coord;
  glm::vec2  pos;
  float      score;
  float      prominence;
};

} // namespace

Array compute_poi_suitability_field(const Array                &heightmap,
                                    PoiType                     type,
                                    const PoiGenerationOptions &options)
{
  if (!validate_non_empty(heightmap))
  {
    return Array();
  }

  int               ir = std::max(1, options.scale_radius);
  const glm::ivec2 &shape = heightmap.shape;
  Array             suitability(shape);

  Array slope = gradient_norm(heightmap);
  Array slope_norm = slope;
  remap(slope_norm, 0.f, 1.f);

  Array elev_norm = heightmap;
  remap(elev_norm, 0.f, 1.f);

  switch (type)
  {
  case PoiType::APEX:
  {
    Array lmax = local_max(heightmap, ir);
    Array lmin = local_min(heightmap, ir);
    Array relief = lmax - lmin;
    remap(relief, 0.f, 1.f);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float diff_from_max = std::max(0.f, lmax(i, j) - heightmap(i, j));
        float peakness = std::exp(-diff_from_max * 30.f);
        float flatness_peak = 1.f -
                              std::clamp(slope_norm(i, j) * 1.5f, 0.f, 1.f);

        suitability(i, j) = 0.5f * peakness + 0.3f * elev_norm(i, j) +
                            0.1f * relief(i, j) + 0.1f * flatness_peak;
      }
    break;
  }

  case PoiType::SINK:
  {
    Array lmax = local_max(heightmap, ir);
    Array lmin = local_min(heightmap, ir);
    Array relief = lmax - lmin;
    remap(relief, 0.f, 1.f);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float diff_from_min = std::max(0.f, heightmap(i, j) - lmin(i, j));
        float sinkness = std::exp(-diff_from_min * 30.f);

        suitability(i, j) = 0.6f * sinkness + 0.2f * (1.f - elev_norm(i, j)) +
                            0.2f * relief(i, j);
      }
    break;
  }

  case PoiType::SADDLE:
  {
    Array lmax = local_max(heightmap, ir);
    Array lmin = local_min(heightmap, ir);
    Array relief = lmax - lmin;
    remap(relief, 0.f, 1.f);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float rel = relief(i, j);
        float mid_rel = 1.f -
                        std::abs((heightmap(i, j) - lmin(i, j)) /
                                     std::max(1e-4f, lmax(i, j) - lmin(i, j)) -
                                 0.5f) *
                            2.f;
        float low_slope = 1.f - slope_norm(i, j);

        suitability(i, j) = 0.5f * mid_rel + 0.3f * rel + 0.2f * low_slope;
      }
    break;
  }

  case PoiType::RIDGE:
  {
    Array lmax = local_max(heightmap, ir);
    Array lmin = local_min(heightmap, ir);
    Array relief = lmax - lmin;
    remap(relief, 0.f, 1.f);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float rel_elev = (heightmap(i, j) - lmin(i, j)) /
                         std::max(1e-4f, lmax(i, j) - lmin(i, j));
        suitability(i, j) = 0.6f * rel_elev + 0.4f * relief(i, j);
      }
    break;
  }

  case PoiType::VALLEY:
  {
    Array lmax = local_max(heightmap, ir);
    Array lmin = local_min(heightmap, ir);
    Array relief = lmax - lmin;
    remap(relief, 0.f, 1.f);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float rel_depth = 1.f - (heightmap(i, j) - lmin(i, j)) /
                                    std::max(1e-4f, lmax(i, j) - lmin(i, j));
        suitability(i, j) = 0.6f * rel_depth + 0.4f * relief(i, j);
      }
    break;
  }

  case PoiType::CLIFF:
  {
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        suitability(i, j) = slope_norm(i, j);
      }
    break;
  }

  case PoiType::PLATEAU:
  {
    Array lmax = local_max(heightmap, ir);
    Array lmin = local_min(heightmap, ir);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float flatness = 1.f - slope_norm(i, j);
        float rel_elev = (heightmap(i, j) - lmin(i, j)) /
                         std::max(1e-4f, lmax(i, j) - lmin(i, j));

        suitability(i, j) = 0.5f * flatness + 0.3f * elev_norm(i, j) +
                            0.2f * rel_elev;
      }
    break;
  }

  case PoiType::FLAT_LOWLAND:
  {
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float flatness = 1.f - slope_norm(i, j);
        suitability(i, j) = 0.6f * flatness + 0.4f * (1.f - elev_norm(i, j));
      }
    break;
  }

  case PoiType::VIEWPOINT:
  {
    Array lmax = local_max(heightmap, ir);
    Array lmin = local_min(heightmap, ir);
    Array relief = lmax - lmin;
    remap(relief, 0.f, 1.f);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float rel_elev = (heightmap(i, j) - lmin(i, j)) /
                         std::max(1e-4f, lmax(i, j) - lmin(i, j));
        suitability(i, j) = 0.4f * elev_norm(i, j) + 0.4f * rel_elev +
                            0.2f * relief(i, j);
      }
    break;
  }

  case PoiType::SETTLEMENT_CANDIDATE:
  {
    Array lmax = local_max(heightmap, ir);
    Array lmin = local_min(heightmap, ir);
    Array relief = lmax - lmin;
    remap(relief, 0.f, 1.f);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float flatness = 1.f - slope_norm(i, j);
        float mid_elev = 1.f - std::abs(elev_norm(i, j) - 0.4f) * 2.f;
        mid_elev = std::clamp(mid_elev, 0.f, 1.f);

        suitability(i, j) = 0.5f * flatness + 0.3f * mid_elev +
                            0.2f * (1.f - relief(i, j));
      }
    break;
  }

  case PoiType::CUSTOM:
  default:
  {
    suitability = elev_norm;
    break;
  }
  }

  remap(suitability, 0.f, 1.f);
  return suitability;
}

Cloud generate_points_of_interest_cloud(const Array                &heightmap,
                                        PoiType                     type,
                                        const PoiGenerationOptions &options)
{
  auto pois = generate_points_of_interest(heightmap, type, options);

  std::vector<Point> points;
  points.reserve(pois.size());

  for (const auto &poi : pois)
  {
    points.emplace_back(poi.position.x, poi.position.y, poi.score);
  }

  return Cloud(std::move(points));
}

std::vector<PointOfInterest> generate_points_of_interest(
    const Array                &heightmap,
    PoiType                     type,
    const PoiGenerationOptions &options)
{
  if (!validate_non_empty(heightmap))
  {
    return {};
  }

  Array score_field = compute_poi_suitability_field(heightmap, type, options);
  return generate_points_from_score_field(score_field,
                                          heightmap,
                                          options,
                                          type);
}

std::vector<glm::vec2> generate_points_of_interest_positions(
    const Array                &heightmap,
    PoiType                     type,
    const PoiGenerationOptions &options)
{
  auto pois = generate_points_of_interest(heightmap, type, options);

  std::vector<glm::vec2> positions;
  positions.reserve(pois.size());

  for (const auto &poi : pois)
  {
    positions.push_back(poi.position);
  }

  return positions;
}

std::vector<PointOfInterest> generate_points_from_score_field(
    const Array                &score_field,
    const Array                &heightmap,
    const PoiGenerationOptions &options,
    PoiType                     type)
{
  if (!validate_non_empty(score_field) || !validate_non_empty(heightmap))
  {
    return {};
  }

  const glm::ivec2 &shape = score_field.shape;
  int               ir = std::max(1, options.scale_radius);

  Array lmin = local_min(heightmap, ir);
  Array lmax = local_max(heightmap, ir);
  Array slope = gradient_norm(heightmap);
  Array aspect_map = gradient_angle(heightmap);

  float xmin = options.bbox.min.x;
  float xmax = options.bbox.max.x;
  float ymin = options.bbox.min.y;
  float ymax = options.bbox.max.y;
  float dx_domain = xmax - xmin;
  float dy_domain = ymax - ymin;

  int border_px_x = static_cast<int>(
      std::ceil(options.border_margin * shape.x));
  int border_px_y = static_cast<int>(
      std::ceil(options.border_margin * shape.y));

  int start_x = std::max(1, border_px_x);
  int end_x = std::min(shape.x - 1, shape.x - border_px_x);
  int start_y = std::max(1, border_px_y);
  int end_y = std::min(shape.y - 1, shape.y - border_px_y);

  if (start_x >= end_x || start_y >= end_y)
  {
    return {};
  }

  std::vector<Candidate> candidates;

  for (int j = start_y; j < end_y; ++j)
    for (int i = start_x; i < end_x; ++i)
    {
      float s = score_field(i, j);

      if (s < options.threshold)
      {
        continue;
      }

      if (options.p_mask != nullptr)
      {
        if ((*options.p_mask)(i, j) <= 0.5f)
        {
          continue;
        }
      }

      // Check if local 3x3 maximum on score field
      bool is_local_max = true;
      for (int dj = -1; dj <= 1 && is_local_max; ++dj)
      {
        for (int di = -1; di <= 1; ++di)
        {
          if (di == 0 && dj == 0) continue;
          if (score_field(i + di, j + dj) > s)
          {
            is_local_max = false;
            break;
          }
        }
      }

      if (!is_local_max)
      {
        continue;
      }

      float prom = 0.f;
      if (type == PoiType::SINK)
      {
        prom = lmax(i, j) - heightmap(i, j);
      }
      else
      {
        prom = heightmap(i, j) - lmin(i, j);
      }

      if (options.prominence > 0.f && prom < options.prominence)
      {
        continue;
      }

      float     u = (shape.x > 1) ? static_cast<float>(i) / (shape.x - 1) : 0.f;
      float     v = (shape.y > 1) ? static_cast<float>(j) / (shape.y - 1) : 0.f;
      glm::vec2 pos = {xmin + u * dx_domain, ymin + v * dy_domain};

      candidates.push_back({glm::ivec2(i, j), pos, s, prom});
    }

  // Sort candidates by score descending
  std::sort(candidates.begin(),
            candidates.end(),
            [](const Candidate &a, const Candidate &b)
            { return a.score > b.score; });

  std::mt19937                          rng(options.seed);
  std::uniform_real_distribution<float> jitter_dist(-0.5f, 0.5f);

  std::vector<PointOfInterest> selected_pois;
  selected_pois.reserve(std::min(candidates.size(), options.max_points));

  float min_dist_sq = options.min_distance * options.min_distance;

  for (const auto &cand : candidates)
  {
    if (selected_pois.size() >= options.max_points)
    {
      break;
    }

    glm::vec2 candidate_pos = cand.pos;

    if (options.jitter > 0.f)
    {
      float jx = jitter_dist(rng) * options.jitter * options.min_distance;
      float jy = jitter_dist(rng) * options.jitter * options.min_distance;
      candidate_pos.x = std::clamp(candidate_pos.x + jx, xmin, xmax);
      candidate_pos.y = std::clamp(candidate_pos.y + jy, ymin, ymax);
    }

    bool suppressed = false;
    for (const auto &poi : selected_pois)
    {
      glm::vec2 d = candidate_pos - poi.position;
      float     dist_sq = d.x * d.x + d.y * d.y;
      if (dist_sq < min_dist_sq)
      {
        suppressed = true;
        break;
      }
    }

    if (suppressed)
    {
      continue;
    }

    PointOfInterest poi{.position = candidate_pos,
                        .score = cand.score,
                        .elevation = heightmap.get_value_bilinear(
                            candidate_pos.x,
                            candidate_pos.y,
                            UvRect::from_vec4(options.bbox.to_vec4())),
                        .slope = slope.get_value_bilinear(
                            candidate_pos.x,
                            candidate_pos.y,
                            UvRect::from_vec4(options.bbox.to_vec4())),
                        .aspect = aspect_map.get_value_bilinear(
                            candidate_pos.x,
                            candidate_pos.y,
                            UvRect::from_vec4(options.bbox.to_vec4())),
                        .type = type};

    selected_pois.push_back(poi);
  }

  return selected_pois;
}

} // namespace hmap
