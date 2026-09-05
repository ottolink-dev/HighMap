/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <vector>

#include "highmap/authoring.hpp"
#include "highmap/geometry/cell_path.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/core.hpp"

namespace hmap
{

namespace
{

// finite sentinel: the library is built with -ffast-math, so infinities
// must not be relied upon
constexpr float T_UNREACHED = std::numeric_limits<float>::max();

// 8-neighbourhood offsets and the corresponding step lengths
constexpr int   DI[8] = {1, -1, 0, 0, 1, 1, -1, -1};
constexpr int   DJ[8] = {0, 0, 1, -1, 1, -1, 1, -1};
constexpr float STEP[8] = {1.f,
                           1.f,
                           1.f,
                           1.f,
                           1.41421356f,
                           1.41421356f,
                           1.41421356f,
                           1.41421356f};

// Contour vertices converted to (float) pixel coordinates, using the same
// mapping as Cloud::to_array.
std::vector<glm::vec2> contour_to_pixel_coords(const Path &path,
                                               glm::ivec2  shape,
                                               glm::vec4   bbox)
{
  const float ai = (shape.x - 1) / (bbox.y - bbox.x);
  const float aj = (shape.y - 1) / (bbox.w - bbox.z);

  std::vector<glm::vec2> pts;
  pts.reserve(path.points.size());
  for (const Point &p : path.points)
    pts.emplace_back(ai * (p.x - bbox.x), aj * (p.y - bbox.z));
  return pts;
}

// Pixels of the closed polygon outline (Bresenham, closing segment included).
std::vector<glm::ivec2> rasterize_outline(const std::vector<glm::vec2> &pts,
                                          glm::ivec2                    shape)
{
  std::vector<glm::ivec2> cells;
  const size_t            n = pts.size();

  for (size_t k = 0; k < n; ++k)
  {
    const glm::vec2 &a = pts[k];
    const glm::vec2 &b = pts[(k + 1) % n];
    add_line_bresenham(cells,
                       {(int)std::round(a.x), (int)std::round(a.y)},
                       {(int)std::round(b.x), (int)std::round(b.y)});
  }

  // keep in-grid cells only
  std::vector<glm::ivec2> out;
  out.reserve(cells.size());
  for (const glm::ivec2 &c : cells)
    if (c.x >= 0 && c.x < shape.x && c.y >= 0 && c.y < shape.y)
      out.push_back(c);
  return out;
}

// Even-odd scanline fill of the closed polygon, evaluated at pixel centres.
std::vector<char> fill_polygon(const std::vector<glm::vec2> &pts,
                               glm::ivec2                    shape)
{
  std::vector<char>  inside(shape.x * shape.y, 0);
  const size_t       n = pts.size();
  std::vector<float> xs;

  for (int j = 0; j < shape.y; ++j)
  {
    const float y = (float)j;
    xs.clear();

    for (size_t k = 0; k < n; ++k)
    {
      const glm::vec2 &a = pts[k];
      const glm::vec2 &b = pts[(k + 1) % n];
      if ((a.y <= y && b.y > y) || (b.y <= y && a.y > y))
        xs.push_back(a.x + (y - a.y) * (b.x - a.x) / (b.y - a.y));
    }

    std::sort(xs.begin(), xs.end());

    for (size_t m = 0; m + 1 < xs.size(); m += 2)
    {
      const int i0 = std::max(0, (int)std::ceil(xs[m]));
      const int i1 = std::min(shape.x - 1, (int)std::floor(xs[m + 1]));
      for (int i = i0; i <= i1; ++i)
        inside[j * shape.x + i] = 1;
    }
  }
  return inside;
}

// Rasterised representation of the contour set.
struct ContourRaster
{
  glm::ivec2       shape;
  std::vector<int> contour_of; // fixed pixel: contour index, else -1
  std::vector<int> zone;       // free pixel: innermost enclosing contour
                               // (-1 outside all), fixed pixel: -2
  std::vector<int> parent;     // per contour: enclosing contour or -1
};

ContourRaster rasterize_contours(const std::vector<Path> &contours,
                                 glm::ivec2               shape,
                                 glm::vec4                bbox)
{
  const int    n = (int)contours.size();
  const size_t npix = (size_t)shape.x * shape.y;

  ContourRaster r;
  r.shape = shape;
  r.contour_of.assign(npix, -1);
  r.zone.assign(npix, -1);
  r.parent.assign(n, -1);

  std::vector<std::vector<char>> inside(n);
  std::vector<int>               area(n, 0);
  std::vector<int>               rep(n, -1); // representative outline pixel

  for (int k = 0; k < n; ++k)
  {
    std::vector<glm::vec2> pts = contour_to_pixel_coords(contours[k],
                                                         shape,
                                                         bbox);

    for (const glm::ivec2 &c : rasterize_outline(pts, shape))
    {
      const int p = c.y * shape.x + c.x;
      r.contour_of[p] = k;
      if (rep[k] < 0) rep[k] = p;
    }

    inside[k] = fill_polygon(pts, shape);
    area[k] = (int)std::count(inside[k].begin(), inside[k].end(), 1);
  }

  // parent = smallest enclosing contour (by area) containing the
  // representative outline pixel
  for (int k = 0; k < n; ++k)
  {
    if (rep[k] < 0) continue;
    for (int j = 0; j < n; ++j)
    {
      if (j == k || area[j] < area[k] || !inside[j][rep[k]]) continue;
      if (r.parent[k] < 0 || area[j] < area[r.parent[k]]) r.parent[k] = j;
    }
  }

  // zones: paint contours from largest to smallest so that the innermost
  // enclosing contour wins
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(),
            order.end(),
            [&area](int a, int b) { return area[a] > area[b]; });

  for (int k : order)
    for (size_t p = 0; p < npix; ++p)
      if (inside[k][p]) r.zone[p] = k;

  for (size_t p = 0; p < npix; ++p)
    if (r.contour_of[p] >= 0) r.zone[p] = -2;

  return r;
}

// Per-pixel passage time for a front: ((1 - r) + r * E) / rate
std::vector<float> passage_times(const Array *p_probability,
                                 bool         descending,
                                 float        randomness,
                                 std::mt19937 &gen,
                                 size_t        npix)
{
  std::exponential_distribution<float> expo(1.f);
  std::vector<float>                   tau(npix);

  for (size_t p = 0; p < npix; ++p)
  {
    float rate = p_probability ? std::clamp((*p_probability)((int)p), 1e-3f, 0.999f)
                               : 0.5f;
    if (descending) rate = 1.f - rate;

    const float e = randomness > 0.f ? expo(gen) : 1.f;
    tau[p] = ((1.f - randomness) + randomness * e) / rate;
  }
  return tau;
}

struct Front
{
  std::vector<float> t;       // arrival time
  std::vector<int>   label;   // seeding contour index
  std::vector<char>  reached; // whether the front reached the pixel
};

// Dijkstra front propagation from the fixed contour pixels into the free
// pixels of a target zone. For the "near" front a contour seeds its own
// interior zone; for the "far" front it seeds the zone of its parent.
Front propagate_front(const ContourRaster      &r,
                      bool                      to_parent,
                      const std::vector<float> &tau)
{
  const glm::ivec2 shape = r.shape;
  const size_t     npix = r.zone.size();

  Front f;
  f.t.assign(npix, T_UNREACHED);
  f.label.assign(npix, -1);
  f.reached.assign(npix, 0);

  std::vector<int> target(npix, -3); // zone a settled pixel expands into

  using Item = std::pair<float, int>;
  std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

  auto relax = [&](int p, int i, int j, int tgt, int lbl, float t0)
  {
    for (int m = 0; m < 8; ++m)
    {
      const int ii = i + DI[m];
      const int jj = j + DJ[m];
      if (ii < 0 || ii >= shape.x || jj < 0 || jj >= shape.y) continue;

      const int q = jj * shape.x + ii;
      if (r.zone[q] != tgt) continue;

      const float t = t0 + STEP[m] * tau[q];
      if (t < f.t[q])
      {
        f.t[q] = t;
        f.label[q] = lbl;
        f.reached[q] = 1;
        target[q] = tgt;
        pq.emplace(t, q);
      }
    }
    (void)p;
  };

  // seeds
  for (size_t p = 0; p < npix; ++p)
  {
    const int k = r.contour_of[p];
    if (k < 0) continue;
    const int tgt = to_parent ? r.parent[k] : k;
    relax((int)p, (int)(p % shape.x), (int)(p / shape.x), tgt, k, 0.f);
  }

  // propagation
  std::vector<char> settled(npix, 0);
  while (!pq.empty())
  {
    const auto [t, p] = pq.top();
    pq.pop();
    if (settled[p] || t > f.t[p]) continue;
    settled[p] = 1;
    relax(p, p % shape.x, p / shape.x, target[p], f.label[p], t);
  }

  return f;
}

// Priority flood of the progress field T inside one zone, seeded from the
// free pixels adjacent to the zone's own contour. Removes pits by replacing T
// with the running maximum along the flood.
void fill_pits(const ContourRaster &r,
               int                  zone_id,
               const std::vector<char> &eligible,
               std::vector<float>      &T)
{
  const glm::ivec2 shape = r.shape;
  const size_t     npix = r.zone.size();

  using Item = std::pair<float, int>;
  std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
  std::vector<char> settled(npix, 0);

  auto push_neighbours = [&](int p, float key)
  {
    const int i = p % shape.x;
    const int j = p / shape.x;
    for (int m = 0; m < 8; ++m)
    {
      const int ii = i + DI[m];
      const int jj = j + DJ[m];
      if (ii < 0 || ii >= shape.x || jj < 0 || jj >= shape.y) continue;
      const int q = jj * shape.x + ii;
      if (settled[q] || r.zone[q] != zone_id || !eligible[q]) continue;
      pq.emplace(std::max(T[q], key), q);
    }
  };

  // seeds: fixed pixels of the zone's own contour
  for (size_t p = 0; p < npix; ++p)
    if (r.contour_of[p] == zone_id) push_neighbours((int)p, 0.f);

  while (!pq.empty())
  {
    const auto [key, p] = pq.top();
    pq.pop();
    if (settled[p]) continue;
    settled[p] = 1;
    T[p] = key;
    push_neighbours(p, key);
  }
}

} // namespace

Array elevation_from_contours(glm::ivec2                shape,
                              const std::vector<Path>  &contours,
                              const std::vector<float> &elevations,
                              const Array              *p_probability,
                              float                     randomness,
                              uint                      seed,
                              float                     peak_ratio,
                              float                     outside_ratio,
                              glm::vec4                 bbox)
{
  // --- validation
  if (!validate_shape(shape)) return Array();

  if (contours.empty())
  {
    log::error("elevation_from_contours: at least one contour is required");
    return Array();
  }

  if (contours.size() != elevations.size())
  {
    log::error("elevation_from_contours: contours ({}) and elevations ({}) "
               "sizes differ",
               contours.size(),
               elevations.size());
    return Array();
  }

  for (size_t k = 0; k < contours.size(); ++k)
    if (contours[k].points.size() < 3)
    {
      log::error("elevation_from_contours: contour {} has fewer than 3 points",
                 k);
      return Array();
    }

  if (p_probability && !validate_same_shape(shape, *p_probability))
    return Array();

  randomness = std::clamp(randomness, 0.f, 1.f);

  const int    n = (int)contours.size();
  const size_t npix = (size_t)shape.x * shape.y;

  // --- rasterization, zones and nesting tree
  const ContourRaster r = rasterize_contours(contours, shape, bbox);

  // mean elevation gap between nested contours
  float spacing = 0.f;
  int   nspacing = 0;
  for (int k = 0; k < n; ++k)
    if (r.parent[k] >= 0)
    {
      spacing += std::abs(elevations[k] - elevations[r.parent[k]]);
      ++nspacing;
    }
  if (nspacing > 0) spacing /= (float)nspacing;
  if (spacing <= 0.f)
  {
    const auto [emin, emax] = std::minmax_element(elevations.begin(),
                                                  elevations.end());
    spacing = *emax - *emin;
  }
  if (spacing <= 0.f) spacing = 1.f;

  const float elev_min = *std::min_element(elevations.begin(),
                                           elevations.end());

  // --- front propagation
  std::mt19937 gen(seed);

  const std::vector<float> tau_near = passage_times(p_probability,
                                                    false,
                                                    randomness,
                                                    gen,
                                                    npix);
  const std::vector<float> tau_far = passage_times(p_probability,
                                                   true,
                                                   randomness,
                                                   gen,
                                                   npix);

  const Front front_near = propagate_front(r, false, tau_near);
  const Front front_far = propagate_front(r, true, tau_far);

  // per-zone maximum arrival times (zone index shifted by one: -1 -> 0)
  std::vector<float> tnear_max(n + 1, 0.f);
  std::vector<float> tfar_max(n + 1, 0.f);
  for (size_t p = 0; p < npix; ++p)
  {
    const int k = r.zone[p];
    if (k < -1) continue;
    if (front_near.reached[p])
      tnear_max[k + 1] = std::max(tnear_max[k + 1], front_near.t[p]);
    if (front_far.reached[p])
      tfar_max[k + 1] = std::max(tfar_max[k + 1], front_far.t[p]);
  }

  // --- progress field
  std::vector<float> T(npix, 0.f);
  std::vector<char>  both(npix, 0); // reached by both fronts

  for (size_t p = 0; p < npix; ++p)
  {
    const int k = r.zone[p];
    if (k < -1) continue;

    const bool  ra = front_near.reached[p];
    const bool  rb = front_far.reached[p];
    const float ta = front_near.t[p];
    const float tb = front_far.t[p];

    if (ra && rb)
    {
      T[p] = ta / (ta + tb);
      both[p] = 1;
    }
    else if (ra)
      T[p] = tnear_max[k + 1] > 0.f ? ta / tnear_max[k + 1] : 0.f;
    else if (rb)
      T[p] = tfar_max[k + 1] > 0.f ? tb / tfar_max[k + 1] : 0.f;
  }

  // pit removal in rising zones (all children higher than the contour)
  for (int k = 0; k < n; ++k)
  {
    bool rising = false;
    bool has_child = false;
    for (int c = 0; c < n; ++c)
      if (r.parent[c] == k)
      {
        has_child = true;
        rising = elevations[c] > elevations[k];
        if (!rising) break;
      }
    if (has_child && rising) fill_pits(r, k, both, T);
  }

  // --- elevation
  Array z(shape, 0.f);

  for (size_t p = 0; p < npix; ++p)
  {
    const int k = r.zone[p];

    if (k == -2) // fixed contour pixel
    {
      z((int)p) = elevations[r.contour_of[p]];
      continue;
    }

    const bool ra = front_near.reached[p];
    const bool rb = front_far.reached[p];
    float      h_near, h_far;

    if (ra && rb)
    {
      h_near = elevations[k];
      h_far = elevations[front_far.label[p]];
    }
    else if (ra) // leaf zone: peak or basin
    {
      const bool up = r.parent[k] < 0 ||
                      elevations[k] >= elevations[r.parent[k]];
      h_near = elevations[k];
      h_far = elevations[k] + (up ? 1.f : -1.f) * peak_ratio * spacing;
    }
    else if (rb) // outside zone (or zone without own front)
    {
      h_near = elevations[front_far.label[p]];
      h_far = k < 0 ? h_near - outside_ratio * spacing : elevations[k];
    }
    else // unreachable pixel
    {
      z((int)p) = k >= 0 ? elevations[k] : elev_min - outside_ratio * spacing;
      continue;
    }

    z((int)p) = lerp(h_near, h_far, T[p]);
  }

  return z;
}

} // namespace hmap
