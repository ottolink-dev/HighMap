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
constexpr float STEP[8] =
    {1.f, 1.f, 1.f, 1.f, 1.41421356f, 1.41421356f, 1.41421356f, 1.41421356f};

// Contour vertices converted to (float) pixel coordinates, using the same
// mapping as Cloud::to_array.
std::vector<glm::vec2> contour_to_pixel_coords(const Path &path,
                                               glm::ivec2  shape,
                                               glm::vec4   bbox)
{
  const float ai = (shape.x - 1) / (bbox.y - bbox.x);
  const float aj = (shape.y - 1) / (bbox.w - bbox.z);

  std::vector<glm::vec2> pts;
  pts.reserve(path.size());
  for (const Point &p : path)
    pts.emplace_back(ai * (p.x - bbox.x), aj * (p.y - bbox.z));
  return pts;
}

// Pixels of the polygon outline (Bresenham, closing segment included only if
// closed).
std::vector<glm::ivec2> rasterize_outline(const std::vector<glm::vec2> &pts,
                                          glm::ivec2                    shape,
                                          bool is_closed = true)
{
  CellPath     cells;
  const size_t n = pts.size();
  if (n == 0) return cells.get_indices();

  const size_t num_segments = is_closed ? n : n - 1;
  for (size_t k = 0; k < num_segments; ++k)
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

// Even-odd scanline fill of the closed polygon, evaluated at pixel centres (j +
// 0.5).
std::vector<char> fill_polygon(const std::vector<glm::vec2> &pts,
                               glm::ivec2                    shape)
{
  std::vector<char>  inside(shape.x * shape.y, 0);
  const size_t       n = pts.size();
  std::vector<float> xs;

  for (int j = 0; j < shape.y; ++j)
  {
    const float y = (float)j + 0.5f;
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
      const int i0 = std::max(0, (int)std::ceil(xs[m] - 0.5f));
      const int i1 = std::min(shape.x - 1, (int)std::floor(xs[m + 1] + 0.5f));
      for (int i = i0; i <= i1; ++i)
        inside[j * shape.x + i] = 1;
    }
  }
  return inside;
}

// Classification of a point on the domain boundary [0, shape.x - 1] x [0,
// shape.y - 1] Returns side: 0=bottom(y=0), 1=right(x=xmax), 2=top(y=ymax),
// 3=left(x=0), -1=interior
int get_border_side(const glm::vec2 &p, glm::ivec2 shape, float eps = 1e-3f)
{
  const float xmax = (float)(shape.x - 1);
  const float ymax = (float)(shape.y - 1);
  if (std::abs(p.y) <= eps) return 0;
  if (std::abs(p.x - xmax) <= eps) return 1;
  if (std::abs(p.y - ymax) <= eps) return 2;
  if (std::abs(p.x) <= eps) return 3;
  return -1;
}

// Perimeter coordinate along rectangle perimeter in clockwise order:
// Bottom edge (y=0, x increasing from 0 to xmax): (x, 0) -> x
// Right edge (x=xmax, y increasing from 0 to ymax): (xmax, y) -> xmax + y
// Top edge (y=ymax, x decreasing from xmax to 0): (x, ymax) -> xmax + ymax +
// (xmax - x) Left edge (x=0, y decreasing from ymax to 0): (0, y) -> 2*xmax +
// ymax + (ymax - y)
float perimeter_pos_cw(const glm::vec2 &p, glm::ivec2 shape)
{
  const float xmax = (float)(shape.x - 1);
  const float ymax = (float)(shape.y - 1);
  const float x = std::clamp(p.x, 0.f, xmax);
  const float y = std::clamp(p.y, 0.f, ymax);

  int side = get_border_side(p, shape);
  if (side == 0) return x;
  if (side == 1) return xmax + y;
  if (side == 2) return xmax + ymax + (xmax - x);
  if (side == 3) return 2.f * xmax + ymax + (ymax - y);

  // If not exactly on border, snap to nearest border
  float d_bot = y;
  float d_right = xmax - x;
  float d_top = ymax - y;
  float d_left = x;
  float min_d = std::min({d_bot, d_right, d_top, d_left});
  if (min_d == d_bot) return x;
  if (min_d == d_right) return xmax + y;
  if (min_d == d_top) return xmax + ymax + (xmax - x);
  return 2.f * xmax + ymax + (ymax - y);
}

// Extended Corner point at perimeter coordinate (slightly outside [0, xmax]x[0,
// ymax] so that boundary pixel centers (e.g. y = ymax + 0.5 - 0.5) are strictly
// enclosed by fill_polygon)
glm::vec2 corner_at_cw(int corner_idx, glm::ivec2 shape, float margin = 1.0f)
{
  const float xmax = (float)(shape.x - 1);
  const float ymax = (float)(shape.y - 1);
  switch (corner_idx)
  {
  case 0: return {xmax + margin, -margin};
  case 1: return {xmax + margin, ymax + margin};
  case 2: return {-margin, ymax + margin};
  case 3: return {-margin, -margin};
  default: return {-margin, -margin};
  }
}

// Close an open polyline along the domain boundary by taking the shorter
// boundary path
std::vector<glm::vec2> close_contour_along_boundary(
    const std::vector<glm::vec2> &pts,
    glm::ivec2                    shape)
{
  if (pts.size() < 2) return pts;

  const glm::vec2 &p_start = pts.front();
  const glm::vec2 &p_end = pts.back();

  // If already closed, return as is
  if (glm::distance(p_start, p_end) < 1e-2f) return pts;

  const float xmax = (float)(shape.x - 1);
  const float ymax = (float)(shape.y - 1);
  const float total_peri = 2.f * (xmax + ymax);

  float s0 = perimeter_pos_cw(p_end, shape);
  float s1 = perimeter_pos_cw(p_start, shape);

  // Helper to extend boundary points slightly outside the grid
  auto extend_point = [xmax, ymax](const glm::vec2 &p,
                                   float            margin = 1.0f) -> glm::vec2
  {
    glm::vec2 ext = p;
    if (ext.x <= 0.01f) ext.x = -margin;
    if (ext.x >= xmax - 0.01f) ext.x = xmax + margin;
    if (ext.y <= 0.01f) ext.y = -margin;
    if (ext.y >= ymax - 0.01f) ext.y = ymax + margin;
    return ext;
  };

  // Corner perimeter positions:
  // 0: (xmax, 0) at s = xmax
  // 1: (xmax, ymax) at s = xmax + ymax
  // 2: (0, ymax) at s = 2*xmax + ymax
  // 3: (0, 0) at s = 0 / total_peri
  const float corner_pos[4] = {xmax, xmax + ymax, 2.f * xmax + ymax, 0.f};

  // Construct CW closing path (from s0 increasing to s1 modulo total_peri)
  std::vector<glm::vec2> pts_cw = pts;
  pts_cw.push_back(extend_point(p_end));
  float dist_cw = s1 >= s0 ? (s1 - s0) : (s1 + total_peri - s0);
  std::vector<std::pair<float, int>> cw_corners;
  for (int c = 0; c < 4; ++c)
  {
    float cp = (c == 3 && s0 > 0.f) ? total_peri : corner_pos[c];
    float d = cp >= s0 ? (cp - s0) : (cp + total_peri - s0);
    if (d > 1e-2f && d < dist_cw - 1e-2f) cw_corners.emplace_back(d, c);
  }
  std::sort(cw_corners.begin(), cw_corners.end());
  for (const auto &p : cw_corners)
    pts_cw.push_back(corner_at_cw(p.second, shape));
  pts_cw.push_back(extend_point(p_start));

  // Construct CCW closing path (from s0 decreasing to s1 modulo total_peri)
  std::vector<glm::vec2> pts_ccw = pts;
  pts_ccw.push_back(extend_point(p_end));
  float       dist_ccw = s0 >= s1 ? (s0 - s1) : (s0 + total_peri - s1);
  const float ccw_corner_pos[4] = {xmax,
                                   xmax + ymax,
                                   2.f * xmax + ymax,
                                   total_peri};
  std::vector<std::pair<float, int>> ccw_corners;
  for (int c = 0; c < 4; ++c)
  {
    float cp = ccw_corner_pos[c];
    // Distance moving backward (decreasing s) from s0 to cp
    float d = s0 >= cp ? (s0 - cp) : (s0 + total_peri - cp);
    if (d > 1e-2f && d < dist_ccw - 1e-2f) ccw_corners.emplace_back(d, c);
  }
  std::sort(ccw_corners.begin(), ccw_corners.end());
  for (const auto &p : ccw_corners)
    pts_ccw.push_back(corner_at_cw(p.second, shape));
  pts_ccw.push_back(extend_point(p_start));

  // Pick whichever creates the smaller enclosed polygon area on grid
  auto poly_area = [](const std::vector<glm::vec2> &p)
  {
    float  a = 0.f;
    size_t m = p.size();
    for (size_t i = 0; i < m; ++i)
    {
      const glm::vec2 &p1 = p[i];
      const glm::vec2 &p2 = p[(i + 1) % m];
      a += (p1.x * p2.y - p2.x * p1.y);
    }
    return std::abs(a) * 0.5f;
  };

  return (poly_area(pts_cw) <= poly_area(pts_ccw)) ? pts_cw : pts_ccw;
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

    for (const glm::ivec2 &c :
         rasterize_outline(pts, shape, contours[k].is_closed()))
    {
      const int p = c.y * shape.x + c.x;
      r.contour_of[p] = k;
      if (rep[k] < 0) rep[k] = p;
    }

    if (!contours[k].is_closed())
    {
      pts = close_contour_along_boundary(pts, shape);
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
std::vector<float> passage_times(const Array  *p_probability,
                                 bool          descending,
                                 float         randomness,
                                 std::mt19937 &gen,
                                 size_t        npix)
{
  std::exponential_distribution<float> expo(1.f);
  std::vector<float>                   tau(npix);

  for (size_t p = 0; p < npix; ++p)
  {
    float rate = p_probability
                     ? std::clamp((*p_probability)((int)p), 1e-3f, 0.999f)
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
  std::vector<char>  reached; // whether the front reached the pixel
  float              tmax = 0.f;
};

// Dijkstra front propagation from the fixed pixels of contour `source` into
// the free pixels of zone `zone_id`.
Front propagate_front(const ContourRaster      &r,
                      int                       zone_id,
                      int                       source,
                      const std::vector<float> &tau)
{
  const glm::ivec2 shape = r.shape;
  const size_t     npix = r.zone.size();

  Front f;
  f.t.assign(npix, T_UNREACHED);
  f.reached.assign(npix, 0);

  using Item = std::pair<float, int>;
  std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

  auto relax = [&](int p, float t0)
  {
    const int i = p % shape.x;
    const int j = p / shape.x;
    for (int m = 0; m < 8; ++m)
    {
      const int ii = i + DI[m];
      const int jj = j + DJ[m];
      if (ii < 0 || ii >= shape.x || jj < 0 || jj >= shape.y) continue;

      const int q = jj * shape.x + ii;
      if (r.zone[q] != zone_id) continue;

      const float t = t0 + STEP[m] * tau[q];
      if (t < f.t[q])
      {
        f.t[q] = t;
        f.reached[q] = 1;
        pq.emplace(t, q);
      }
    }
  };

  for (size_t p = 0; p < npix; ++p)
    if (r.contour_of[p] == source) relax((int)p, 0.f);

  std::vector<char> settled(npix, 0);
  while (!pq.empty())
  {
    const auto [t, p] = pq.top();
    pq.pop();
    if (settled[p] || t > f.t[p]) continue;
    settled[p] = 1;
    f.tmax = std::max(f.tmax, t);
    relax(p, t);
  }

  return f;
}

// Priority flood of the elevation inside one zone, seeded from the free
// pixels adjacent to the zone's own contour. Removes pits by replacing the
// elevation with the running maximum along the flood.
void fill_pits(const ContourRaster &r, int zone_id, std::vector<float> &z)
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
      if (settled[q] || r.zone[q] != zone_id) continue;
      pq.emplace(std::max(z[q], key), q);
    }
  };

  for (size_t p = 0; p < npix; ++p)
    if (r.contour_of[p] == zone_id) push_neighbours((int)p, z[p]);

  while (!pq.empty())
  {
    const auto [key, p] = pq.top();
    pq.pop();
    if (settled[p]) continue;
    settled[p] = 1;
    z[p] = key;
    push_neighbours(p, key);
  }
}

/// Core synthesis from rasterised contours
Array synthesize_from_contour_raster(const ContourRaster      &r,
                                     const std::vector<float> &elevations,
                                     const Array              *p_probability,
                                     float                     randomness,
                                     std::uint32_t             seed,
                                     float                     peak_ratio,
                                     float                     outside_ratio,
                                     bool                      smoothstep)
{
  const glm::ivec2 shape = r.shape;
  const int        n = (int)elevations.size();
  const size_t     npix = (size_t)shape.x * shape.y;

  std::vector<std::vector<int>> children(n);
  std::vector<int>              roots;
  for (int k = 0; k < n; ++k)
  {
    if (r.parent[k] >= 0)
      children[r.parent[k]].push_back(k);
    else
      roots.push_back(k);
  }

  // pixels of each zone (index shifted by one: outside zone -1 -> 0)
  std::vector<std::vector<int>> zone_pixels(n + 1);
  for (size_t p = 0; p < npix; ++p)
    if (r.zone[p] >= -1) zone_pixels[r.zone[p] + 1].push_back((int)p);

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
  const float floor = elev_min - outside_ratio * spacing;

  // --- passage times
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

  // --- per-zone synthesis
  std::vector<float> z(npix, 0.f);

  for (size_t p = 0; p < npix; ++p)
    if (r.zone[p] == -2) z[p] = elevations[r.contour_of[p]];

  constexpr float t_eps = 1e-6f;

  for (int k = -1; k < n; ++k)
  {
    const std::vector<int> &pixels = zone_pixels[k + 1];
    if (pixels.empty()) continue;

    const std::vector<int> &kids = k < 0 ? roots : children[k];

    // one front per source: the zone's own contour (if any) and its children
    Front              own;
    std::vector<Front> far(kids.size());

    if (k >= 0) own = propagate_front(r, k, k, tau_near);
    for (size_t c = 0; c < kids.size(); ++c)
      far[c] = propagate_front(r, k, kids[c], tau_far);

    for (int p : pixels)
    {
      if (k < 0)
      {
        // outside zone: inverse-travel-time blend of the root contours, each
        // falling off with distance towards the floor
        float num = 0.f, den = 0.f;
        for (size_t c = 0; c < kids.size(); ++c)
          if (far[c].reached[p])
          {
            const float t = far[c].t[p];
            const float w = 1.f / std::max(t, t_eps);
            const float u = far[c].tmax > 0.f ? t / far[c].tmax : 0.f;
            const float factor = smoothstep ? smoothstep3(u) : u;
            const float zc = elevations[kids[c]] -
                             outside_ratio * spacing * factor;
            num += w * zc;
            den += w;
          }
        z[p] = den > 0.f ? num / den : floor;
      }
      else if (kids.empty())
      {
        // leaf zone: rise to a peak (or sink to a basin bottom)
        const bool up = r.parent[k] < 0 ||
                        elevations[k] >= elevations[r.parent[k]];
        const float delta = (up ? 1.f : -1.f) * peak_ratio * spacing;
        z[p] = elevations[k];
        if (own.reached[p] && own.tmax > 0.f)
        {
          const float u = own.t[p] / own.tmax;
          const float factor = smoothstep ? smoothstep3(u) : u;
          z[p] += delta * factor;
        }
      }
      else
      {
        // ring zone: blend own contour with child contours
        float num = 0.f, den = 0.f;
        if (own.reached[p])
        {
          const float w = 1.f / std::max(own.t[p], t_eps);
          num += w * elevations[k];
          den += w;
        }
        for (size_t c = 0; c < kids.size(); ++c)
          if (far[c].reached[p])
          {
            const float w = 1.f / std::max(far[c].t[p], t_eps);
            num += w * elevations[kids[c]];
            den += w;
          }

        float val = den > 0.f ? num / den : elevations[k];
        if (smoothstep && kids.size() == 1 && own.reached[p] &&
            far[0].reached[p])
        {
          const float t0 = own.t[p];
          const float t1 = far[0].t[p];
          const float tot = t0 + t1;
          const float u = tot > 0.f
                              ? t0 / tot
                              : 0.f; // u=0 at outer (own), u=1 at inner (child)
          val = (1.f - smoothstep3(u)) * elevations[k] +
                smoothstep3(u) * elevations[kids[0]];
        }
        z[p] = val;
      }
    }

    // pit removal in rising ring zones (all children higher than the contour)
    if (k >= 0 && !kids.empty())
    {
      bool rising = true;
      for (int c : kids)
        rising = rising && elevations[c] > elevations[k];
      if (rising) fill_pits(r, k, z);
    }
  }

  Array out(shape);
  out.vector = z;
  return out;
}

// Build ContourRaster from an Array with rasterized contours where contour
// pixels have their elevation and non-contour pixels are 0
ContourRaster rasterize_contours_from_array(const Array        &contours,
                                            std::vector<float> &elevations)
{
  const glm::ivec2 shape = contours.shape;
  const size_t     npix = (size_t)shape.x * shape.y;

  // 1. Group 8-connected non-zero pixels into distinct contour components
  std::vector<int>   contour_id(npix, -1);
  std::vector<float> comp_elev;
  int                num_contours = 0;

  for (size_t p = 0; p < npix; ++p)
  {
    if (contours((int)p) == 0.f || contour_id[p] >= 0) continue;

    const float target_elev = contours((int)p);
    const int   id = num_contours++;
    comp_elev.push_back(target_elev);

    std::vector<int> q;
    contour_id[p] = id;
    q.push_back((int)p);

    size_t head = 0;
    while (head < q.size())
    {
      int curr = q[head++];
      int cx = curr % shape.x;
      int cy = curr / shape.x;

      for (int m = 0; m < 8; ++m)
      {
        int nx = cx + DI[m];
        int ny = cy + DJ[m];
        if (nx < 0 || nx >= shape.x || ny < 0 || ny >= shape.y) continue;

        int nidx = ny * shape.x + nx;
        if (contours(nidx) != 0.f && contour_id[nidx] < 0)
        {
          contour_id[nidx] = id;
          q.push_back(nidx);
        }
      }
    }
  }

  elevations = comp_elev;
  const int n = num_contours;

  ContourRaster r;
  r.shape = shape;
  r.contour_of.assign(npix, -1);
  r.zone.assign(npix, -1);
  r.parent.assign(n, -1);

  if (n == 0) return r;

  for (size_t p = 0; p < npix; ++p)
    if (contour_id[p] >= 0) r.contour_of[p] = contour_id[p];

  // 2. Interior filling for each contour component: find connected components
  // of pixels != k. The component containing the most boundary pixels (or
  // largest area) is the exterior; all other components are inside.
  std::vector<std::vector<char>> inside(n, std::vector<char>(npix, 0));
  std::vector<int>               area(n, 0);
  std::vector<int>               rep(n, -1);

  for (size_t p = 0; p < npix; ++p)
  {
    int id = contour_id[p];
    if (id >= 0 && rep[id] < 0) rep[id] = (int)p;
  }

  for (int k = 0; k < n; ++k)
  {
    std::vector<int> comp(npix, -1);
    std::vector<int> comp_sizes;
    std::vector<int> comp_border_counts;
    int              num_comps = 0;

    for (size_t p = 0; p < npix; ++p)
    {
      if (contour_id[p] == k || comp[p] >= 0) continue;

      int cid = num_comps++;
      int sz = 0;
      int border_cnt = 0;

      std::vector<int> q;
      comp[p] = cid;
      q.push_back((int)p);

      size_t head = 0;
      while (head < q.size())
      {
        int curr = q[head++];
        ++sz;
        int cx = curr % shape.x;
        int cy = curr / shape.x;

        if (cx == 0 || cx == shape.x - 1 || cy == 0 || cy == shape.y - 1)
          ++border_cnt;

        constexpr int d4x[4] = {1, -1, 0, 0};
        constexpr int d4y[4] = {0, 0, 1, -1};

        for (int m = 0; m < 4; ++m)
        {
          int nx = cx + d4x[m];
          int ny = cy + d4y[m];
          if (nx < 0 || nx >= shape.x || ny < 0 || ny >= shape.y) continue;

          int nidx = ny * shape.x + nx;
          if (contour_id[nidx] != k && comp[nidx] < 0)
          {
            comp[nidx] = cid;
            q.push_back(nidx);
          }
        }
      }

      comp_sizes.push_back(sz);
      comp_border_counts.push_back(border_cnt);
    }

    if (num_comps <= 1)
    {
      // No enclosed interior or only a single connected background
      area[k] = 0;
      continue;
    }

    // Exterior component is the one touching the most domain boundary pixels,
    // breaking ties with maximum size
    int ext_comp = 0;
    for (int c = 1; c < num_comps; ++c)
    {
      if (comp_border_counts[c] > comp_border_counts[ext_comp] ||
          (comp_border_counts[c] == comp_border_counts[ext_comp] &&
           comp_sizes[c] > comp_sizes[ext_comp]))
      {
        ext_comp = c;
      }
    }

    int in_count = 0;
    for (size_t p = 0; p < npix; ++p)
    {
      if (comp[p] >= 0 && comp[p] != ext_comp)
      {
        inside[k][p] = 1;
        ++in_count;
      }
    }
    area[k] = in_count;
  }

  // 3. Parent identification by smallest containing enclosing contour
  for (int k = 0; k < n; ++k)
  {
    int smallest_parent = -1;
    int min_parent_area = std::numeric_limits<int>::max();
    int r_idx = rep[k];

    for (int j = 0; j < n; ++j)
    {
      if (j == k) continue;
      if (inside[j][r_idx])
      {
        if (area[j] < min_parent_area)
        {
          min_parent_area = area[j];
          smallest_parent = j;
        }
      }
    }
    r.parent[k] = smallest_parent;
  }

  // 4. Assign each pixel to innermost contour containing it, or -1 if outside
  // all
  for (size_t p = 0; p < npix; ++p)
  {
    if (contour_id[p] >= 0)
    {
      r.zone[p] = -2;
      continue;
    }

    int zone_id = -1;
    int min_area = std::numeric_limits<int>::max();
    for (int k = 0; k < n; ++k)
    {
      if (inside[k][p] && area[k] < min_area)
      {
        min_area = area[k];
        zone_id = k;
      }
    }
    r.zone[p] = zone_id;
  }

  return r;
}

} // namespace

// --- Public API

Array elevation_from_contours(glm::ivec2                shape,
                              const std::vector<Path>  &contours,
                              const std::vector<float> &elevations,
                              const Array              *p_probability,
                              float                     randomness,
                              std::uint32_t             seed,
                              float                     peak_ratio,
                              float                     outside_ratio,
                              bool                      smoothstep,
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
    if (contours[k].size() < 3)
    {
      log::error("elevation_from_contours: contour {} has fewer than 3 points",
                 k);
      return Array();
    }

  if (p_probability && !validate_same_shape(shape, *p_probability))
    return Array();

  randomness = std::clamp(randomness, 0.f, 1.f);

  // --- rasterization, zones and nesting tree
  const ContourRaster r = rasterize_contours(contours, shape, bbox);

  return synthesize_from_contour_raster(r,
                                        elevations,
                                        p_probability,
                                        randomness,
                                        seed,
                                        peak_ratio,
                                        outside_ratio,
                                        smoothstep);
}

Array elevation_from_contours(const Array  &contours,
                              const Array  *p_probability,
                              float         randomness,
                              std::uint32_t seed,
                              float         peak_ratio,
                              float         outside_ratio,
                              bool          smoothstep)
{
  if (!validate_non_empty(contours)) return Array();

  if (p_probability && !validate_same_shape(contours.shape, *p_probability))
    return Array();

  randomness = std::clamp(randomness, 0.f, 1.f);

  std::vector<float>  elevations;
  const ContourRaster r = rasterize_contours_from_array(contours, elevations);

  if (elevations.empty())
  {
    log::error("elevation_from_contours: no contours found in input array");
    return Array();
  }

  return synthesize_from_contour_raster(r,
                                        elevations,
                                        p_probability,
                                        randomness,
                                        seed,
                                        peak_ratio,
                                        outside_ratio,
                                        smoothstep);
}

} // namespace hmap
