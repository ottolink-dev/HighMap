/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/hydrology/drainage_basin_cell_based.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/random.hpp"

#include <format>
#include <unordered_map>

namespace hmap
{

DrainageBasinCellBased::DrainageBasinCellBased(const Array &z_) : z(z_)
{
  const glm::ivec2 &shape = z.shape;

  // outlets == convex hull by default
  this->outlets_mask = Mat<int>(shape, 0);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      if (i == 0 || i == shape.x - 1 || j == 0 || j == shape.y - 1)
        this->outlets_mask(i, j) = 1;
}

void DrainageBasinCellBased::accumulate_area_by_outlet(Array &acc) const
{
  if (!validate_non_empty(acc)) return;

  for (const glm::ivec2 &o : this->get_outlets())
  {
    const auto &traversal = this->traversals.at(o);

    // downstream => upstream
    for (auto it = traversal.begin(); it != traversal.end(); ++it)
    {
      const glm::ivec2 &v = *it;

      acc(v) = 1.f;

      for (const glm::ivec2 &child : this->children(v))
        acc(v) += acc(child);
    }
  }
}

void DrainageBasinCellBased::compute_receivers(unsigned int seed,
                                               float        noise_strength)
{
  const glm::ivec2 &shape = z.shape;
  const int         rows = shape.x;
  const int         cols = shape.y;

  this->receivers = Mat<glm::ivec2>(shape);

  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
    {
      if (this->outlets_mask(i, j))
      {
        this->receivers(i, j) = glm::ivec2(i, j); // self-receiver
        continue;
      }

      float      best_score = 0.f;
      glm::ivec2 best_ij(i, j);

      for (int k = 0; k < 8; ++k)
      {
        int ni = i + di[k];
        int nj = j + dj[k];

        if (ni < 0 || ni >= rows || nj < 0 || nj >= cols) continue;

        float dz = z(i, j) - z(ni, nj);
        if (dz <= 0.f) continue;

        float slope = dz * cd[k];

        // deterministic noise in [-1, 1)
        float noise = 2.f * fast_hash32_to_unit_float(seed,
                                                      (i * cols + j) ^
                                                          (ni * cols + nj)) -
                      1.f;
        float score = slope * (1.f + noise_strength * noise);

        if (score > best_score)
        {
          best_score = score;
          best_ij = glm::ivec2(ni, nj);
        }
      }

      this->receivers(i, j) = best_ij;
    }
}

void DrainageBasinCellBased::compute_receivers_priority_flood()
{
  const glm::ivec2 shape = z.shape;

  this->receivers = Mat<glm::ivec2>(shape, this->null_cell);
  Mat<int> basin_id(shape, -1);

  struct Node
  {
    int   i, j;
    float z;
  };

  struct NodeCmp
  {
    bool operator()(const Node &a, const Node &b) const
    {
      return a.z > b.z; // min-heap
    }
  };

  std::priority_queue<Node, std::vector<Node>, NodeCmp> pq;

  int basin_counter = 0;

  // --- initialize with outlets

  for (const auto &o : this->get_outlets())
  {
    basin_id(o) = basin_counter++;
    pq.push({o.x, o.y, z(o)});

    // boundary drains to itself
    this->receivers(o) = o;
  }

  // --- flood inward

  while (!pq.empty())
  {
    Node c = pq.top();
    pq.pop();

    int bc = basin_id(c.i, c.j);

    for (int k = 0; k < 8; ++k)
    {
      int ni = c.i + di[k];
      int nj = c.j + dj[k];

      if (ni < 0 || nj < 0 || ni >= shape.x || nj >= shape.y) continue;

      if (basin_id(ni, nj) == -1)
      {
        basin_id(ni, nj) = bc;

        float zn = z(ni, nj);
        float zf = std::max(zn, c.z); // flood level

        pq.push({ni, nj, zf});

        // FLOW: neighbor drains to current cell
        this->receivers(ni, nj) = glm::ivec2(c.i, c.j);
      }
    }
  }
}

Array DrainageBasinCellBased::compute_response_times(const Array &area_acc,
                                                     const Array &erodibility,
                                                     float        m_exp) const
{
  const glm::ivec2 &shape = z.shape;
  Array             response_times(shape, 0.f);

  for (const auto &[outlet, traversal] : traversals)
  {
    // traversal is stored upstream -> downstream (rbegin is outlet -> upstream)
    for (auto it = traversal.rbegin(); it != traversal.rend(); ++it)
    {
      const glm::ivec2 &i = *it;
      const glm::ivec2 &j = receivers(i);

      if (j != i)
      {
        int   dx = i.x - j.x;
        int   dy = i.y - j.y;
        float distance = (dx != 0 && dy != 0) ? M_SQRT2 : 1.f;

        float celerity = erodibility(i) *
                         std::pow(std::max(area_acc(i), 1e-8f), m_exp);

        response_times(i) = response_times(j) + distance / celerity;
      }
      else
      {
        response_times(i) = 0.f; // outlet
      }
    }
  }

  return response_times;
}

std::vector<std::vector<glm::ivec2>> DrainageBasinCellBased::
    compute_upstream_traversals()
{
  const glm::ivec2 &shape = z.shape;
  const int         rows = shape.x;
  const int         cols = shape.y;

  std::vector<std::vector<glm::ivec2>> traversals;

  for (const glm::ivec2 &outlet : this->get_outlets())
  {
    std::vector<glm::ivec2> traversal;
    traversal.reserve(rows * cols); // upper bound

    // BFS from this outlet
    traversal.push_back(outlet);

    size_t i = 0;
    while (i < traversal.size())
    {
      const glm::ivec2 &node = traversal[i];

      for (const glm::ivec2 &child : this->children(node.x, node.y))
        traversal.push_back(child);

      ++i;
    }

    traversals.push_back(std::move(traversal));
  }

  return traversals;
}

std::pair<Mat<glm::ivec2>, bool> DrainageBasinCellBased::find_subroots()
{
  const glm::ivec2 &shape = z.shape;
  const int         rows = shape.x;
  const int         cols = shape.y;

  Mat<glm::ivec2> subroot(shape, this->null_cell);
  bool            has_lake = false;

  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
    {
      if (subroot(i, j) != this->null_cell) continue;

      std::vector<glm::ivec2> path;
      glm::ivec2              p(i, j);

      while (subroot(p) == this->null_cell && receivers(p) != p &&
             receivers(p) != this->null_cell && receivers(p).x >= 0 &&
             receivers(p).x < rows && receivers(p).y >= 0 &&
             receivers(p).y < cols)
      {
        path.push_back(p);
        p = receivers(p);
      }

      glm::ivec2 root;

      if (subroot(p) != this->null_cell)
      {
        // already assigned
        root = subroot(p);
      }
      else if (outlets_mask(p))
      {
        // outlet
        root = p;
      }
      else
      {
        // true lake
        has_lake = true;
        root = p;
      }

      for (const auto &v : path)
        subroot(v) = root;

      subroot(p) = root;
    }

  return {subroot, has_lake};
}

void DrainageBasinCellBased::flow_breach()
{
  for (const auto &[outlet, traversal] : traversals)
  {
    if (traversal.empty()) continue;

    // traversal is upstream -> downstream
    glm::ivec2 p = traversal.front();
    float      zc = z(p);

    // follow receivers downstream to outlet
    while (true)
    {
      const glm::ivec2 &r = receivers(p.x, p.y);

      if (r == p) break; // reached outlet

      if (z(p) < zc)
        zc = z(p);
      else
        z(p) = zc;

      p = r;
    }
  }
}

std::vector<std::vector<glm::ivec2>> DrainageBasinCellBased::get_main_channels()
    const
{
  std::vector<std::vector<glm::ivec2>> channels;
  channels.reserve(traversals.size());

  for (const auto &[outlet, traversal] : traversals)
  {
    if (traversal.empty())
    {
      channels.emplace_back();
      continue;
    }

    std::vector<glm::ivec2> channel;
    channel.reserve(traversal.size());

    // traversal is upstream -> downstream
    glm::ivec2 p = traversal.front();
    channel.push_back(p);

    // follow receivers downstream to outlet
    while (true)
    {
      const glm::ivec2 &r = receivers(p.x, p.y);

      if (r == p) break; // reached outlet

      p = r;
      channel.push_back(p);
    }

    channels.push_back(std::move(channel));
  }

  return channels;
}

std::vector<glm::ivec2> DrainageBasinCellBased::get_outlets() const
{
  const glm::ivec2 &shape = z.shape;

  std::vector<glm::ivec2> outlet_indices;
  outlet_indices.reserve(shape.x * shape.y); // overkill

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      if (this->outlets_mask(i, j)) outlet_indices.push_back({i, j});

  return outlet_indices;
}

const Array &DrainageBasinCellBased::get_z() const
{
  return this->z;
}

void DrainageBasinCellBased::remove_lakes(const Mat<glm::ivec2> &subroot)
{
  const glm::ivec2 &shape = z.shape;
  const int         rows = shape.x;
  const int         cols = shape.y;

  Mat<uint8_t> visited(shape, 0);
  Mat<float>   dist(shape, std::numeric_limits<float>::max());
  this->roots = Mat<glm::ivec2>(shape, this->null_cell);

  struct HeapNode
  {
    float      dist;
    glm::ivec2 cell;

    bool operator>(const HeapNode &o) const noexcept
    {
      return dist > o.dist;
    }
  };

  std::vector<HeapNode> heap_storage;
  heap_storage.reserve(rows * cols * 2);

  std::priority_queue<HeapNode, std::vector<HeapNode>, std::greater<HeapNode>>
      heap(std::greater<HeapNode>(), std::move(heap_storage));

  size_t remaining_lakes = 0;
  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
      if (subroot(i, j) != this->null_cell) remaining_lakes++;

  // initialize outlets
  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
    {
      if (outlets_mask(i, j))
      {
        roots(i, j) = glm::ivec2(i, j);
        dist(i, j) = 0.f;
        heap.push({0.f, glm::ivec2(i, j)});
      }
    }

  while (!heap.empty())
  {
    HeapNode top = heap.top();
    heap.pop();

    glm::ivec2 p = top.cell;
    if (visited(p)) continue;

    visited(p) = 1;

    for (int k = 0; k < 8; ++k)
    {
      int ni = p.x + di[k];
      int nj = p.y + dj[k];

      if (ni < 0 || ni >= rows || nj < 0 || nj >= cols) continue;
      if (visited(ni, nj)) continue;

      glm::ivec2 sr = subroot(ni, nj);

      if (sr != this->null_cell && roots(sr) == this->null_cell)
      {
        // redirect lake receivers to outlet along the path
        glm::ivec2 k_cell(ni, nj);
        glm::ivec2 nk = p;

        while (receivers(k_cell) != k_cell)
        {
          glm::ivec2 tmp = receivers(k_cell);
          receivers(k_cell) = nk;
          nk = k_cell;
          k_cell = tmp;
        }
        receivers(k_cell) = nk;

        roots(sr) = roots(subroot(p));

        if (--remaining_lakes == 0) return;
      }

      // Euclidean distance in grid units
      float new_dist = top.dist + cd[k] + std::abs(z(ni, nj) - z(p)) * cd[k];
      if (new_dist < dist(ni, nj))
      {
        dist(ni, nj) = new_dist;
        heap.push({new_dist, glm::ivec2(ni, nj)});
      }
    }

    glm::ivec2 sri = subroot(p);
    if (sri != this->null_cell) roots(p) = roots(sri);
  }
}

void DrainageBasinCellBased::set_outlets(
    const std::vector<glm::ivec2> &outlet_indices)
{
  const glm::ivec2 &shape = z.shape;
  this->outlets_mask = Mat<int>(shape, 0);

  for (const auto &o : outlet_indices)
    this->outlets_mask(o) = 1;
}

float DrainageBasinCellBased::update_elevations(const Array &response_times,
                                                float        uplift_rate,
                                                const Array &max_slope)
{
  const glm::ivec2 &shape = z.shape;

  float delta_sum = 0.f;

  // compute z range
  glm::vec2 zr = this->z.range();
  float     zmin = zr.x;
  float     zmax = zr.y;
  float     zptp = zmax - zmin;

  // iterate upstream (outlet first, using reverse traversal since traversal is
  // upstream->downstream)
  for (const auto &[outlet, traversal] : traversals)
  {
    for (auto it = traversal.rbegin(); it != traversal.rend(); ++it)
    {
      const glm::ivec2 &i = *it;
      const glm::ivec2 &j = receivers(i);

      if (j == i) continue; // outlet

      float dt = std::max(response_times(i) - response_times(outlet), 0.f);
      float new_elevation = z(outlet) + uplift_rate * dt;

      // slope limiting
      int   dx = i.x - j.x;
      int   dy = i.y - j.y;
      float distance = (dx != 0 && dy != 0) ? M_SQRT2 : 1.f;

      float slope = (new_elevation - z(j)) / distance;
      float max_slope_n = max_slope(i) * zptp;

      if (slope > max_slope_n) new_elevation = z(j) + max_slope_n * distance;

      float delta = std::abs(new_elevation - z(i));
      delta_sum += delta;

      z(i) = new_elevation;
    }
  }

  delta_sum /= float(shape.x * shape.y);

  return delta_sum;
}

void DrainageBasinCellBased::update_stream_tree(unsigned int seed,
                                                float        noise_strength)
{
  this->compute_receivers(seed, noise_strength);
  auto [subroots, has_lake] = this->find_subroots();

  if (has_lake) this->remove_lakes(subroots);

  this->update_traversals();
}

void DrainageBasinCellBased::update_stream_tree()
{
  this->update_stream_tree(0, 0.f);
}

void DrainageBasinCellBased::update_traversals()
{
  this->traversals.clear();
  this->children = invert_receiver_map(this->receivers);

  for (const glm::ivec2 &o : this->get_outlets())
  {
    std::vector<glm::ivec2> traversal;
    traversal.reserve(z.shape.x * z.shape.y);
    traversal.push_back(o);

    size_t i = 0;
    while (i < traversal.size())
    {
      const glm::ivec2 &node = traversal[i];

      // add upstream nodes
      for (const glm::ivec2 &child : this->children(node.x, node.y))
        traversal.push_back(child);

      ++i;
    }

    std::reverse(traversal.begin(), traversal.end());
    this->traversals[o] = std::move(traversal);
  }
}

// --- FUNCTIONS

Mat<std::vector<glm::ivec2>> invert_receiver_map(
    const Mat<glm::ivec2> &receivers)
{
  const glm::ivec2 &shape = receivers.shape;
  const int         rows = shape.x;
  const int         cols = shape.y;

  Mat<std::vector<glm::ivec2>> children(shape);

  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
    {
      const glm::ivec2 &r = receivers(i, j);

      if (r.x >= 0 && r.x < rows && r.y >= 0 && r.y < cols)
      {
        if (r != glm::ivec2(i, j))
          children(r.x, r.y).push_back(glm::ivec2(i, j));
      }
    }

  return children;
}

} // namespace hmap
