/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "hmm/src/heightmap.h"
#include "hmm/src/triangulator.h"

#include "highmap/array.hpp"
#include "highmap/hydrology/drainage_basin.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/random.hpp"
#include "highmap/terrain_tri_mesh.hpp"

#include <unordered_map>

namespace hmap
{

DrainageBasin::DrainageBasin(std::vector<glm::vec3> xyz_)
    : mesh(std::move(xyz_))
{
  // outlets == convex hull by default
  this->outlets_mask.assign(this->mesh.size(), false);

  for (const auto &idx : this->mesh.get_convex_hull())
    this->outlets_mask[idx] = true;
}

void DrainageBasin::accumulate_area_by_outlet(const std::vector<float> &area,
                                              std::vector<float> &acc) const
{
  const auto &outlets = this->get_outlets();
  const int   n_outlets = static_cast<int>(outlets.size());

#pragma omp parallel for schedule(static)
  for (int oi = 0; oi < n_outlets; ++oi)
  {
    size_t o = outlets[oi];
    auto   it_t = this->traversals.find(o);
    if (it_t == this->traversals.end()) continue;

    const auto &traversal = it_t->second;

    for (size_t v : traversal)
    {
      acc[v] = area[v];

      for (size_t child : this->children[v])
        acc[v] += acc[child];
    }
  }
}

std::vector<bool> DrainageBasin::compute_is_ridge_node() const
{
  std::vector<bool>                   is_ridge(this->mesh.size(), false);
  const TerrainTriMesh::NeighborData &nbrs_data = this->mesh.get_neighbors();

  if (this->roots.size() != this->mesh.size()) return is_ridge;

  for (size_t k = 0; k < this->mesh.size(); ++k)
  {
    for (const auto &nb : nbrs_data.adjacency[k])
    {
      size_t n = nb.index;

      if (this->roots[k] != this->roots[n])
      {
        is_ridge[k] = true;
      }
    }
  }

  return is_ridge;
}

void DrainageBasin::compute_receivers()
{
  const size_t n = this->mesh.size();
  this->receivers.resize(n);

  const TerrainTriMesh::NeighborData &nbrs_data = this->mesh.get_neighbors();
  const auto                         &points = this->mesh.get_points();

#pragma omp parallel for schedule(static)
  for (int k = 0; k < int(n); ++k)
  {
    if (this->outlets_mask[k])
    {
      this->receivers[k] = static_cast<size_t>(k);
      continue;
    }

    const float zk = points[k].z;
    float       best_dz = 0.f;
    float       best_dist = 1.f;
    size_t      best_k = static_cast<size_t>(k);

    for (const auto &nb : nbrs_data.adjacency[k])
    {
      const float dz = zk - points[nb.index].z;
      if (dz > 0.f && dz * best_dist > best_dz * nb.distance2d)
      {
        best_dz = dz;
        best_dist = nb.distance2d;
        best_k = nb.index;
      }
    }

    this->receivers[k] = best_k;
  }
}

void DrainageBasin::compute_receivers(unsigned int seed, float noise_strength)
{
  const size_t n = this->mesh.size();
  this->receivers.resize(n);

  const TerrainTriMesh::NeighborData &nbrs_data = this->mesh.get_neighbors();
  const auto                         &points = this->mesh.get_points();

#pragma omp parallel for schedule(static)
  for (int k = 0; k < int(n); ++k)
  {
    if (this->outlets_mask[k])
    {
      this->receivers[k] = static_cast<size_t>(k);
      continue;
    }

    float       best_score = -std::numeric_limits<float>::infinity();
    size_t      best_k = static_cast<size_t>(k);
    const float zk = points[k].z;

    for (const auto &nb : nbrs_data.adjacency[k])
    {
      const float dz = zk - points[nb.index].z;

      if (dz > 0.f)
      {
        const float slope = dz / nb.distance2d;
        const float noise = 2.f * fast_hash32_to_unit_float(
                                      seed,
                                      static_cast<uint32_t>(k ^ nb.index)) -
                            1.f;
        const float score = slope * (1.f + noise_strength * noise);

        if (score > best_score)
        {
          best_score = score;
          best_k = nb.index;
        }
      }
    }

    this->receivers[k] = best_k;
  }
}

std::vector<size_t> DrainageBasin::compute_strahler_order() const
{
  const size_t        n = this->receivers.size();
  std::vector<size_t> order(n, 1);

  for (size_t o : this->get_outlets())
  {
    auto it_t = this->traversals.find(o);
    if (it_t == this->traversals.end()) continue;

    const auto &traversal = it_t->second;

    // traversal is upstream -> downstream (leaves first)
    for (size_t i : traversal)
    {
      if (this->children[i].empty())
      {
        order[i] = 1;
        continue;
      }

      size_t max_child_order = 0;
      size_t count_max = 0;

      for (size_t child : this->children[i])
      {
        size_t child_order = order[child];
        if (child_order > max_child_order)
        {
          max_child_order = child_order;
          count_max = 1;
        }
        else if (child_order == max_child_order)
        {
          count_max++;
        }
      }

      if (count_max > 1)
        order[i] = max_child_order + 1;
      else
        order[i] = max_child_order;
    }
  }

  return order;
}

std::vector<float> DrainageBasin::compute_response_times(
    const std::vector<float> &area_acc,
    const std::vector<float> &erodibility,
    float                     m_exp) const
{
  const size_t       n = this->receivers.size();
  std::vector<float> response_times(n, 0.f);

  const glm::vec3 dref = this->mesh.get_reference_lengths();
  const auto     &outlets = this->get_outlets();
  const int       n_outlets = static_cast<int>(outlets.size());
  const auto     &points = this->mesh.get_points();

#pragma omp parallel for schedule(static)
  for (int oi = 0; oi < n_outlets; ++oi)
  {
    size_t o = outlets[oi];
    auto   it_t = this->traversals.find(o);
    if (it_t == this->traversals.end()) continue;

    const auto &traversal = it_t->second;

    // downstream => upstream
    for (auto it = traversal.rbegin(); it != traversal.rend(); ++it)
    {
      size_t i = *it;
      size_t j = this->receivers[i];

      if (j != i)
      {
        const auto &pi = points[i];
        const auto &pj = points[j];
        float       vx = (pi.x - pj.x) / dref.x;
        float       vy = (pi.y - pj.y) / dref.y;
        float       distance = std::sqrt(vx * vx + vy * vy);

        float celerity = erodibility[i] *
                         std::pow(std::max(area_acc[i], 1e-8f), m_exp);

        response_times[i] = response_times[j] + distance / celerity;
      }
      else
      {
        response_times[i] = 0.f; // outlet
      }
    }
  }

  return response_times;
}

std::vector<float> DrainageBasin::compute_vertex_areas() const
{
  return this->mesh.get_vertex_areas(false);
}

std::pair<std::vector<size_t>, bool> DrainageBasin::find_subroots()
{
  const size_t n = this->receivers.size();

  std::vector<size_t> subroot(n, this->invalid_index);
  bool                has_lake = false;

  std::vector<size_t> path;
  path.reserve(64);

  for (size_t i = 0; i < n; ++i)
  {
    if (subroot[i] != this->invalid_index) continue;

    path.clear();
    size_t p = i;

    while (subroot[p] == this->invalid_index && this->receivers[p] != p)
    {
      path.push_back(p);
      p = this->receivers[p];
    }

    size_t root;

    if (subroot[p] != this->invalid_index)
    {
      // already assigned
      root = subroot[p];
    }
    else if (this->outlets_mask[p])
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

    for (size_t v : path)
      subroot[v] = root;

    subroot[p] = root;
  }

  return {subroot, has_lake};
}

void DrainageBasin::flow_breach()
{
  auto &pts = this->mesh.get_points();

  for (const auto &[outlet, traversal] : this->traversals)
  {
    if (traversal.empty()) continue;

    // traversal is upstream -> downstream
    size_t p = traversal.front();
    float  zc = pts[p].z; // current

    // follow receivers downstream to outlet
    while (true)
    {
      const size_t &r = this->receivers[p];

      if (r == p) break; // reached outlet

      if (pts[p].z < zc)
        zc = pts[p].z;
      else
        pts[p].z = zc;

      p = r;
    }
  }
}

std::vector<std::vector<glm::vec3>> DrainageBasin::flow_breach_paths()
{
  auto &pts = this->mesh.get_points();

  std::vector<std::vector<glm::vec3>> paths;
  paths.reserve(this->traversals.size());

  for (const auto &[outlet, traversal] : this->traversals)
  {
    if (traversal.empty())
    {
      paths.emplace_back();
      continue;
    }

    std::vector<glm::vec3> path;
    path.reserve(traversal.size());

    // traversal is upstream -> downstream
    size_t p = traversal.front();
    float  zc = pts[p].z; // current

    path.push_back(pts[p]);

    // follow receivers downstream to outlet
    while (true)
    {
      const size_t &r = this->receivers[p];

      if (r == p) break; // reached outlet

      if (pts[p].z < zc)
        zc = pts[p].z;
      else
        pts[p].z = zc;

      path.push_back(pts[p]);

      p = r;
    }

    paths.push_back(std::move(path));
  }

  return paths;
}

const std::vector<size_t> &DrainageBasin::for_each_upstream(size_t outlet) const
{
  return this->traversals.at(outlet);
}

std::vector<std::vector<size_t>> DrainageBasin::get_main_channels() const
{
  std::vector<std::vector<size_t>> channels;
  channels.reserve(this->traversals.size());

  for (const auto &[outlet, traversal] : this->traversals)
  {
    if (traversal.empty())
    {
      channels.emplace_back();
      continue;
    }

    std::vector<size_t> channel;
    channel.reserve(traversal.size());

    // traversal is upstream -> downstream
    size_t p = traversal.front();
    channel.push_back(p);

    // follow receivers downstream to outlet
    while (true)
    {
      const size_t &r = this->receivers[p];

      if (r == p) break; // reached outlet

      p = r;
      channel.push_back(p);
    }

    channels.push_back(std::move(channel));
  }

  return channels;
}

const TerrainTriMesh &DrainageBasin::get_mesh() const
{
  return this->mesh;
}

TerrainTriMesh &DrainageBasin::get_mesh()
{
  return this->mesh;
}

const std::vector<size_t> &DrainageBasin::get_outlets() const
{
  if (this->outlets_dirty)
  {
    this->cached_outlets.clear();

    for (size_t k = 0; k < this->mesh.size(); ++k)
      if (this->outlets_mask[k]) this->cached_outlets.push_back(k);

    this->outlets_dirty = false;
  }
  return this->cached_outlets;
}

const std::vector<size_t> &DrainageBasin::get_receivers() const
{
  return this->receivers;
}

const std::vector<glm::vec3> &DrainageBasin::get_xyz() const
{
  return this->mesh.get_points();
}

void DrainageBasin::invert_receiver_map()
{
  const size_t n = this->receivers.size();

  // rebuild in-place, reusing existing capacity
  for (auto &c : this->children)
    c.clear();
  if (this->children.size() != n) this->children.resize(n);

  for (size_t v = 0; v < n; ++v)
  {
    size_t r = this->receivers[v];
    if (r != v) this->children[r].push_back(v);
  }
}

void DrainageBasin::remap(float vmin, float vmax)
{
  if (this->mesh.get_points().empty()) return;

  float zmin = this->mesh.get_points()[0].z;
  float zmax = this->mesh.get_points()[0].z;

  for (const auto &p : this->mesh.get_points())
  {
    zmin = std::min(zmin, p.z);
    zmax = std::max(zmax, p.z);
  }

  float dz = zmax - zmin;

  if (dz == 0.f)
  {
    float mid = 0.5f * (vmin + vmax);
    for (auto &p : this->mesh.get_points())
      p.z = mid;
    return;
  }

  float scale = (vmax - vmin) / dz;

  for (auto &p : this->mesh.get_points())
    p.z = vmin + (p.z - zmin) * scale;
}

void DrainageBasin::remove_lakes(const std::vector<size_t> &subroot)
{
  const size_t n = this->receivers.size();

  std::vector<uint8_t> visited(n, 0);
  std::vector<float>   dist(n, std::numeric_limits<float>::max());

  this->roots.assign(n, this->invalid_index);

  struct HeapNode
  {
    float  dist;
    size_t node;

    bool operator>(const HeapNode &o) const noexcept
    {
      return dist > o.dist;
    }
  };

  std::vector<HeapNode> heap_storage;
  heap_storage.reserve(n);

  std::priority_queue<HeapNode, std::vector<HeapNode>, std::greater<HeapNode>>
      heap(std::greater<HeapNode>(), std::move(heap_storage));

  // count distinct lake roots (depressions)
  size_t remaining_lakes = 0;
  for (size_t i = 0; i < n; ++i)
  {
    if (subroot[i] == i && !this->outlets_mask[i]) remaining_lakes++;
  }

  // initialize outlets
  for (size_t k = 0; k < n; ++k)
  {
    if (this->outlets_mask[k])
    {
      this->roots[k] = k;
      dist[k] = 0.f;
      heap.push({0.f, k});
    }
  }

  const size_t                       *subroot_ptr = subroot.data();
  const TerrainTriMesh::NeighborData &nbrs_data = this->mesh.get_neighbors();

  while (!heap.empty())
  {
    HeapNode top = heap.top();
    heap.pop();

    size_t i = top.node;
    if (visited[i]) continue;

    visited[i] = 1;

    const auto &nbrs = nbrs_data.adjacency[i];

    for (const auto &nb : nbrs)
    {
      size_t j = nb.index;
      if (visited[j]) continue;

      size_t sr = subroot_ptr[j];

      if (sr != this->invalid_index && this->roots[sr] == this->invalid_index)
      {
        size_t k_node = j;
        size_t nk = i;

        while (this->receivers[k_node] != k_node)
        {
          size_t tmp = this->receivers[k_node];
          this->receivers[k_node] = nk;
          nk = k_node;
          k_node = tmp;
        }

        this->receivers[k_node] = nk;

        this->roots[subroot[j]] = this->roots[subroot[i]];

        if (--remaining_lakes == 0)
        {
          // assign remaining unresolved node roots from their subroots
          for (size_t v = 0; v < n; ++v)
          {
            if (this->roots[v] == this->invalid_index &&
                subroot[v] != this->invalid_index)
              this->roots[v] = this->roots[subroot[v]];
          }
          return;
        }
      }

      float new_dist = top.dist + nb.distance2d;

      if (new_dist < dist[j])
      {
        dist[j] = new_dist;
        heap.push({new_dist, j});
      }
    }

    size_t sri = subroot_ptr[i];
    if (sri != this->invalid_index && this->roots[i] == this->invalid_index)
      this->roots[i] = this->roots[sri];
  }

  // finalize any remaining node roots
  for (size_t v = 0; v < n; ++v)
  {
    if (this->roots[v] == this->invalid_index &&
        subroot[v] != this->invalid_index)
      this->roots[v] = this->roots[subroot[v]];
  }
}

void DrainageBasin::set_outlets(const std::vector<size_t> &outlet_indices)
{
  this->outlets_mask = std::vector<bool>(this->size(), false);

  for (const auto &idx : outlet_indices)
    this->outlets_mask[idx] = true;

  this->outlets_dirty = true;
}

size_t DrainageBasin::size() const
{
  return this->mesh.get_points().size();
}

void DrainageBasin::to_csv(const std::string &filename) const
{
  std::ofstream f(filename, std::ios::out);
  if (!f.is_open()) return;

  auto order = this->compute_strahler_order();
  auto area = this->compute_vertex_areas();
  auto is_ridge = this->compute_is_ridge_node();

  std::vector<float> acc(this->size(), 0.f);
  std::vector<float> flow(this->size(), 1.f);
  this->accumulate_area_by_outlet(area, acc);
  this->accumulate_area_by_outlet(area, flow);

  f << "# vertices\n";
  f << "# "
       "vertex_id,x,y,z,is_outlet,receiver,root,order,area,area_acc,flow_"
       "acc,is_ridge\n";

  for (size_t i = 0; i < static_cast<size_t>(this->mesh.get_points().size());
       ++i)
  {
    size_t r = (i < this->receivers.size()) ? this->receivers[i]
                                            : this->invalid_index;
    size_t rt = (i < this->roots.size()) ? this->roots[i] : this->invalid_index;

    f << i << "," << this->mesh.get_points()[i].x << ","
      << this->mesh.get_points()[i].y << "," << this->mesh.get_points()[i].z
      << "," << (this->outlets_mask[i] ? 1 : 0) << "," << r << "," << rt << ","
      << order[i] << "," << area[i] << "," << acc[i] << "," << flow[i] << ","
      << (is_ridge[i] ? 1 : 0) << "\n";
  }

  f.close();
}

float DrainageBasin::update_elevations(const std::vector<float> &response_times,
                                       float                     uplift_rate,
                                       const std::vector<float> &max_slope)
{
  const auto  &outlets = this->get_outlets();
  const size_t n_outlets = outlets.size();

  const glm::vec2         zr = this->mesh.get_range_z();
  const float             zptp = zr.y - zr.x;
  std::vector<glm::vec3> &points = this->mesh.get_points();

  float delta_sum = 0.f;

  // safe to parallelize over outlets: each basin owns a disjoint set
  // of nodes, so points[i].z writes never conflict across threads.

#pragma omp parallel for schedule(static) reduction(+ : delta_sum)
  for (int oi = 0; oi < int(n_outlets); ++oi)
  {
    const size_t outlet = outlets[oi];
    auto         it_t = this->traversals.find(outlet);
    if (it_t == this->traversals.end()) continue;

    const auto &traversal = it_t->second;

    // hoist per-traversal constants out of inner loop
    const float outlet_z = points[outlet].z; // not modified (outlet skipped)
    const float outlet_rt = response_times[outlet];

    // iterate downstream -> upstream (outlet outwards to leaves)
    for (auto it = traversal.rbegin(); it != traversal.rend(); ++it)
    {
      const size_t i = *it;
      const size_t j = this->receivers[i];
      if (j == i) continue; // skip outlet node

      const float dt = std::max(response_times[i] - outlet_rt, 0.f);
      float       new_elevation = outlet_z + uplift_rate * dt;

      // slope limiting
      const float dx = points[i].x - points[j].x;
      const float dy = points[i].y - points[j].y;
      const float distance = std::sqrt(dx * dx + dy * dy);
      const float slope = (new_elevation - points[j].z) / distance;
      const float max_sl = max_slope[i] * zptp;

      if (slope > max_sl) new_elevation = points[j].z + max_sl * distance;

      delta_sum += std::abs(new_elevation - points[i].z);
      points[i].z = new_elevation;
    }
  }

  return delta_sum / float(this->size());
}

void DrainageBasin::update_stream_tree(unsigned int seed, float noise_strength)
{
  this->compute_receivers(seed, noise_strength);

  auto [subroots, has_lake] = this->find_subroots();
  if (has_lake)
  {
    this->remove_lakes(subroots);
  }
  else
  {
    this->roots = std::move(subroots);
  }

  this->invert_receiver_map();

  this->update_traversals();
}

void DrainageBasin::update_stream_tree()
{
  // deactivate noise for receivers
  float         noise_strength = 0.f;
  std::uint32_t seed = 0; // dummy value

  this->update_stream_tree(seed, noise_strength);
}

void DrainageBasin::update_traversals()
{
  const std::vector<size_t> &outlets = this->get_outlets();

  this->traversals.clear();
  const size_t n_outlets = outlets.size();
  const size_t reserve_size = n_outlets > 0 ? this->mesh.size() / n_outlets
                                            : this->mesh.size();

  // pre-insert with empty vector to make sure all keys exist
  for (size_t o : outlets)
    this->traversals[o];

#pragma omp parallel for schedule(static)
  for (int oi = 0; oi < int(n_outlets); ++oi)
  {
    size_t              o = outlets[oi];
    std::vector<size_t> traversal;
    traversal.reserve(reserve_size);
    traversal.push_back(o);

    size_t i = 0;
    while (i < traversal.size())
    {
      size_t node = traversal[i];

      // add upstream nodes
      for (size_t child : this->children[node])
        traversal.push_back(child);

      ++i;
    }

    std::reverse(traversal.begin(), traversal.end());
    this->traversals[o] = std::move(traversal);
  }
}

// --- FUNCTIONS

std::vector<size_t> find_border_minima(const std::vector<glm::vec3> &xyz,
                                       float                         eps)
{
  const size_t invalid = std::numeric_limits<size_t>::max();

  size_t xmin = invalid;
  size_t xmax = invalid;
  size_t ymin = invalid;
  size_t ymax = invalid;

  float zxmin = std::numeric_limits<float>::max();
  float zxmax = std::numeric_limits<float>::max();
  float zymin = std::numeric_limits<float>::max();
  float zymax = std::numeric_limits<float>::max();

  for (size_t i = 0; i < xyz.size(); ++i)
  {
    const auto &p = xyz[i];

    if (std::abs(p.x) < eps)
    {
      if (p.z < zxmin)
      {
        zxmin = p.z;
        xmin = i;
      }
    }

    if (std::abs(p.x - 1.f) < eps)
    {
      if (p.z < zxmax)
      {
        zxmax = p.z;
        xmax = i;
      }
    }

    if (std::abs(p.y) < eps)
    {
      if (p.z < zymin)
      {
        zymin = p.z;
        ymin = i;
      }
    }

    if (std::abs(p.y - 1.f) < eps)
    {
      if (p.z < zymax)
      {
        zymax = p.z;
        ymax = i;
      }
    }
  }

  return {xmin, xmax, ymin, ymax};
}

std::vector<size_t> find_border_sinks(const TerrainTriMesh &mesh, float eps)
{
  const auto &pts = mesh.get_points();
  const auto &nbrs_data = mesh.get_neighbors();
  const auto &bbox = mesh.get_bbox();

  std::vector<size_t> sinks;

  for (size_t i = 0; i < pts.size(); ++i)
  {
    const auto &p = pts[i];

    // border test
    bool is_border = (std::abs(p.x - bbox.min.x) < eps) ||
                     (std::abs(p.x - bbox.max.x) < eps) ||
                     (std::abs(p.y - bbox.min.y) < eps) ||
                     (std::abs(p.y - bbox.max.y) < eps);

    if (!is_border) continue;

    const auto &nbrs = nbrs_data.adjacency[i];

    if (nbrs.empty()) continue;

    bool is_sink = true;

    for (const auto &nb : nbrs)
    {
      size_t j = nb.index;

      if (j >= pts.size()) continue;

      // strictly lower neighbor -> not a sink
      if (pts[j].z < p.z - eps)
      {
        is_sink = false;
        break;
      }
    }

    if (is_sink) sinks.push_back(i);
  }

  return sinks;
}

std::vector<glm::vec3> heightmap_retopology(const Array &z,
                                            float        max_error,
                                            int          max_triangles,
                                            int          max_points)
{
  if (!validate_non_empty(z)) return {};

  const glm::ivec2 &shape = z.shape;

  const auto   p_hmap = std::make_shared<Heightmap>(shape.y,
                                                  shape.x,
                                                  z.get_vector());
  Triangulator tri(p_hmap);
  tri.Run(max_error, max_triangles, max_points);

  const auto &points = tri.Points(1.f);

  // x, y normalization coefficients
  const float ax = 1.f / float(shape.y - 1);
  const float ay = 1.f / float(shape.x - 1);

  std::vector<glm::vec3> xyz;
  xyz.reserve(points.size());

  for (const auto &p : points)
    xyz.push_back({ay * p.y, ax * p.x, p.z});

  return xyz;
}

std::vector<size_t> sample_border_points(const std::vector<glm::vec3> &xyz,
                                         size_t                        nb)
{
  std::vector<size_t> indices;
  indices.reserve(nb);

  if (xyz.empty() || nb == 0) return indices;

  float perimeter = 4.f;
  float step = perimeter / float(nb);

  for (size_t k = 0; k < nb; ++k)
  {
    float s = k * step;

    float x, y;

    if (s < 1.f) // bottom edge
    {
      x = s;
      y = 0.f;
    }
    else if (s < 2.f) // right edge
    {
      x = 1.f;
      y = s - 1.f;
    }
    else if (s < 3.f) // top edge
    {
      x = 3.f - s;
      y = 1.f;
    }
    else // left edge
    {
      x = 0.f;
      y = 4.f - s;
    }

    float  best_d2 = std::numeric_limits<float>::max();
    size_t best_i = 0;

    for (size_t i = 0; i < xyz.size(); ++i)
    {
      float dx = xyz[i].x - x;
      float dy = xyz[i].y - y;
      float d2 = dx * dx + dy * dy;

      if (d2 < best_d2)
      {
        best_d2 = d2;
        best_i = i;
      }
    }

    indices.push_back(best_i);
  }

  return indices;
}

} // namespace hmap
