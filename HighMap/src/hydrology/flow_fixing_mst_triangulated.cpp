/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <future>
#include <iterator>
#include <limits>
#include <mutex>
#include <queue>
#include <utility>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/carving.hpp"
#include "highmap/filters.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/constants.hpp"
#include "highmap/math/profiles.hpp"
#include "highmap/shortest_path.hpp"
#include "highmap/terrain_tri_mesh.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <format>
#include <unordered_map>

namespace hmap
{
std::vector<Path> flow_fixing_mst_paths(const TerrainTriMesh &mesh,
                                        float                 riverbed_talus,
                                        float                 elevation_ratio,
                                        float                 distance_exponent,
                                        float upward_penalization,
                                        float minimum_depth)
{
  if (mesh.size() < 3) return {};

  const auto  &points = mesh.get_points();
  const auto  &nbrs_data = mesh.get_neighbors();
  const auto  &convex_hull = mesh.get_convex_hull();
  const size_t n_vertices = points.size();

  std::vector<bool> is_boundary(n_vertices, false);
  for (size_t idx : convex_hull)
    is_boundary[idx] = true;

  // --- Identify interior sinks

  std::vector<size_t> sinks;
  for (size_t i = 0; i < n_vertices; ++i)
  {
    if (is_boundary[i]) continue;

    bool is_sink = true;
    for (const auto &nb : nbrs_data.adjacency[i])
    {
      if (points[nb.index].z < points[i].z)
      {
        is_sink = false;
        break;
      }
    }
    if (is_sink) sinks.push_back(i);
  }

  if (sinks.empty()) return {};

  int n_sinks = static_cast<int>(sinks.size());
  int boundary_src_id = n_sinks;
  int total_sources = n_sinks + 1;

  glm::vec2 z_range = mesh.get_range_z();
  float     z_span = std::max(z_range.y - z_range.x, 1e-6f);

  // --- Multi-source Dijkstra expansion on mesh

  struct MeshMSTEdge
  {
    float               cost;
    int                 u;
    int                 v;
    std::vector<size_t> path; // vertex indices from source u to source v

    bool operator<(const MeshMSTEdge &other) const
    {
      return cost < other.cost;
    }
  };

  struct MeshDijkstraNode
  {
    float  dist;
    size_t u;

    bool operator>(const MeshDijkstraNode &o) const
    {
      return dist > o.dist;
    }
  };

  std::unordered_map<int64_t, MeshMSTEdge> candidate_edges;
  auto make_key = [](int u, int v) -> int64_t
  {
    if (u > v) std::swap(u, v);
    return (static_cast<int64_t>(u) << 32) | static_cast<int64_t>(v);
  };

  std::vector<float>  dist_map(n_vertices, std::numeric_limits<float>::max());
  std::vector<int>    owner_map(n_vertices, -1);
  std::vector<size_t> prev_node(n_vertices, size_t(-1));

  std::priority_queue<MeshDijkstraNode,
                      std::vector<MeshDijkstraNode>,
                      std::greater<MeshDijkstraNode>>
      pq;

  // initialize boundary outlets
  for (size_t idx : convex_hull)
  {
    float norm_z = (points[idx].z - z_range.x) / z_span;
    float init_d = 5.f * norm_z;
    dist_map[idx] = init_d;
    owner_map[idx] = boundary_src_id;
    prev_node[idx] = idx;
    pq.push({init_d, idx});
  }

  // initialize interior sinks
  for (int s = 0; s < n_sinks; ++s)
  {
    size_t u = sinks[s];
    dist_map[u] = 0.f;
    owner_map[u] = s;
    prev_node[u] = u;
    pq.push({0.f, u});
  }

  while (!pq.empty())
  {
    MeshDijkstraNode top = pq.top();
    pq.pop();

    size_t ci = top.u;
    if (top.dist > dist_map[ci]) continue;

    int cur_owner = owner_map[ci];

    for (const auto &nb : nbrs_data.adjacency[ci])
    {
      size_t ni = nb.index;
      int    nb_owner = owner_map[ni];

      // when meeting a different source region, record candidate bridge edge
      if (nb_owner != -1 && nb_owner != cur_owner)
      {
        float   total_cost = dist_map[ci] + dist_map[ni] + nb.distance2d;
        int64_t key = make_key(cur_owner, nb_owner);

        if (candidate_edges.find(key) == candidate_edges.end() ||
            total_cost < candidate_edges[key].cost)
        {
          std::vector<size_t> p1;
          size_t              curr = ci;
          while (true)
          {
            p1.push_back(curr);
            size_t nxt = prev_node[curr];
            if (nxt == curr) break;
            curr = nxt;
          }
          std::reverse(p1.begin(), p1.end()); // from cur_owner to ci

          std::vector<size_t> p2;
          curr = ni;
          while (true)
          {
            p2.push_back(curr);
            size_t nxt = prev_node[curr];
            if (nxt == curr) break;
            curr = nxt;
          } // from ni to nb_owner

          p1.insert(p1.end(), p2.begin(), p2.end());
          candidate_edges[key] = {total_cost,
                                  cur_owner,
                                  nb_owner,
                                  std::move(p1)};
        }
      }

      // transition cost
      float dz = points[ni].z - points[ci].z;
      float cost_step = (1.f - elevation_ratio) * nb.distance2d;

      if (dz > 0.f)
        cost_step += upward_penalization * std::pow(dz, distance_exponent);
      else
        cost_step += std::abs(dz);

      cost_step += elevation_ratio * std::max(0.f, points[ni].z);

      float new_dist = dist_map[ci] + cost_step;

      if (new_dist < dist_map[ni])
      {
        dist_map[ni] = new_dist;
        owner_map[ni] = cur_owner;
        prev_node[ni] = ci;
        pq.push({new_dist, ni});
      }
    }
  }

  // --- Build Kruskal Minimum Spanning Tree across all sinks + boundary outlet

  struct MeshDSU
  {
    std::vector<int> parent;
    MeshDSU(int n) : parent(n)
    {
      for (int i = 0; i < n; ++i)
        parent[i] = i;
    }
    int find(int i)
    {
      if (parent[i] == i) return i;
      return parent[i] = find(parent[i]);
    }
    bool unite(int i, int j)
    {
      int root_i = find(i);
      int root_j = find(j);
      if (root_i != root_j)
      {
        parent[root_i] = root_j;
        return true;
      }
      return false;
    }
  };

  std::vector<MeshMSTEdge> edge_list;
  edge_list.reserve(candidate_edges.size());

  for (auto &[key, edge] : candidate_edges)
    edge_list.push_back(edge);

  std::sort(edge_list.begin(), edge_list.end());

  MeshDSU                  dsu(total_sources);
  std::vector<MeshMSTEdge> mst_edges;

  for (const auto &edge : edge_list)
  {
    if (dsu.unite(edge.u, edge.v)) mst_edges.push_back(edge);
  }

  // Fallback: If any sink component is not connected to the boundary outlet,
  // connect it directly
  int boundary_root = dsu.find(boundary_src_id);
  for (int s = 0; s < n_sinks; ++s)
  {
    if (dsu.find(s) != boundary_root)
    {
      float       best_cost = std::numeric_limits<float>::max();
      MeshMSTEdge best_edge;
      bool        found = false;

      for (const auto &edge : edge_list)
      {
        if ((dsu.find(edge.u) == dsu.find(s) &&
             dsu.find(edge.v) == boundary_root) ||
            (dsu.find(edge.v) == dsu.find(s) &&
             dsu.find(edge.u) == boundary_root))
        {
          if (edge.cost < best_cost)
          {
            best_cost = edge.cost;
            best_edge = edge;
            found = true;
          }
        }
      }

      if (found && dsu.unite(best_edge.u, best_edge.v))
      {
        mst_edges.push_back(best_edge);
        boundary_root = dsu.find(boundary_src_id);
      }
    }
  }

  // --- Build directed adjacency tree rooted at the boundary outlet

  std::vector<std::vector<std::pair<int, std::vector<size_t>>>> adj(
      total_sources);
  for (const auto &edge : mst_edges)
  {
    adj[edge.u].push_back({edge.v, edge.path});
    std::vector<size_t> rev_path = edge.path;
    std::reverse(rev_path.begin(), rev_path.end());
    adj[edge.v].push_back({edge.u, std::move(rev_path)});
  }

  // BFS from boundary_src_id inward to orient all edges towards boundary
  std::vector<bool> visited(total_sources, false);
  std::vector<int>  bfs_queue;
  bfs_queue.push_back(boundary_src_id);
  visited[boundary_src_id] = true;

  struct DirectedMeshPath
  {
    int                 child;
    int                 parent;
    std::vector<size_t> path; // from upstream child to downstream parent
  };
  std::vector<DirectedMeshPath> directed_paths;

  size_t qhead = 0;
  while (qhead < bfs_queue.size())
  {
    int u = bfs_queue[qhead++];

    for (const auto &[v, path_u_to_v] : adj[u])
    {
      if (!visited[v])
      {
        visited[v] = true;
        bfs_queue.push_back(v);

        std::vector<size_t> path_v_to_u = path_u_to_v;
        std::reverse(path_v_to_u.begin(), path_v_to_u.end());
        directed_paths.push_back({v, u, std::move(path_v_to_u)});
      }
    }
  }

  std::reverse(directed_paths.begin(), directed_paths.end());

  // --- Enforce monotonic downstream elevations and construct Path objects

  std::vector<float> node_z(n_vertices);
  for (size_t i = 0; i < n_vertices; ++i)
    node_z[i] = points[i].z;

  float min_d = std::max(minimum_depth, 0.f);

  for (const auto &dp : directed_paths)
  {
    const auto &path = dp.path;
    if (path.size() < 2) continue;

    node_z[path.front()] = std::min(node_z[path.front()],
                                    points[path.front()].z - min_d);
    float current_z = node_z[path.front()];

    for (size_t idx = 1; idx < path.size(); ++idx)
    {
      size_t curr = path[idx];
      size_t prev = path[idx - 1];
      float  dx = points[curr].x - points[prev].x;
      float  dy = points[curr].y - points[prev].y;
      float  dist = std::hypot(dx, dy);

      current_z -= std::max(riverbed_talus, 1e-6f) * dist;
      float target_z = std::min(current_z, points[curr].z - min_d);

      if (node_z[curr] > target_z) node_z[curr] = target_z;
      current_z = node_z[curr];
    }
  }

  std::vector<Path> output_paths;
  output_paths.reserve(directed_paths.size());

  for (const auto &dp : directed_paths)
  {
    const auto &path_indices = dp.path;
    if (path_indices.size() < 2) continue;

    std::vector<Point> pts;
    pts.reserve(path_indices.size());
    for (size_t v_idx : path_indices)
    {
      pts.push_back(Point(points[v_idx].x, points[v_idx].y, node_z[v_idx]));
    }
    output_paths.emplace_back(std::move(pts));
  }

  return output_paths;
}

std::vector<Path> flow_fixing_mst_paths(const Array  &z,
                                        size_t        control_points_count,
                                        std::uint32_t seed,
                                        float         riverbed_talus,
                                        float         elevation_ratio,
                                        float         distance_exponent,
                                        float         upward_penalization,
                                        float         minimum_depth)
{
  if (!validate_non_empty(z)) return {};

  TerrainTriMesh mesh = generate_terrain_tri_mesh_from_heightmap_random(
      z,
      control_points_count,
      seed);

  return flow_fixing_mst_paths(mesh,
                               riverbed_talus,
                               elevation_ratio,
                               distance_exponent,
                               upward_penalization,
                               minimum_depth);
}

Array flow_fixing_mst_triangulated(const Array  &z,
                                   size_t        control_points_count,
                                   std::uint32_t seed,
                                   float         riverbed_talus,
                                   float         elevation_ratio,
                                   float         distance_exponent,
                                   float         upward_penalization,
                                   float         minimum_depth,
                                   float         merging_distance,
                                   RadialProfile radial_profile,
                                   float         radial_profile_parameter,
                                   const Array  *p_noise_r)
{
  if (!validate_non_empty(z)) return Array();
  if (p_noise_r && !validate_same_shape(z, *p_noise_r)) return Array();

  std::vector<Path> paths = flow_fixing_mst_paths(z,
                                                  control_points_count,
                                                  seed,
                                                  riverbed_talus,
                                                  elevation_ratio,
                                                  distance_exponent,
                                                  upward_penalization,
                                                  minimum_depth);

  Array zb = z;
  if (paths.empty()) return zb;

  float trench_width = merging_distance / float(z.shape.x);

  trench(zb,
         paths,
         trench_width,
         /* enable_width_depth_scaling */ false,
         /* enable_width_distance_scaling */ false,
         /* enable_width_curvature_scaling */ false,
         /* curvature_radius_min */ 1.f,
         /* curv_width_ratio_min */ 0.5f,
         /* curv_width_ratio_max */ 2.f,
         radial_profile,
         radial_profile_parameter,
         ElevationLongitudinalProfile::ELP_DECREASING,
         /* elevation_shift */ 0.f,
         /* shift_ramp_start_ratio */ 0.f,
         /* shift_ramp_end_ratio */ 0.f,
         /* min_slope */ std::max(riverbed_talus, 1e-4f),
         /* k_neighbors */ 4,
         /* p_noise_r */ p_noise_r);

  return zb;
}

} // namespace hmap

namespace hmap::va
{

std::vector<Path> flow_fixing_mst_paths(const ComputeMode  &cm,
                                        const VirtualArray &z,
                                        size_t        control_points_count,
                                        std::uint32_t seed,
                                        float         riverbed_talus,
                                        float         elevation_ratio,
                                        float         distance_exponent,
                                        float         upward_penalization,
                                        float         minimum_depth)
{
  TerrainTriMesh mesh = va::generate_terrain_tri_mesh_from_heightmap_random(
      cm,
      z,
      control_points_count,
      seed);

  return flow_fixing_mst_paths(mesh,
                               riverbed_talus,
                               elevation_ratio,
                               distance_exponent,
                               upward_penalization,
                               minimum_depth);
}

VirtualArray flow_fixing_mst_triangulated(const ComputeMode  &cm,
                                          const VirtualArray &z,
                                          size_t        control_points_count,
                                          std::uint32_t seed,
                                          float         riverbed_talus,
                                          float         elevation_ratio,
                                          float         distance_exponent,
                                          float         upward_penalization,
                                          float         minimum_depth,
                                          float         merging_distance,
                                          RadialProfile radial_profile,
                                          float radial_profile_parameter,
                                          const VirtualArray *p_noise_r)
{
  std::vector<Path> paths = flow_fixing_mst_paths(cm,
                                                  z,
                                                  control_points_count,
                                                  seed,
                                                  riverbed_talus,
                                                  elevation_ratio,
                                                  distance_exponent,
                                                  upward_penalization,
                                                  minimum_depth);

  VirtualArray zb;
  zb.copy_from(z, cm, /* copy_src_data */ true);

  if (paths.empty()) return zb;

  float trench_width = merging_distance / float(z.shape.x);

  hmap::for_each_tile(
      {p_noise_r},
      {&zb},
      [&](std::vector<const hmap::Array *> p_arrays_in,
          std::vector<hmap::Array *>       p_arrays_out,
          const hmap::TileRegion          &region)
      {
        auto [pa_noise_r] = unpack<1>(p_arrays_in);
        auto [pa_zb] = unpack<1>(p_arrays_out);

        trench(*pa_zb,
               paths,
               trench_width,
               /* enable_width_depth_scaling */ false,
               /* enable_width_distance_scaling */ false,
               /* enable_width_curvature_scaling */ false,
               /* curvature_radius_min */ 1.f,
               /* curv_width_ratio_min */ 0.5f,
               /* curv_width_ratio_max */ 2.f,
               radial_profile,
               radial_profile_parameter,
               ElevationLongitudinalProfile::ELP_DECREASING,
               /* elevation_shift */ 0.f,
               /* shift_ramp_start_ratio */ 0.f,
               /* shift_ramp_end_ratio */ 0.f,
               /* min_slope */ std::max(riverbed_talus, 1e-4f),
               /* k_neighbors */ 4,
               /* p_noise_r */ pa_noise_r,
               /* p_bending_mask */ nullptr,
               /* tile_size */ 32,
               /* bbox */ region.bbox);
      },
      cm);

  return zb;
}

} // namespace hmap::va
