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
#include "highmap/math/profiles.hpp"
#include "highmap/shortest_path.hpp"
#include "highmap/terrain_tri_mesh.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <format>
#include <unordered_map>

namespace hmap
{

Array flow_fixing_mst(const Array  &z,
                      float         riverbed_talus,
                      float         elevation_ratio,
                      float         distance_exponent,
                      float         upward_penalization,
                      float         valley_affinity,
                      int           prefilter_ir,
                      float         minimum_depth,
                      bool          carve_riverbed,
                      float         merging_distance,
                      RadialProfile radial_profile,
                      float         radial_profile_parameter,
                      const Array  *p_noise_r,
                      bool          use_midpoint,
                      float         offset_ratio)
{
  if (!validate_non_empty(z)) return Array();
  if (p_noise_r && !validate_same_shape(z, *p_noise_r)) return Array();

  const glm::ivec2 shape = z.shape;
  Array            zb = z;

  const int   di[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
  const int   dj[8] = {0, 1, -1, 0, -1, 1, -1, 1};
  const float cd[8] = {1.f, 1.f, 1.f, 1.f, M_SQRT2, M_SQRT2, M_SQRT2, M_SQRT2};

  auto is_inside = [&shape](int i, int j)
  { return i >= 0 && i < shape.x && j >= 0 && j < shape.y; };

  // --- Graph structures

  // Disjoint Set (Union-Find) structure for Kruskal's MST
  struct DSU
  {
    std::vector<int> parent;
    DSU(int n) : parent(n)
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

  struct MSTEdge
  {
    float                   cost;
    int                     u;
    int                     v;
    std::vector<glm::ivec2> path; // from u to v

    bool operator<(const MSTEdge &other) const
    {
      return cost < other.cost;
    }
  };

  struct DijkstraNode
  {
    float dist;
    int   i;
    int   j;

    bool operator>(const DijkstraNode &o) const
    {
      return dist > o.dist;
    }
  };

  Array zf = zb;
  if (prefilter_ir > 0) smooth_cpulse(zf, prefilter_ir);

  std::vector<glm::ivec2> sinks = find_flow_sinks(zf);
  if (sinks.empty()) return zb;

  // Compute terrain concavity/valley affinity field via Laplace
  Array valley_field(shape, 0.f);
  if (valley_affinity > 0.f)
  {
    valley_field = zf;
    laplace(valley_field); // positive in concave valleys / troughs
    float max_v = valley_field.max();
    float min_v = valley_field.min();
    float span = std::max(max_v - min_v, 1e-6f);
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
        valley_field(i, j) = (valley_field(i, j) - min_v) /
                             span; // in [0, 1], 1 = deepest valley
  }

  // Source 0 ... N_sinks-1: interior sinks and Source N_sinks: virtual
  // outlet node (representing all boundary pixels)
  int n_sinks = static_cast<int>(sinks.size());
  int boundary_src_id = n_sinks;
  int total_sources = n_sinks + 1;

  glm::vec2 z_range = zb.range();
  float     z_span = std::max(z_range.y - z_range.x, 1e-6f);

  // Store candidate edges between meeting sources: key = (min(u,v), max(u,v))
  std::unordered_map<int64_t, MSTEdge> candidate_edges;

  auto make_key = [](int u, int v) -> int64_t
  {
    if (u > v) std::swap(u, v);
    return (static_cast<int64_t>(u) << 32) | static_cast<int64_t>(v);
  };

  if (!use_midpoint)
  {
    // --- Multi-source Dijkstra expansion

    Mat<float>      dist_map(shape, std::numeric_limits<float>::max());
    Mat<int>        owner_map(shape, -1);
    Mat<glm::ivec2> prev_cell(shape, {-1, -1});

    std::priority_queue<DijkstraNode,
                        std::vector<DijkstraNode>,
                        std::greater<DijkstraNode>>
        pq;

    // initialize boundary outlets (Source ID = n_sinks)
    auto init_boundary_cell = [&](int i, int j)
    {
      float norm_z = (zb(i, j) - z_range.x) / z_span;
      float init_d = 5.f *
                     norm_z; // higher border cells have a starting cost penalty
      dist_map(i, j) = init_d;
      owner_map(i, j) = boundary_src_id;
      prev_cell(i, j) = {i, j};
      pq.push({init_d, i, j});
    };

    for (int i = 0; i < shape.x; ++i)
    {
      init_boundary_cell(i, 0);
      init_boundary_cell(i, shape.y - 1);
    }

    for (int j = 1; j < shape.y - 1; ++j)
    {
      init_boundary_cell(0, j);
      init_boundary_cell(shape.x - 1, j);
    }

    // initialize sinks (Source ID = 0 ... n_sinks - 1)
    for (int s = 0; s < n_sinks; ++s)
    {
      glm::ivec2 p = sinks[s];
      dist_map(p) = 0.f;
      owner_map(p) = s;
      prev_cell(p) = p;
      pq.push({0.f, p.x, p.y});
    }

    while (!pq.empty())
    {
      DijkstraNode top = pq.top();
      pq.pop();

      int ci = top.i;
      int cj = top.j;

      if (top.dist > dist_map(ci, cj)) continue;

      int cur_owner = owner_map(ci, cj);

      for (int k = 0; k < 8; ++k)
      {
        int ni = ci + di[k];
        int nj = cj + dj[k];

        if (!is_inside(ni, nj)) continue;

        int nb_owner = owner_map(ni, nj);

        // when meeting a different source region, record candidate bridge edge
        if (nb_owner != -1 && nb_owner != cur_owner)
        {
          float   total_cost = dist_map(ci, cj) + dist_map(ni, nj) + cd[k];
          int64_t key = make_key(cur_owner, nb_owner);

          if (candidate_edges.find(key) == candidate_edges.end() ||
              total_cost < candidate_edges[key].cost)
          {
            // reconstruct path from cur_owner to nb_owner through (ci, cj) -
            // (ni, nj)
            std::vector<glm::ivec2> p1;
            glm::ivec2              curr = {ci, cj};
            while (true)
            {
              p1.push_back(curr);
              glm::ivec2 nxt = prev_cell(curr);
              if (nxt == curr) break;
              curr = nxt;
            }
            std::reverse(p1.begin(), p1.end()); // from source to (ci, cj)

            std::vector<glm::ivec2> p2;
            curr = {ni, nj};
            while (true)
            {
              p2.push_back(curr);
              glm::ivec2 nxt = prev_cell(curr);
              if (nxt == curr) break;
              curr = nxt;
            } // from (ni, nj) to other source

            p1.insert(p1.end(), p2.begin(), p2.end());
            candidate_edges[key] = {total_cost,
                                    cur_owner,
                                    nb_owner,
                                    std::move(p1)};
          }
        }

        // Dijkstra transition cost from (ci, cj) to (ni, nj)
        float dz = (zb(ni, nj) - zb(ci, cj)) * cd[k];
        float cost_step = (1.f - elevation_ratio) * cd[k];

        if (dz > 0.f)
          cost_step += upward_penalization * std::pow(dz, distance_exponent);
        else
          cost_step += std::abs(dz);

        cost_step += elevation_ratio * std::max(0.f, zb(ni, nj));

        // valley / concavity affinity: reduce cost in natural valleys and
        // depressions
        if (valley_affinity > 0.f)
        {
          float val_factor = 1.f - valley_affinity * valley_field(ni, nj);
          cost_step *= std::max(0.1f, val_factor);
        }

        float new_dist = dist_map(ci, cj) + cost_step;

        if (new_dist < dist_map(ni, nj))
        {
          dist_map(ni, nj) = new_dist;
          owner_map(ni, nj) = cur_owner;
          prev_cell(ni, nj) = {ci, cj};
          pq.push({new_dist, ni, nj});
        }
      }
    }
  }
  else
  {
    // --- Midpoint displacement pathfinding

    auto compute_path_cost = [&](const std::vector<glm::ivec2> &path) -> float
    {
      float cost = 0.f;
      for (size_t k = 0; k + 1 < path.size(); ++k)
      {
        glm::ivec2 curr = path[k];
        glm::ivec2 nxt = path[k + 1];
        float      dx = float(nxt.x - curr.x);
        float      dy = float(nxt.y - curr.y);
        float      step_len = std::hypot(dx, dy);
        float      dz = (zb(nxt) - zb(curr));
        float      cost_step = (1.f - elevation_ratio) * step_len;

        if (dz > 0.f)
          cost_step += upward_penalization * std::pow(dz, distance_exponent);
        else
          cost_step += std::abs(dz);

        cost_step += elevation_ratio * std::max(0.f, zb(nxt));

        if (valley_affinity > 0.f)
        {
          float val_factor = 1.f - valley_affinity * valley_field(nxt);
          cost_step *= std::max(0.1f, val_factor);
        }

        cost += cost_step;
      }
      return cost;
    };

    // candidate paths between k-nearest sink pairs
    const int max_neighbors = std::min(n_sinks - 1, 8);
    for (int u = 0; u < n_sinks; ++u)
    {
      std::vector<std::pair<float, int>> neighbors;
      neighbors.reserve(n_sinks - 1);
      for (int v = 0; v < n_sinks; ++v)
      {
        if (u == v) continue;
        float d2 = float((sinks[u].x - sinks[v].x) * (sinks[u].x - sinks[v].x) +
                         (sinks[u].y - sinks[v].y) * (sinks[u].y - sinks[v].y));
        neighbors.push_back({d2, v});
      }
      std::partial_sort(neighbors.begin(),
                        neighbors.begin() + max_neighbors,
                        neighbors.end());

      for (int i = 0; i < max_neighbors; ++i)
      {
        int v = neighbors[i].second;
        if (u >= v) continue;
        int64_t key = make_key(u, v);
        if (candidate_edges.find(key) != candidate_edges.end()) continue;

        std::vector<glm::ivec2> path = find_path_midpoint(zb,
                                                          sinks[u],
                                                          sinks[v],
                                                          offset_ratio,
                                                          0,  // max_it
                                                          4); // steps
        if (!path.empty())
        {
          float cost = compute_path_cost(path);
          candidate_edges[key] = {cost, u, v, std::move(path)};
        }
      }
    }

    // find lowest boundary cells on each edge
    glm::ivec2 min_b = {0, 0};
    glm::ivec2 min_t = {0, shape.y - 1};
    glm::ivec2 min_l = {0, 0};
    glm::ivec2 min_r = {shape.x - 1, 0};

    for (int i = 0; i < shape.x; ++i)
    {
      if (zb(i, 0) < zb(min_b)) min_b = {i, 0};
      if (zb(i, shape.y - 1) < zb(min_t)) min_t = {i, shape.y - 1};
    }
    for (int j = 0; j < shape.y; ++j)
    {
      if (zb(0, j) < zb(min_l)) min_l = {0, j};
      if (zb(shape.x - 1, j) < zb(min_r)) min_r = {shape.x - 1, j};
    }

    // candidate paths from each sink to boundary
    for (int u = 0; u < n_sinks; ++u)
    {
      glm::ivec2 p = sinks[u];

      // candidate boundary targets: orthogonal projections + lowest border
      // points
      std::vector<glm::ivec2> b_candidates = {{p.x, 0},
                                              {p.x, shape.y - 1},
                                              {0, p.y},
                                              {shape.x - 1, p.y},
                                              min_b,
                                              min_t,
                                              min_l,
                                              min_r};

      float                   best_cost = std::numeric_limits<float>::max();
      std::vector<glm::ivec2> best_path;

      for (const auto &b_pt : b_candidates)
      {
        std::vector<glm::ivec2> path = find_path_midpoint(zb,
                                                          p,
                                                          b_pt,
                                                          offset_ratio,
                                                          0,  // max_it
                                                          4); // steps
        if (path.empty()) continue;

        float norm_z = (zb(b_pt) - z_range.x) / z_span;
        float total_cost = compute_path_cost(path) + 5.f * norm_z;

        if (total_cost < best_cost)
        {
          best_cost = total_cost;
          best_path = std::move(path);
        }
      }

      if (!best_path.empty())
      {
        int64_t key = make_key(u, boundary_src_id);
        candidate_edges[key] = {best_cost,
                                u,
                                boundary_src_id,
                                std::move(best_path)};
      }
    }
  }

  // --- Build Kruskal Minimum Spanning Tree across all sinks + boundary outlet

  std::vector<MSTEdge> edge_list;
  edge_list.reserve(candidate_edges.size());

  for (auto &[key, edge] : candidate_edges)
    edge_list.push_back(edge);

  std::sort(edge_list.begin(), edge_list.end());

  DSU                  dsu(total_sources);
  std::vector<MSTEdge> mst_edges;

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
      // Find best edge connecting sink component s to boundary
      float   best_cost = std::numeric_limits<float>::max();
      MSTEdge best_edge;
      bool    found = false;

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

  std::vector<std::vector<std::pair<int, std::vector<glm::ivec2>>>> adj(
      total_sources);
  for (const auto &edge : mst_edges)
  {
    adj[edge.u].push_back({edge.v, edge.path});
    // Reverse path for other direction
    std::vector<glm::ivec2> rev_path = edge.path;
    std::reverse(rev_path.begin(), rev_path.end());
    adj[edge.v].push_back({edge.u, std::move(rev_path)});
  }

  // BFS from boundary_src_id inward to orient all edges towards the boundary
  std::vector<bool> visited(total_sources, false);
  std::vector<int>  bfs_queue;
  bfs_queue.push_back(boundary_src_id);
  visited[boundary_src_id] = true;

  struct DirectedPath
  {
    int                     child;
    int                     parent;
    std::vector<glm::ivec2> path; // from upstream child to downstream parent
  };
  std::vector<DirectedPath> directed_paths;

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

        // Path from child v to parent u (towards boundary)
        std::vector<glm::ivec2> path_v_to_u = path_u_to_v;
        std::reverse(path_v_to_u.begin(), path_v_to_u.end());
        directed_paths.push_back({v, u, std::move(path_v_to_u)});
      }
    }
  }

  std::reverse(directed_paths.begin(), directed_paths.end());

  for (const auto &dp : directed_paths)
  {
    const auto &path = dp.path;
    if (path.size() < 2) continue;

    float min_d = std::max(minimum_depth, 0.f);
    zb(path.front()) = std::min(zb(path.front()), z(path.front()) - min_d);
    float current_z = zb(path.front());

    for (size_t idx = 1; idx < path.size(); ++idx)
    {
      glm::ivec2 curr = path[idx];
      glm::ivec2 prev = path[idx - 1];
      int        dx = curr.x - prev.x;
      int        dy = curr.y - prev.y;

      // decrement elevation along downstream flow direction
      float dist = std::hypot(float(dx), float(dy));
      current_z -= std::max(riverbed_talus, 1e-6f) * dist;

      float target_z = std::min(current_z, z(curr) - min_d);

      if (zb(curr) > target_z) zb(curr) = target_z;

      current_z = zb(curr);
    }
  }

  // --- Optional riverbed carving and smoothing with trench

  if (carve_riverbed)
  {
    float trench_width = merging_distance / float(shape.x);

    std::vector<Path> river_paths;
    river_paths.reserve(directed_paths.size());

    for (const auto &dp : directed_paths)
    {
      const auto &path_cells = dp.path;
      if (path_cells.size() < 2) continue;

      std::vector<Point> pts;
      pts.reserve(path_cells.size());
      for (const auto &p : path_cells)
      {
        float x = (float(p.x) + 0.5f) / float(shape.x);
        float y = (float(p.y) + 0.5f) / float(shape.y);
        pts.push_back(Point(x, y, zb(p)));
      }

      river_paths.push_back(Path(pts));
    }

    trench(zb,
           river_paths,
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
  }

  return zb;
}

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

  const glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f};
  Cloud           cloud = random_cloud_jittered(control_points_count,
                                                {0.5f, 0.5f},
                                                {0.f, 0.f},
                                      seed,
                                      bbox);
  cloud.snap_points_to_bounding_box(bbox);
  cloud.set_values_from_array(z, bbox);
  auto mesh = TerrainTriMesh(cloud.to_vec3());

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
  const glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f};
  Cloud           cloud = random_cloud_jittered(control_points_count,
                                                {0.5f, 0.5f},
                                                {0.f, 0.f},
                                      seed,
                                      bbox);
  cloud.snap_points_to_bounding_box(bbox);

  // sample the control points tile after tile in a fixed order: all tiles
  // write into the same cloud, so distributed sampling races on the points in
  // the overlaps and the resulting sinks change from one run to the next
  ComputeMode cm_sampling = cm;
  cm_sampling.mode = ForEachMode::VA_SEQUENTIAL;

  hmap::for_each_tile(
      {&z},
      {},
      [&](std::vector<const hmap::Array *> p_arrays_in,
          std::vector<hmap::Array *>,
          const hmap::TileRegion &region)
      {
        auto [pa_z] = unpack<1>(p_arrays_in);
        cloud.set_values_from_array(*pa_z, region.bbox);
      },
      cm_sampling);

  auto mesh = TerrainTriMesh(cloud.to_vec3());

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
               /* bbox */ region.bbox);
      },
      cm);

  return zb;
}

} // namespace hmap::va
