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
#include "highmap/virtual_array/communicator.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <format>
#include <unordered_map>

namespace hmap
{

Array flow_fixing_mst(const Array      &z,
                      float             riverbed_talus,
                      float             elevation_ratio,
                      float             distance_exponent,
                      float             upward_penalization,
                      float             valley_affinity,
                      int               prefilter_ir,
                      float             minimum_depth,
                      bool              carve_riverbed,
                      float             merging_distance,
                      RadialProfile     radial_profile,
                      float             radial_profile_parameter,
                      const Array      *p_noise_r,
                      bool              use_midpoint,
                      float             offset_ratio,
                      const TileRegion &region,
                      ICommunicator    &comm)
{
  if (!validate_non_empty(z)) return Array();
  if (p_noise_r && !validate_same_shape(z, *p_noise_r)) return Array();

  // --- Distributed tile-communicated execution

  if (comm.size() > 1)
  {
    Array zb = z;

    const float z_min = comm.all_reduce(z.min(), ReduceOp::Min);
    const float z_max = comm.all_reduce(z.max(), ReduceOp::Max);
    const float z_span = std::max(z_max - z_min, 1e-6f);

    Array zf = z;
    if (prefilter_ir > 0) smooth_cpulse(zf, prefilter_ir);

    std::vector<glm::ivec2> sinks = find_flow_sinks(zf);

    Array valley_field(z.shape, 0.f);
    if (valley_affinity > 0.f)
    {
      valley_field = zf;
      laplace(valley_field);
      float max_v = valley_field.max();
      float min_v = valley_field.min();
      float span = std::max(max_v - min_v, 1e-6f);
      for (int j = 0; j < z.shape.y; ++j)
        for (int i = 0; i < z.shape.x; ++i)
          valley_field(i, j) = (valley_field(i, j) - min_v) / span;
    }

    auto to_global = [&region](int i, int j) -> glm::ivec2
    {
      int gi = region.key.tx * region.tile_shape.x - region.halo.x + i;
      int gj = region.key.ty * region.tile_shape.y - region.halo.z + j;
      return {gi, gj};
    };

    auto pack_global_id = [&region](int gi, int gj) -> int64_t
    {
      return static_cast<int64_t>(gj) *
                 static_cast<int64_t>(region.global_shape.x) +
             gi;
    };

    auto compute_path_cost = [&](const std::vector<glm::ivec2> &path) -> float
    {
      if (path.size() < 2) return std::numeric_limits<float>::max();

      float total_cost = 0.f;
      float z_start = z(path.front());

      for (size_t k = 1; k < path.size(); ++k)
      {
        glm::ivec2 p0 = path[k - 1];
        glm::ivec2 p1 = path[k];
        float      step_len = glm::length(glm::vec2(p1 - p0));
        float      dz = z(p1) - z_start;
        float      up_pen = (dz > 0.f) ? upward_penalization * dz : 0.f;
        float      elev_cost = elevation_ratio * ((z(p1) - z_min) / z_span);
        float      valley_bonus = valley_affinity * valley_field(p1);
        float step_cost = std::pow(step_len, distance_exponent) + elev_cost +
                          up_pen - valley_bonus;
        total_cost += std::max(step_cost, 1e-5f);
      }

      return total_cost;
    };

    // serialize local edges: [cost, u_id, v_id, z_start, z_end, n_pts, x0, y0,
    // ...]
    std::vector<float> local_packed;
    auto               push_edge = [&](float                          cost,
                         int64_t                        u_id,
                         int64_t                        v_id,
                         const std::vector<glm::ivec2> &path)
    {
      local_packed.push_back(cost);
      local_packed.push_back(static_cast<float>(u_id));
      local_packed.push_back(static_cast<float>(v_id));
      local_packed.push_back(z(path.front()));
      local_packed.push_back(z(path.back()));
      local_packed.push_back(static_cast<float>(path.size()));
      for (const auto &p : path)
      {
        glm::vec2 c = region.cell_center(p.x, p.y);
        local_packed.push_back(c.x);
        local_packed.push_back(c.y);
      }
    };

    const int n_sinks = static_cast<int>(sinks.size());
    const int k_max = std::min(6, std::max(0, n_sinks - 1));
    const int steps = use_midpoint ? 4 : 6;

    for (int u = 0; u < n_sinks; ++u)
    {
      std::vector<std::pair<float, int>> neighbors;
      for (int v = 0; v < n_sinks; ++v)
      {
        if (u == v) continue;
        float d2 = float((sinks[u].x - sinks[v].x) * (sinks[u].x - sinks[v].x) +
                         (sinks[u].y - sinks[v].y) * (sinks[u].y - sinks[v].y));
        neighbors.push_back({d2, v});
      }

      if (k_max > 0)
      {
        std::partial_sort(neighbors.begin(),
                          neighbors.begin() + k_max,
                          neighbors.end());
        for (int i = 0; i < k_max; ++i)
        {
          int v = neighbors[i].second;
          if (u >= v) continue;
          auto path = find_path_midpoint(z,
                                         sinks[u],
                                         sinks[v],
                                         offset_ratio,
                                         0,
                                         steps);
          if (path.size() >= 2)
          {
            glm::ivec2 gu = to_global(sinks[u].x, sinks[u].y);
            glm::ivec2 gv = to_global(sinks[v].x, sinks[v].y);
            push_edge(compute_path_cost(path),
                      pack_global_id(gu.x, gu.y),
                      pack_global_id(gv.x, gv.y),
                      path);
          }
        }
      }

      glm::ivec2              p = sinks[u];
      glm::ivec2              gu = to_global(p.x, p.y);
      int64_t                 u_id = pack_global_id(gu.x, gu.y);
      std::vector<glm::ivec2> b_candidates = {{p.x, 0},
                                              {p.x, z.shape.y - 1},
                                              {0, p.y},
                                              {z.shape.x - 1, p.y}};

      for (const auto &b_pt : b_candidates)
      {
        auto path = find_path_midpoint(z, p, b_pt, offset_ratio, 0, steps);
        if (path.size() < 2) continue;

        glm::ivec2 gb = to_global(b_pt.x, b_pt.y);
        bool       on_global_boundary = (gb.x <= 0 ||
                                   gb.x >= region.global_shape.x - 1 ||
                                   gb.y <= 0 ||
                                   gb.y >= region.global_shape.y - 1);
        int64_t    v_id = on_global_boundary ? -1 : pack_global_id(gb.x, gb.y);
        float      norm_z = (z(b_pt) - z_min) / z_span;
        push_edge(compute_path_cost(path) + 5.f * norm_z, u_id, v_id, path);
      }
    }

    std::vector<float> all_packed = comm.all_union(local_packed);

    struct DistEdge
    {
      float                  cost;
      int64_t                u;
      int64_t                v;
      float                  z_start;
      float                  z_end;
      std::vector<glm::vec2> coords;
      bool                   operator<(const DistEdge &o) const
      {
        return cost < o.cost;
      }
    };

    std::vector<DistEdge> edges;
    size_t                idx = 0;
    while (idx + 6 <= all_packed.size())
    {
      DistEdge e;
      e.cost = all_packed[idx++];
      e.u = static_cast<int64_t>(all_packed[idx++]);
      e.v = static_cast<int64_t>(all_packed[idx++]);
      e.z_start = all_packed[idx++];
      e.z_end = all_packed[idx++];
      size_t n_pts = static_cast<size_t>(all_packed[idx++]);
      if (idx + 2 * n_pts > all_packed.size()) break;
      e.coords.reserve(n_pts);
      for (size_t k = 0; k < n_pts; ++k)
      {
        float cx = all_packed[idx++];
        float cy = all_packed[idx++];
        e.coords.emplace_back(cx, cy);
      }
      edges.push_back(std::move(e));
    }

    std::sort(edges.begin(), edges.end());

    std::unordered_map<int64_t, int64_t> parent;
    auto find_root = [&parent](int64_t x) -> int64_t
    {
      int64_t root = x;
      while (parent.count(root) && parent[root] != root)
        root = parent[root];
      return root;
    };

    std::vector<Path> mst_paths;
    for (const auto &e : edges)
    {
      if (!parent.count(e.u)) parent[e.u] = e.u;
      if (!parent.count(e.v)) parent[e.v] = e.v;
      int64_t ru = find_root(e.u);
      int64_t rv = find_root(e.v);
      if (ru != rv)
      {
        parent[ru] = rv;
        std::vector<Point> pts;
        pts.reserve(e.coords.size());
        for (size_t k = 0; k < e.coords.size(); ++k)
        {
          float t = float(k) / float(std::max<size_t>(1, e.coords.size() - 1));
          pts.emplace_back(e.coords[k].x,
                           e.coords[k].y,
                           glm::mix(e.z_start, e.z_end, t) - minimum_depth);
        }
        mst_paths.emplace_back(pts);
      }
    }

    if (carve_riverbed && !mst_paths.empty())
    {
      float trench_width = merging_distance / float(region.global_shape.x);
      trench(zb,
             mst_paths,
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
             /* p_noise_r */ p_noise_r,
             /* p_bending_mask */ nullptr,
             /* tile_size */ 32,
             /* bbox */ region.bbox);
    }

    comm.sync_overlap_buffers(zb, region, SyncOperation::Min);
    return zb;
  }

  const glm::ivec2 shape = z.shape;
  Array            zb = z;

  const int   di[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
  const int   dj[8] = {0, 1, -1, 0, -1, 1, -1, 1};
  const float cd[8] =
      {1.f, 1.f, 1.f, 1.f, cte::sqrt2, cte::sqrt2, cte::sqrt2, cte::sqrt2};

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

} // namespace hmap

namespace hmap::va
{

VirtualArray flow_fixing_mst(const ComputeMode  &cm,
                             const VirtualArray &z,
                             float               riverbed_talus,
                             float               elevation_ratio,
                             float               distance_exponent,
                             float               upward_penalization,
                             float               valley_affinity,
                             int                 prefilter_ir,
                             float               minimum_depth,
                             bool                carve_riverbed,
                             float               merging_distance,
                             RadialProfile       radial_profile,
                             float               radial_profile_parameter,
                             const VirtualArray *p_noise_r,
                             bool                use_midpoint,
                             float               offset_ratio)
{
  VirtualArray zb;
  zb.copy_from(z, cm, /* copy_src_data */ true);

  if (z.shape.x == 0 || z.shape.y == 0) return zb;

  // --- Global elevation range

  glm::vec2 z_range = {z.min(cm), z.max(cm)};
  float     z_span = std::max(z_range.y - z_range.x, 1e-6f);

  // --- Generic path-finder abstraction (midpoint by default)

  using PathFinder = std::function<
      std::vector<glm::ivec2>(const Array &, glm::ivec2, glm::ivec2)>;

  PathFinder path_finder;
  if (use_midpoint)
  {
    path_finder =
        [offset_ratio](const Array &arr, glm::ivec2 start, glm::ivec2 end)
    {
      return find_path_midpoint(arr,
                                start,
                                end,
                                offset_ratio,
                                0,  // max_it
                                4); // steps
    };
  }
  else
  {
    path_finder =
        [offset_ratio](const Array &arr, glm::ivec2 start, glm::ivec2 end)
    {
      // fallback to midpoint with finer steps when dijkstra is not used
      return find_path_midpoint(arr,
                                start,
                                end,
                                offset_ratio,
                                0,  // max_it
                                6); // steps
    };
  }

  // --- Graph structures for global MST across tiles

  struct GlobalEdge
  {
    float                  cost;
    int64_t                u;
    int64_t                v;
    std::vector<glm::vec2> path_coords;
    float                  z_start;
    float                  z_end;

    bool operator<(const GlobalEdge &other) const
    {
      return cost < other.cost;
    }
  };

  auto pack_global_id = [&z](int gi, int gj) -> int64_t
  { return static_cast<int64_t>(gj) * static_cast<int64_t>(z.shape.x) + gi; };

  const int64_t boundary_node_id = -1;

  std::mutex              edges_mutex;
  std::vector<GlobalEdge> all_edges;

  // --- Per-tile sink detection & local portal graph extraction

  hmap::for_each_tile(
      {&z},
      {},
      [&](std::vector<const hmap::Array *> p_arrays_in,
          std::vector<hmap::Array *>,
          const hmap::TileRegion &region)
      {
        auto [pa_z] = unpack<1>(p_arrays_in);
        const glm::ivec2 tshape = pa_z->shape;

        Array zf = *pa_z;
        if (prefilter_ir > 0) smooth_cpulse(zf, prefilter_ir);

        std::vector<glm::ivec2> sinks = find_flow_sinks(zf);
        if (sinks.empty()) return;

        // valley affinity field
        Array valley_field(tshape, 0.f);
        if (valley_affinity > 0.f)
        {
          valley_field = zf;
          laplace(valley_field);
          float max_v = valley_field.max();
          float min_v = valley_field.min();
          float span = std::max(max_v - min_v, 1e-6f);
          for (int j = 0; j < tshape.y; ++j)
            for (int i = 0; i < tshape.x; ++i)
              valley_field(i, j) = (valley_field(i, j) - min_v) / span;
        }

        auto to_global = [&region](int i, int j) -> glm::ivec2
        {
          int gi = region.key.tx * region.tile_shape.x - region.halo.x + i;
          int gj = region.key.ty * region.tile_shape.y - region.halo.z + j;
          return {gi, gj};
        };

        auto compute_path_cost =
            [&](const std::vector<glm::ivec2> &path) -> float
        {
          if (path.size() < 2) return std::numeric_limits<float>::max();

          float total_cost = 0.f;
          float z_start = (*pa_z)(path.front());

          for (size_t k = 1; k < path.size(); ++k)
          {
            glm::ivec2 p0 = path[k - 1];
            glm::ivec2 p1 = path[k];
            float      step_len = glm::length(glm::vec2(p1 - p0));
            float      dz = (*pa_z)(p1)-z_start;
            float      up_pen = (dz > 0.f) ? upward_penalization * dz : 0.f;
            float      elev_cost = elevation_ratio *
                              (((*pa_z)(p1)-z_range.x) / z_span);
            float valley_bonus = valley_affinity * valley_field(p1);
            float step_cost = std::pow(step_len, distance_exponent) +
                              elev_cost + up_pen - valley_bonus;
            total_cost += std::max(step_cost, 1e-5f);
          }

          return total_cost;
        };

        auto to_coords = [&region](const std::vector<glm::ivec2> &path)
            -> std::vector<glm::vec2>
        {
          std::vector<glm::vec2> coords;
          coords.reserve(path.size());
          for (const auto &p : path)
            coords.push_back(region.cell_center(p.x, p.y));
          return coords;
        };

        std::vector<GlobalEdge> local_edges;

        // connect each sink to k nearest local sinks
        const int n_sinks = static_cast<int>(sinks.size());
        const int k_max = std::min(6, n_sinks - 1);

        for (int u = 0; u < n_sinks; ++u)
        {
          std::vector<std::pair<float, int>> neighbors;
          neighbors.reserve(n_sinks - 1);
          for (int v = 0; v < n_sinks; ++v)
          {
            if (u == v) continue;
            float d2 = float(
                (sinks[u].x - sinks[v].x) * (sinks[u].x - sinks[v].x) +
                (sinks[u].y - sinks[v].y) * (sinks[u].y - sinks[v].y));
            neighbors.push_back({d2, v});
          }

          if (k_max > 0)
          {
            std::partial_sort(neighbors.begin(),
                              neighbors.begin() + k_max,
                              neighbors.end());
            for (int i = 0; i < k_max; ++i)
            {
              int v = neighbors[i].second;
              if (u >= v) continue;

              std::vector<glm::ivec2> path = path_finder(*pa_z,
                                                         sinks[u],
                                                         sinks[v]);
              if (path.size() >= 2)
              {
                glm::ivec2 gu = to_global(sinks[u].x, sinks[u].y);
                glm::ivec2 gv = to_global(sinks[v].x, sinks[v].y);
                local_edges.push_back({compute_path_cost(path),
                                       pack_global_id(gu.x, gu.y),
                                       pack_global_id(gv.x, gv.y),
                                       to_coords(path),
                                       (*pa_z)(sinks[u]),
                                       (*pa_z)(sinks[v])});
              }
            }
          }
        }

        // find lowest boundary cells on each tile edge
        glm::ivec2 min_b = {0, 0};
        glm::ivec2 min_t = {0, tshape.y - 1};
        glm::ivec2 min_l = {0, 0};
        glm::ivec2 min_r = {tshape.x - 1, 0};

        for (int i = 0; i < tshape.x; ++i)
        {
          if ((*pa_z)(i, 0) < (*pa_z)(min_b)) min_b = {i, 0};
          if ((*pa_z)(i, tshape.y - 1) < (*pa_z)(min_t))
            min_t = {i, tshape.y - 1};
        }
        for (int j = 0; j < tshape.y; ++j)
        {
          if ((*pa_z)(0, j) < (*pa_z)(min_l)) min_l = {0, j};
          if ((*pa_z)(tshape.x - 1, j) < (*pa_z)(min_r))
            min_r = {tshape.x - 1, j};
        }

        // connect each sink to tile boundary portals (or global outlet)
        for (int u = 0; u < n_sinks; ++u)
        {
          glm::ivec2 p = sinks[u];
          glm::ivec2 gu = to_global(p.x, p.y);
          int64_t    u_id = pack_global_id(gu.x, gu.y);

          std::vector<glm::ivec2> b_candidates = {{p.x, 0},
                                                  {p.x, tshape.y - 1},
                                                  {0, p.y},
                                                  {tshape.x - 1, p.y},
                                                  min_b,
                                                  min_t,
                                                  min_l,
                                                  min_r};

          for (const auto &b_pt : b_candidates)
          {
            std::vector<glm::ivec2> path = path_finder(*pa_z, p, b_pt);
            if (path.size() < 2) continue;

            glm::ivec2 gb = to_global(b_pt.x, b_pt.y);
            bool on_global_boundary = (gb.x <= 0 || gb.x >= z.shape.x - 1 ||
                                       gb.y <= 0 || gb.y >= z.shape.y - 1);

            int64_t v_id = on_global_boundary ? boundary_node_id
                                              : pack_global_id(gb.x, gb.y);

            float norm_z = ((*pa_z)(b_pt)-z_range.x) / z_span;
            float cost = compute_path_cost(path) + 5.f * norm_z;

            local_edges.push_back(
                {cost, u_id, v_id, to_coords(path), (*pa_z)(p), (*pa_z)(b_pt)});
          }
        }

        if (!local_edges.empty())
        {
          std::lock_guard<std::mutex> lock(edges_mutex);
          all_edges.insert(all_edges.end(),
                           std::make_move_iterator(local_edges.begin()),
                           std::make_move_iterator(local_edges.end()));
        }
      },
      cm);

  if (all_edges.empty()) return zb;

  // --- Build global Kruskal MST across all sinks, portals & boundary outlet

  std::sort(all_edges.begin(), all_edges.end());

  std::unordered_map<int64_t, int64_t> parent;
  auto find_root = [&parent](int64_t x) -> int64_t
  {
    int64_t root = x;
    while (parent.count(root) && parent[root] != root)
      root = parent[root];
    int64_t curr = x;
    while (curr != root)
    {
      int64_t nxt = parent[curr];
      parent[curr] = root;
      curr = nxt;
    }
    return root;
  };

  std::vector<Path> mst_paths;
  for (const auto &edge : all_edges)
  {
    if (!parent.count(edge.u)) parent[edge.u] = edge.u;
    if (!parent.count(edge.v)) parent[edge.v] = edge.v;

    int64_t ru = find_root(edge.u);
    int64_t rv = find_root(edge.v);
    if (ru != rv)
    {
      parent[ru] = rv;

      std::vector<float> z_vals(edge.path_coords.size());
      for (size_t k = 0; k < edge.path_coords.size(); ++k)
      {
        float t = float(k) /
                  float(std::max<size_t>(1, edge.path_coords.size() - 1));
        z_vals[k] = glm::mix(edge.z_start, edge.z_end, t) - minimum_depth;
      }

      std::vector<Point> pts;
      pts.reserve(edge.path_coords.size());
      for (size_t k = 0; k < edge.path_coords.size(); ++k)
        pts.emplace_back(edge.path_coords[k].x,
                         edge.path_coords[k].y,
                         z_vals[k]);

      mst_paths.emplace_back(pts);
    }
  }

  if (mst_paths.empty() || !carve_riverbed) return zb;

  // --- Per-tile riverbed carving

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
               mst_paths,
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
