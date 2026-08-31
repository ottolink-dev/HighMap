/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <queue>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/carving.hpp"
#include "highmap/filters.hpp"
#include "highmap/hydrology/drainage_basin_cell_based.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/morphology.hpp"
#include "highmap/random.hpp"
#include "highmap/transform.hpp"

#include <unordered_map>

namespace hmap
{

Array flow_fixing(const Array &z,
                  float        riverbed_talus,
                  int          iterations,
                  int          prefilter_ir,
                  bool         carve_riverbed,
                  float        merging_distance,
                  const Array *p_noise_r)
{
  if (!validate_non_empty(z)) return Array();
  if (p_noise_r && !validate_same_shape(z, *p_noise_r)) return Array();

  // local node type for heap queues
  struct Node
  {
    float h;
    int   i;
    int   j;

    bool operator<(const Node &other) const
    {
      return h > other.h; // min-heap
    }
  };

  struct NodePath
  {
    float      h;  // elevation at end point
    glm::ivec2 p0; // start
    glm::ivec2 p1; // end

    bool operator<(const NodePath &other) const
    {
      return h < other.h; // max-heap
    }
  };

  //
  const glm::ivec2 shape = z.shape;
  Array            zb = z;
  size_t           n_sinks = 0;

  // neighbor search
  const int di[8] = {1, 1, 0, -1, -1, -1, 0, 1};
  const int dj[8] = {0, 1, 1, 1, 0, -1, -1, -1};

  auto is_inside = [&shape](int i, int j)
  { return i >= 0 && i < shape.x && j >= 0 && j < shape.y; };

  std::unordered_map<glm::ivec4, std::vector<glm::ivec2>, IVec4Hash, IVec4Eq>
      breach_history;

  // --- main loop

  for (int it = 0; it < iterations; ++it)
  {
    Array zf = zb;
    smooth_cpulse(zf, prefilter_ir);

    std::vector<glm::ivec2> sinks = find_flow_sinks(zf);
    Mat<int>                is_sink(shape, 0);

    for (const auto &p : sinks)
      is_sink(p) = 1;

    if (sinks.size() == n_sinks)
      break;
    else
      n_sinks = sinks.size();

    // --- flow breaching: 1st pass

    std::vector<Node> queue;
    queue.reserve(shape.x * shape.y);

    Mat<int>        visited(shape, 0);
    Mat<glm::ivec2> flow_map(shape, {0, 0});

    // --- initialize heap queue with the lowest cell of each border

    for (int i : {0, shape.x - 1})
    {
      float vmin = 1e30f;
      int   jmin = 0;
      for (int j = 0; j < shape.y; ++j)
      {
        if (zb(i, j) < vmin)
        {
          vmin = zb(i, j);
          jmin = j;
        }
      }
      queue.push_back({vmin, i, jmin});
    }

    for (int j : {0, shape.y - 1})
    {
      float vmin = 1e30f;
      int   imin = 0;
      for (int i = 0; i < shape.x; ++i)
      {
        if (zb(i, j) < vmin)
        {
          vmin = zb(i, j);
          imin = i;
        }
      }
      queue.push_back({vmin, imin, j});
    }

    // --- traverse the queue

    std::make_heap(queue.begin(), queue.end());

    while (!queue.empty())
    {
      std::pop_heap(queue.begin(), queue.end());
      const Node c = queue.back();
      queue.pop_back();

      for (int k = 0; k < 8; ++k)
      {
        int ni = c.i + di[k];
        int nj = c.j + dj[k];

        if (is_inside(ni, nj) && visited(ni, nj) == 0)
        {
          // store flow direction
          flow_map(ni, nj) = {c.i, c.j};
          visited(ni, nj) = 1;
          // queue.push_back({zb(ni, nj), ni, nj});
          queue.push_back({zb(ni, nj) + 1.f, ni, nj});
          // queue.push_back({std::abs(zb(ni, nj) - zb(c.i, c.j)) + c.h, ni,
          // nj});
          std::push_heap(queue.begin(), queue.end());

          // if the current cell is a sink, "breach" the
          // heightmap by following the reverse flow direction
          // in order to connect this sink to another sink, or
          // to connect this connect to the domain border
          if (is_sink(ni, nj))
          {
            int bi = ni;
            int bj = nj;

            bool                    keep_breaching = true;
            std::vector<glm::ivec2> path = {{bi, bj}};

            // stop on a boundary or at a sink
            while (keep_breaching &&
                   (bi > 0 && bi < shape.x - 1 && bj > 0 && bj < shape.y - 1))
            {
              glm::ivec2 tmp = flow_map(bi, bj);
              bi = tmp.x;
              bj = tmp.y;
              path.push_back({bi, bj});
              if (is_sink(bi, bj)) keep_breaching = false;
            }

            // after the breaching path has been identifier,
            // follow the this path and make sure the
            // elevation is monotonic along this path
            if (path.size() > 2)
            {
              glm::ivec2 p0 = path.front();
              glm::ivec2 p1 = path.back();
              if (p0 != p1)
              {
                // store the breaching path for the second
                // pass of the algorithm
                glm::ivec4 key = {p0.x, p0.y, p1.x, p1.y};
                breach_history[key] = path;

                for (size_t r = 0; r < path.size() - 1; ++r)
                {
                  if (zb(path[r + 1]) > zb(path[r]))
                    zb(path[r + 1]) = zb(path[r]) - riverbed_talus;
                  else
                    zb(path[r + 1]) -= riverbed_talus;
                }
              }
            }
          }
        } // if visited
      } // neighbors k-loop
    } // queue

    // --- 2nd pass

    // breach again from top to bottom (hence the -z in the cost) to
    // ensure overall elevations are coherent between the sinks
    std::vector<NodePath> queue_path;
    queue.reserve(breach_history.size());

    for (const auto &[key, path] : breach_history)
    {
      glm::ivec2 p0 = {key.x, key.y};
      glm::ivec2 p1 = {key.z, key.w};
      queue_path.push_back({z(p1), p0, p1});
    }
    std::make_heap(queue_path.begin(), queue_path.end());

    while (!queue.empty())
    {
      std::pop_heap(queue_path.begin(), queue_path.end());
      const NodePath current = queue_path.back();
      queue_path.pop_back();

      glm::ivec4 key = {current.p0.x, current.p0.y, current.p1.x, current.p1.y};
      const std::vector<glm::ivec2> &path = breach_history[key];

      // breach again
      for (size_t r = 0; r < path.size() - 1; ++r)
      {
        if (zb(path[r + 1]) > zb(path[r]))
          zb(path[r + 1]) = zb(path[r]) - riverbed_talus;
        else
          zb(path[r + 1]) -= riverbed_talus;
      }
    }

  } // main it

  // --- carve the river

  if (carve_riverbed)
  {
    float trench_width = merging_distance / float(shape.x);
    for (const auto &[key, path_cells] : breach_history)
    {
      if (path_cells.size() < 2) continue;
      std::vector<Point> pts;
      pts.reserve(path_cells.size());
      for (const auto &p : path_cells)
      {
        float x = (float(p.x) + 0.5f) / float(shape.x);
        float y = (float(p.y) + 0.5f) / float(shape.y);
        pts.push_back(Point(x, y, zb(p)));
      }
      Path river_path(pts);
      trench(zb,
             river_path,
             trench_width,
             /* enable_width_depth_scaling */ true,
             /* enable_width_distance_scaling */ false,
             /* enable_width_curvature_scaling */ false,
             /* curvature_radius_min */ 1.f,
             /* curv_width_ratio_min */ 0.5f,
             /* curv_width_ratio_max */ 2.f,
             RadialProfile::RP_SMOOTHSTEP_UPPER,
             /* radial_profile_parameter */ 2.f,
             ElevationLongitudinalProfile::ELP_DECREASING,
             /* elevation_shift */ 0.f,
             /* shift_ramp_start_ratio */ 0.f,
             /* shift_ramp_end_ratio */ 0.f,
             /* min_slope */ std::max(riverbed_talus, 1e-4f),
             /* k_neighbors */ 4,
             /* p_noise_r */ p_noise_r);
    }
  }

  return zb;
}

Array flow_fixing_drainage_basin(const Array        &z,
                                 FlowDirectionMethod fd_method,
                                 float               riverbed_talus,
                                 int                 iterations,
                                 bool                carve_riverbed,
                                 float               talus_riverbank,
                                 float               merging_distance,
                                 std::uint32_t       seed,
                                 float               noise_strength,
                                 const Array        *p_noise_x,
                                 const Array        *p_noise_y)
{
  if (!validate_non_empty(z)) return Array();
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return Array();

  Array zb = z;

  for (int it = 0; it < iterations; ++it)
  {
    DrainageBasinCellBased db(zb);

    if (fd_method == FlowDirectionMethod::FDM_D8)
    {
      db.compute_receivers(seed + it, noise_strength);
      auto [subroots, has_lake] = db.find_subroots();
      if (has_lake) db.remove_lakes(subroots);
    }
    else
    {
      db.compute_receivers_priority_flood();
    }

    db.update_traversals();

    // Correct upslopes by moving downstream from the highest headwater cells to
    // the outlets. Since db.traversals stores each outlet's tree in
    // upstream->downstream order (nodes ordered from leaves to outlet),
    // traversing in normal order visits every cell before its downstream
    // receiver.
    for (const auto &[outlet, traversal] : db.traversals)
    {
      for (const glm::ivec2 &i : traversal)
      {
        const glm::ivec2 &j = db.receivers(i);
        if (j == i) continue; // outlet cell

        int   dx = i.x - j.x;
        int   dy = i.y - j.y;
        float dist = (dx != 0 && dy != 0) ? M_SQRT2 : 1.f;
        float max_receiver_z = zb(i) - riverbed_talus * dist;

        // If the downstream receiver is higher than the current node, carve the
        // receiver down
        if (zb(j) > max_receiver_z) zb(j) = max_receiver_z;
      }
    }
  }

  // --- optional riverbed carving
  if (carve_riverbed)
  {
    return hmap::carve_riverbed(z,
                                zb,
                                talus_riverbank,
                                true, // smooth_river_bottom
                                merging_distance,
                                seed,
                                0.f, // riverbank_noise_ratio
                                p_noise_x,
                                p_noise_y);
  }
  else
  {
    return zb;
  }
}

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
                      int           fractalize_iterations,
                      float         fractalize_sigma,
                      int           decimate_target_points,
                      std::uint32_t fractalize_seed)
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

  // --- Dijkstra structures

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

  Mat<float>      dist_map(shape, std::numeric_limits<float>::max());
  Mat<int>        owner_map(shape, -1);
  Mat<glm::ivec2> prev_cell(shape, {-1, -1});

  std::priority_queue<DijkstraNode,
                      std::vector<DijkstraNode>,
                      std::greater<DijkstraNode>>
      pq;

  // ---- Initialize boundary outlets (Source ID = n_sinks)

  // Weight initial distance by border elevation so flow actively steers
  // towards lowest natural exit saddles/valleys along the boundary rather than
  // perpendicular lines.
  glm::vec2 z_range = zb.range();
  float     z_span = std::max(z_range.y - z_range.x, 1e-6f);

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

  // --- Initialize sinks (Source ID = 0 ... n_sinks - 1)

  for (int s = 0; s < n_sinks; ++s)
  {
    glm::ivec2 p = sinks[s];
    dist_map(p) = 0.f;
    owner_map(p) = s;
    prev_cell(p) = p;
    pq.push({0.f, p.x, p.y});
  }

  // Store candidate edges between meeting sources: key = (min(u,v), max(u,v))
  std::unordered_map<int64_t, MSTEdge> candidate_edges;

  auto make_key = [](int u, int v) -> int64_t
  {
    if (u > v) std::swap(u, v);
    return (static_cast<int64_t>(u) << 32) | static_cast<int64_t>(v);
  };

  // --- Multi-source Dijkstra expansion

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

      // When meeting a different source region, record candidate bridge edge
      if (nb_owner != -1 && nb_owner != cur_owner)
      {
        float   total_cost = dist_map(ci, cj) + dist_map(ni, nj) + cd[k];
        int64_t key = make_key(cur_owner, nb_owner);

        if (candidate_edges.find(key) == candidate_edges.end() ||
            total_cost < candidate_edges[key].cost)
        {
          // Reconstruct path from cur_owner to nb_owner through (ci, cj) - (ni,
          // nj)
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

      // Valley / concavity affinity: reduce cost in natural valleys and
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
  // --- (boundary_src_id) to guarantee all paths flow towards
  // --- domain boundaries.

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

  // We collect directed edges from child -> parent (pointing towards boundary)
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

  // Process paths in post-order (from deepest inner leaf sinks down towards the
  // boundary outlet) by reversing directed_paths (since BFS discovered
  // top-near-boundary first). This ensures that when a downstream sink's
  // elevation is lowered, all upstream sinks feeding into it have already
  // established their flow constraints or inherit downstream lowering. We do a
  // two-way pass (or iterative relaxation) along the tree branches to guarantee
  // a strictly decreasing profile along the entire chain of sinks to the
  // boundary.
  std::reverse(directed_paths.begin(), directed_paths.end());

  // Store the relaxed centerline elevations along each discrete cell path
  // without modifying zb directly if carving with trench on pristine z.
  Mat<float> cell_target_z(shape, std::numeric_limits<float>::max());
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      cell_target_z(i, j) = z(i, j);

  for (const auto &dp : directed_paths)
  {
    const auto &path = dp.path;
    if (path.size() < 2) continue;

    float min_d = std::max(minimum_depth, 0.f);
    cell_target_z(path.front()) = std::min(cell_target_z(path.front()),
                                           z(path.front()) - min_d);
    float current_z = cell_target_z(path.front());

    for (size_t idx = 1; idx < path.size(); ++idx)
    {
      glm::ivec2 curr = path[idx];
      glm::ivec2 prev = path[idx - 1];
      int        dx = curr.x - prev.x;
      int        dy = curr.y - prev.y;
      float      dist = (dx != 0 && dy != 0) ? M_SQRT2 : 1.f;

      // Decrement elevation along downstream flow direction
      current_z -= std::max(riverbed_talus, 1e-6f) * dist;

      // Ensure the elevation is strictly lower than the initial terrain
      // elevation
      float target_z = std::min(current_z, z(curr) - min_d);
      if (cell_target_z(curr) > target_z) cell_target_z(curr) = target_z;

      current_z = cell_target_z(curr);
    }
  }

  // --- Optional continuous riverbed carving with trench, or discrete incision
  // fallback
  if (carve_riverbed)
  {
    // Start carving on pristine terrain to avoid carving over an already
    // notched 1-pixel line
    zb = z;
    float trench_width = merging_distance / float(shape.x);

    // Carve riverbed using continuous trench along each path from inner sinks
    // to outlets
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
        pts.push_back(Point(x, y, cell_target_z(p)));
      }

      Path river_path(pts);

      if (decimate_target_points > 0 &&
          river_path.size() > static_cast<size_t>(decimate_target_points))
      {
        river_path = decimate_vw(river_path, decimate_target_points);
      }

      if (fractalize_iterations > 0)
      {
        river_path = fractalize_uniform(river_path,
                                        fractalize_iterations,
                                        fractalize_seed,
                                        fractalize_sigma);
      }

      trench(zb,
             river_path,
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
  }
  else
  {
    // Apply 1D discrete line lowering directly to zb
    for (const auto &dp : directed_paths)
    {
      for (const auto &p : dp.path)
      {
        if (zb(p) > cell_target_z(p)) zb(p) = cell_target_z(p);
      }
    }
  }

  return zb;
}

} // namespace hmap
