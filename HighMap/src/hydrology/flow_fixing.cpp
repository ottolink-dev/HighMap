/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <queue>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/carving.hpp"
#include "highmap/filters.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/hydrology/drainage_basin_cell_based.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/profiles.hpp"

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

} // namespace hmap
