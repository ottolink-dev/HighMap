/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>
#include <memory>
#include <vector>

#include "hmm/src/heightmap.h"
#include "hmm/src/triangulator.h"

#include "highmap/array.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/spatial/uv_tile.hpp"
#include "highmap/terrain_tri_mesh.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

TerrainTriMesh generate_terrain_tri_mesh_from_heightmap(const Array &z,
                                                        float        max_error,
                                                        int max_triangles,
                                                        int max_points)
{
  if (!validate_non_empty(z)) return TerrainTriMesh();

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

  return TerrainTriMesh(xyz);
}

TerrainTriMesh generate_terrain_tri_mesh_from_heightmap_random(
    const Array   &z,
    int            control_points_count,
    std::uint32_t  seed,
    const UvTile  &region,
    ICommunicator *p_comm)
{
  if (!validate_non_empty(z)) return TerrainTriMesh();

  SingleArrayComm default_comm;
  ICommunicator  &comm = p_comm ? *p_comm : default_comm;

  const BBox2D bbox = BBox2D::unit();

  Cloud cloud = random_cloud_jittered(control_points_count,
                                      {0.5f, 0.5f},
                                      {0.f, 0.f},
                                      seed,
                                      bbox);
  cloud.snap_points_to_bounding_box(bbox);

  const BBox2D core_bbox = region.core_uv_rect().to_bbox2d();

  Cloud tile_cloud = cloud;
  tile_cloud.set_values_from_array(z, region.uv_rect.to_bbox2d());
  const std::vector<glm::vec3> tile_xyz = tile_cloud.to_vec3();

  std::vector<glm::vec3> xyz = cloud.to_vec3();

  const bool is_single = (comm.size() == 1);

  for (size_t i = 0; i < xyz.size(); ++i)
  {
    const glm::vec3 &pt = xyz[i];

    bool owned = is_single;
    if (!owned)
    {
      const bool in_x = (region.is_left_border()
                             ? pt.x >= core_bbox.min.x - 1e-6f
                             : pt.x >= core_bbox.min.x) &&
                        (region.is_right_border()
                             ? pt.x <= core_bbox.max.x + 1e-6f
                             : pt.x < core_bbox.max.x);
      const bool in_y = (region.is_bottom_border()
                             ? pt.y >= core_bbox.min.y - 1e-6f
                             : pt.y >= core_bbox.min.y) &&
                        (region.is_top_border()
                             ? pt.y <= core_bbox.max.y + 1e-6f
                             : pt.y < core_bbox.max.y);
      owned = in_x && in_y;
    }

    const float local_z = owned ? tile_xyz[i].z : 0.f;
    xyz[i].z = comm.all_reduce(local_z, ReduceOp::Sum);
  }

  if (comm.rank() == 0) return TerrainTriMesh(xyz);

  return TerrainTriMesh();
}

} // namespace hmap
