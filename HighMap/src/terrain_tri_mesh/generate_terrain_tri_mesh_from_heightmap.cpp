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
    const Array  &z,
    int           control_points_count,
    std::uint32_t seed)
{
  if (!validate_non_empty(z)) return TerrainTriMesh();

  const BBox2D bbox = BBox2D::unit();

  Cloud cloud = random_cloud_jittered(control_points_count,
                                      {0.5f, 0.5f},
                                      {0.f, 0.f},
                                      seed,
                                      bbox);
  cloud.snap_points_to_bounding_box(bbox);
  cloud.set_values_from_array(z, bbox);

  return TerrainTriMesh(cloud.to_vec3());
}

} // namespace hmap

namespace hmap::va
{

TerrainTriMesh generate_terrain_tri_mesh_from_heightmap_random(
    const ComputeMode  &cm,
    const VirtualArray &z,
    int                 control_points_count,
    std::uint32_t       seed)
{
  const BBox2D bbox = BBox2D::unit();

  Cloud cloud = random_cloud_jittered(control_points_count,
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
          const hmap::UvTile &region)
      {
        auto [pa_z] = unpack<1>(p_arrays_in);
        cloud.set_values_from_array(*pa_z, region.uv_rect.to_bbox2d());
      },
      cm_sampling);

  return TerrainTriMesh(cloud.to_vec3());
}

} // namespace hmap::va
