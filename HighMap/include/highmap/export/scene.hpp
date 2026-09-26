/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file scene.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for 3D composite scene export functionalities (USD, etc.).
 *
 * This header declares functions and types related to exporting structured 3D
 * scenes composed of terrain meshes, point clouds, and path splines/curves.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <string>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/export/asset.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/path.hpp"

namespace hmap
{

/**
 * @brief Exports a composite scene (terrain, clouds, and paths) to Universal
 * Scene Description (USD).
 *
 * Generates a USD stage (.usda, .usdc, or .usdz) containing:
 * - Terrain mesh from `elevation` (with optional UV coordinates and elevation
 * scaling)
 * - Point collections from `clouds` as UsdGeomPoints
 * - Polylines/splines from `paths` as UsdGeomBasisCurves
 *
 * @param fname             Output file name (.usda, .usdc, or .usdz).
 * @param elevation         The heightmap array representing terrain elevation.
 * @param clouds            Vector of Cloud objects to export as points
 * primitives.
 * @param paths             Vector of Path objects to export as curve
 * primitives.
 * @param mesh_type         The type of mesh to generate for the terrain.
 * @param elevation_scaling Elevation scaling factor applied to the terrain and
 * 3D points.
 * @param texture_fname     Optional diffuse texture filename.
 * @param normal_map_fname  Optional normal map filename.
 * @param max_error         Maximum error for optimized Delaunay triangulation.
 * @param fit_boundaries    Scale domain coordinates to [0, 1].
 * @return                  `true` if export succeeded, `false` otherwise.
 *
 * **Example**
 * @include ex_export_usd.cpp
 */
bool export_usd(const std::string        &fname,
                const Array              &elevation,
                const std::vector<Cloud> &clouds = {},
                const std::vector<Path>  &paths = {},
                MeshType                  mesh_type = MeshType::TRI,
                float                     elevation_scaling = 0.2f,
                const std::string        &texture_fname = "",
                const std::string        &normal_map_fname = "",
                float                     max_error = 5e-4f,
                bool                      fit_boundaries = false);

} // namespace hmap
