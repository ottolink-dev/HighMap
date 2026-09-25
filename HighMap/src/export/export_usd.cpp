/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "hmm/src/heightmap.h"
#include "hmm/src/triangulator.h"

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/export.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/operator.hpp"

#include "lightusd.hh"
#include "stage.hh"
#include "usdGeom.hh"
#include "usda-writer.hh"
#include "usdc-writer.hh"

namespace hmap
{

namespace
{

lightusd::GeomBasisCurves build_usd_curves(const Path &path,
                                           size_t      path_idx,
                                           float       elevation_scaling,
                                           float       x_max,
                                           float       y_max)
{
  lightusd::GeomBasisCurves curves_prim;
  curves_prim.name = "Path_" + std::to_string(path_idx);
  curves_prim.type = lightusd::GeomBasisCurves::Type::Linear;

  const auto &x = path.get_x();
  const auto &y = path.get_y();
  const auto &v = path.get_values();

  std::vector<lightusd::value::point3f> pts(path.size());
  for (size_t i = 0; i < path.size(); ++i)
  {
    float val = (i < v.size()) ? v[i] : 0.f;
    pts[i] = {(1.f - y[i]) * y_max, elevation_scaling * val, x[i] * x_max};
  }

  curves_prim.points.set_value(pts);
  curves_prim.curveVertexCounts.set_value(
      std::vector<int>{static_cast<int>(path.size())});

  return curves_prim;
}

lightusd::GeomPoints build_usd_points(const Cloud &cloud,
                                      size_t       cloud_idx,
                                      float        elevation_scaling,
                                      float        x_max,
                                      float        y_max)
{
  lightusd::GeomPoints points_prim;
  points_prim.name = "Cloud_" + std::to_string(cloud_idx);

  const auto &x = cloud.get_x();
  const auto &y = cloud.get_y();
  const auto &v = cloud.get_values();

  std::vector<lightusd::value::point3f> pts(cloud.size());
  for (size_t i = 0; i < cloud.size(); ++i)
  {
    float val = (i < v.size()) ? v[i] : 0.f;
    pts[i] = {(1.f - y[i]) * y_max, elevation_scaling * val, x[i] * x_max};
  }

  points_prim.points.set_value(pts);
  return points_prim;
}

lightusd::GeomMesh build_usd_terrain_mesh(const Array &array,
                                          MeshType     mesh_type,
                                          float        elevation_scaling,
                                          float        max_error,
                                          bool         fit_boundaries)
{
  lightusd::GeomMesh mesh;
  mesh.name = "Terrain";

  float x_max = fit_boundaries
                    ? 1.f
                    : (array.shape.x > 0 ? 1.f - 1.f / (float)array.shape.x
                                         : 1.f);
  float y_max = fit_boundaries
                    ? 1.f
                    : (array.shape.y > 0 ? 1.f - 1.f / (float)array.shape.y
                                         : 1.f);

  std::vector<float> x = linspace(0.f, x_max, array.shape.x);
  std::vector<float> y = linspace(y_max, 0.f, array.shape.y);

  if (mesh_type == MeshType::TRI_OPTIMIZED)
  {
    auto p_hmap = std::make_shared<Heightmap>(array.shape.y,
                                              array.shape.x,
                                              array.get_vector());

    Triangulator tri(p_hmap);
    hmap::log::trace("remeshing USD terrain (Delaunay)");
    tri.Run(max_error, 0, 0);

    auto points = tri.Points(elevation_scaling);
    auto triangles = tri.Triangles();

    float ax = fit_boundaries
                   ? (array.shape.y > 1 ? 1.f / (float)(array.shape.y - 1)
                                        : 1.f)
                   : 1.f / (float)array.shape.y;
    float ay = fit_boundaries
                   ? (array.shape.x > 1 ? 1.f / (float)(array.shape.x - 1)
                                        : 1.f)
                   : 1.f / (float)array.shape.x;

    std::vector<lightusd::value::point3f> usd_points(points.size());
    for (size_t k = 0; k < points.size(); ++k)
    {
      usd_points[k] = {ay * points[k].y, points[k].z, ax * points[k].x};
    }

    std::vector<int32_t> face_vertex_counts(triangles.size(), 3);
    std::vector<int32_t> face_vertex_indices(triangles.size() * 3);

    for (size_t k = 0; k < triangles.size(); ++k)
    {
      face_vertex_indices[3 * k + 0] = static_cast<int32_t>(triangles[k].x);
      face_vertex_indices[3 * k + 1] = static_cast<int32_t>(triangles[k].y);
      face_vertex_indices[3 * k + 2] = static_cast<int32_t>(triangles[k].z);
    }

    mesh.points.set_value(usd_points);
    mesh.faceVertexCounts.set_value(face_vertex_counts);
    mesh.faceVertexIndices.set_value(face_vertex_indices);
  }
  else
  {
    // Regular TRI grid mesh
    std::vector<lightusd::value::point3f> usd_points(array.size());
    for (int j = 0; j < array.shape.y; ++j)
    {
      for (int i = 0; i < array.shape.x; ++i)
      {
        int k = array.linear_index(i, j);
        usd_points[k] = {y[j], elevation_scaling * array(i, j), x[i]};
      }
    }

    size_t               num_quads = (array.shape.x - 1) * (array.shape.y - 1);
    std::vector<int32_t> face_vertex_counts(num_quads * 2, 3);
    std::vector<int32_t> face_vertex_indices(num_quads * 6);

    size_t idx = 0;
    for (int j = 0; j < array.shape.y - 1; ++j)
    {
      for (int i = 0; i < array.shape.x - 1; ++i)
      {
        face_vertex_indices[idx++] = array.linear_index(i, j);
        face_vertex_indices[idx++] = array.linear_index(i, j + 1);
        face_vertex_indices[idx++] = array.linear_index(i + 1, j);

        face_vertex_indices[idx++] = array.linear_index(i + 1, j);
        face_vertex_indices[idx++] = array.linear_index(i, j + 1);
        face_vertex_indices[idx++] = array.linear_index(i + 1, j + 1);
      }
    }

    mesh.points.set_value(usd_points);
    mesh.faceVertexCounts.set_value(face_vertex_counts);
    mesh.faceVertexIndices.set_value(face_vertex_indices);
  }

  return mesh;
}

} // namespace

bool export_usd(const std::string        &fname,
                const Array              &elevation,
                const std::vector<Cloud> &clouds,
                const std::vector<Path>  &paths,
                MeshType                  mesh_type,
                float                     elevation_scaling,
                const std::string        &texture_fname,
                const std::string        &normal_map_fname,
                float                     max_error,
                bool                      fit_boundaries)
{
  if (!validate_non_empty(elevation)) return false;

  (void)texture_fname;
  (void)normal_map_fname;

  hmap::log::trace("exporting USD scene to [{}]", fname);

  float x_max = fit_boundaries ? 1.f
                               : (elevation.shape.x > 0
                                      ? 1.f - 1.f / (float)elevation.shape.x
                                      : 1.f);
  float y_max = fit_boundaries ? 1.f
                               : (elevation.shape.y > 0
                                      ? 1.f - 1.f / (float)elevation.shape.y
                                      : 1.f);

  lightusd::Stage stage;

  // --- Terrain mesh
  lightusd::GeomMesh terrain = build_usd_terrain_mesh(elevation,
                                                      mesh_type,
                                                      elevation_scaling,
                                                      max_error,
                                                      fit_boundaries);
  stage.add_root_prim(lightusd::Prim(terrain));

  // --- Clouds
  for (size_t i = 0; i < clouds.size(); ++i)
  {
    if (clouds[i].size() > 0)
    {
      lightusd::GeomPoints pts = build_usd_points(clouds[i],
                                                  i,
                                                  elevation_scaling,
                                                  x_max,
                                                  y_max);
      stage.add_root_prim(lightusd::Prim(pts));
    }
  }

  // --- Paths
  for (size_t i = 0; i < paths.size(); ++i)
  {
    if (paths[i].size() > 0)
    {
      lightusd::GeomBasisCurves cur = build_usd_curves(paths[i],
                                                       i,
                                                       elevation_scaling,
                                                       x_max,
                                                       y_max);
      stage.add_root_prim(lightusd::Prim(cur));
    }
  }

  std::filesystem::path fpath(fname);
  std::string           ext = fpath.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  std::string warn, err;
  bool        ok = false;

  if (ext == ".usdc")
  {
    ok = lightusd::usdc::SaveAsUSDCToFile(fname, stage, &warn, &err);
  }
  else if (ext == ".usdz")
  {
    std::map<std::string, std::vector<uint8_t>> asset_map;
    ok = lightusd::SaveAsUSDZToFile(fname, stage, asset_map, &warn, &err);
  }
  else
  {
    // Default to ASCII USDA format
    std::string out_fname = fname;
    if (ext != ".usda" && ext != ".usd")
    {
      out_fname += ".usda";
    }
    ok = lightusd::usda::SaveAsUSDA(out_fname, stage, &warn, &err);
  }

  if (!warn.empty())
  {
    hmap::log::warn("USD export warning: {}", warn);
  }

  if (!ok || !err.empty())
  {
    hmap::log::error("USD export failed: {}", err);
    return false;
  }

  return true;
}

} // namespace hmap
