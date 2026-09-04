/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "hmm/src/heightmap.h"
#include "hmm/src/triangulator.h"

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/export.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/operator.hpp"

#include <assimp/Exporter.hpp>
#include <assimp/material.h>
#include <assimp/mesh.h>
#include <assimp/metadata.h>
#include <assimp/scene.h>
#include <assimp/types.h>
#include <assimp/vector3.h>
#include <unordered_map>

namespace hmap
{

void helper_build_mesh_masked(aiMesh      *p_mesh,
                              const Array &array,
                              const Array &mask,
                              float        elevation_scaling)
{
  std::vector<float> x = linspace(1.f, 0.f, array.shape.x);
  std::vector<float> y = linspace(0.f, 1.f, array.shape.y);

  {
    std::uint32_t n_vertices = count_non_zero(mask);

    p_mesh->mVertices = new aiVector3D[n_vertices];
    p_mesh->mNumVertices = n_vertices;

    p_mesh->mTextureCoords[0] = new aiVector3D[n_vertices];
    p_mesh->mNumUVComponents[0] = 2;

    // --- generate vertices

    std::unordered_map<glm::ivec2, std::uint32_t, IVec2Hash, IVec2Eq> index_map;

    std::uint32_t k = 0;

    for (int j = 0; j < array.shape.y; j++)
      for (int i = 0; i < array.shape.x; i++)
      {
        if (mask(i, j))
        {
          p_mesh->mVertices[k] = aiVector3D(1.f - y[j],
                                            elevation_scaling * array(i, j),
                                            1.f - x[i]);

          p_mesh->mTextureCoords[0][k] = aiVector3D(1.f - x[i], y[j], 0.f);

          // store association
          index_map[{i, j}] = k;
          k++;
        }
      }

    // --- generate faces

    std::uint32_t n_faces = 0;

    for (int j = 0; j < array.shape.y - 1; j++)
      for (int i = 0; i < array.shape.x - 1; i++)
      {
        if (mask(i, j) && mask(i, j + 1) && mask(i + 1, j)) n_faces++;
        if (mask(i + 1, j) && mask(i, j + 1) && mask(i + 1, j + 1)) n_faces++;
      }

    p_mesh->mNumFaces = n_faces;
    p_mesh->mFaces = new aiFace[n_faces];
    p_mesh->mPrimitiveTypes = aiPrimitiveType_TRIANGLE;

    k = 0;

    for (int j = 0; j < array.shape.y - 1; j++)
      for (int i = 0; i < array.shape.x - 1; i++)
      {
        if (mask(i, j) && mask(i, j + 1) && mask(i + 1, j))
        {
          p_mesh->mFaces[k].mNumIndices = 3;
          p_mesh->mFaces[k].mIndices = new unsigned[3]{index_map[{i, j}],
                                                       index_map[{i, j + 1}],
                                                       index_map[{i + 1, j}]};
          k++;
        }

        if (mask(i + 1, j) && mask(i, j + 1) && mask(i + 1, j + 1))
        {
          p_mesh->mFaces[k].mNumIndices = 3;
          p_mesh->mFaces[k].mIndices = new unsigned[3]{
              index_map[{i + 1, j}],
              index_map[{i, j + 1}],
              index_map[{i + 1, j + 1}]};
          k++;
        }
      }
  }
}

void helper_build_mesh(aiMesh      *p_mesh,
                       const Array &array,
                       MeshType     mesh_type,
                       float        elevation_scaling,
                       float        max_error)
{
  std::vector<float> x = linspace(1.f, 0.f, array.shape.x);
  std::vector<float> y = linspace(0.f, 1.f, array.shape.y);

  switch (mesh_type)
  {
  case hmap::MeshType::TRI:
  {
    std::uint32_t n_vertices = array.size();
    std::uint32_t n_faces = 2 * (array.shape.x - 1) * (array.shape.y - 1);

    p_mesh->mVertices = new aiVector3D[n_vertices];
    p_mesh->mNumVertices = n_vertices;

    p_mesh->mTextureCoords[0] = new aiVector3D[n_vertices];
    p_mesh->mNumUVComponents[0] = 2;

    for (int j = 0; j < array.shape.y; j++)
      for (int i = 0; i < array.shape.x; i++)
      {
        int k = array.linear_index(i, j);

        p_mesh->mVertices[k] = aiVector3D(1.f - y[j],
                                          elevation_scaling * array(i, j),
                                          1.f - x[i]);

        p_mesh->mTextureCoords[0][k] = aiVector3D(1.f - x[i], y[j], 0.f);
      }

    p_mesh->mNumFaces = n_faces;
    p_mesh->mFaces = new aiFace[n_faces];
    p_mesh->mPrimitiveTypes = aiPrimitiveType_TRIANGLE;

    int k = 0;

    for (int j = 0; j < array.shape.y - 1; j++)
      for (int i = 0; i < array.shape.x - 1; i++)
      {
        p_mesh->mFaces[k].mNumIndices = 3;
        p_mesh->mFaces[k].mIndices = new unsigned[3]{
            (std::uint32_t)array.linear_index(i, j),
            (std::uint32_t)array.linear_index(i, j + 1),
            (std::uint32_t)array.linear_index(i + 1, j)};
        k++;

        p_mesh->mFaces[k].mNumIndices = 3;
        p_mesh->mFaces[k].mIndices = new unsigned[3]{
            (std::uint32_t)array.linear_index(i + 1, j),
            (std::uint32_t)array.linear_index(i, j + 1),
            (std::uint32_t)array.linear_index(i + 1, j + 1)};
        k++;
      }
  }
  break;

  case hmap::MeshType::TRI_OPTIMIZED:
  {
    auto p_hmap = std::make_shared<Heightmap>(array.shape.y,
                                              array.shape.x,
                                              array.get_vector());

    Triangulator tri(p_hmap);

    const int max_triangles = 0;
    const int max_points = 0;

    hmap::log::trace("remeshing (Delaunay)");
    tri.Run(max_error, max_triangles, max_points);

    auto points = tri.Points(elevation_scaling);
    auto triangles = tri.Triangles();

    float ax = 1.f / (float)array.shape.y;
    float ay = 1.f / (float)array.shape.x;

    std::uint32_t n_vertices = points.size();
    std::uint32_t n_faces = triangles.size();

    p_mesh->mVertices = new aiVector3D[n_vertices];
    p_mesh->mNumVertices = n_vertices;

    p_mesh->mTextureCoords[0] = new aiVector3D[n_vertices];
    p_mesh->mNumUVComponents[0] = 2;

    for (size_t k = 0; k < points.size(); k++)
    {
      p_mesh->mVertices[k] = aiVector3D(ay * points[k].y,
                                        points[k].z,
                                        ax * points[k].x);

      p_mesh->mTextureCoords[0][k] = aiVector3D(ax * points[k].x,
                                                1.f - ay * points[k].y,
                                                0.f);
    }

    p_mesh->mNumFaces = n_faces;
    p_mesh->mFaces = new aiFace[n_faces];
    p_mesh->mPrimitiveTypes = aiPrimitiveType_TRIANGLE;

    for (size_t k = 0; k < triangles.size(); k++)
    {
      p_mesh->mFaces[k].mNumIndices = 3;
      p_mesh->mFaces[k].mIndices = new unsigned[3]{
          (std::uint32_t)triangles[k].x,
          (std::uint32_t)triangles[k].y,
          (std::uint32_t)triangles[k].z};
    }
  }
  break;

  default: hmap::log::error("unknown mesh type"); break;
  }
}

bool export_asset(const std::string &fname,
                  const Array       &array,
                  MeshType           mesh_type,
                  AssetExportFormat  export_format,
                  float              elevation_scaling,
                  const std::string &texture_fname,
                  const std::string &normal_map_fname,
                  float              max_error)
{
  if (!validate_non_empty(array)) return false;

  hmap::log::trace("exporting asset, format [{}] aka [{}]",
                   asset_export_format_as_string.at(export_format)[0],
                   asset_export_format_as_string.at(export_format)[1]);

  aiScene *p_scene = new aiScene();
  aiMesh  *p_mesh = new aiMesh();

  // --- generate mesh

  std::vector<float> x = linspace(1.f, 0.f, array.shape.x);
  std::vector<float> y = linspace(0.f, 1.f, array.shape.y);

  aiMaterial *p_material = new aiMaterial();
  aiNode     *p_root = new aiNode();
  p_root->mNumMeshes = 1;
  p_root->mMeshes = new unsigned[]{0};

  p_scene->mNumMeshes = 1;
  p_scene->mMeshes = new aiMesh *[] { p_mesh };
  p_scene->mNumMaterials = 1;
  p_scene->mMaterials = new aiMaterial *[] { p_material };
  p_scene->mRootNode = p_root;
  p_scene->mMetaData = new aiMetadata();

  if (texture_fname != "")
  {
    aiString fname = aiString(texture_fname);
    p_material->AddProperty(&fname, AI_MATKEY_TEXTURE_DIFFUSE(0));
  }

  if (normal_map_fname != "")
  {
    aiString fname = aiString(normal_map_fname);
    p_material->AddProperty(&fname, AI_MATKEY_TEXTURE_NORMALS(0));
  }

  // --- build mesh

  helper_build_mesh(p_mesh, array, mesh_type, elevation_scaling, max_error);

  // --- export

  Assimp::Exporter exporter;

  std::string format_id = asset_export_format_as_string.at(export_format)[1];
  std::string format_extension = "." + asset_export_format_as_string.at(
                                           export_format)[2];

  aiReturn err = exporter.Export(p_scene,
                                 format_id,
                                 fname + format_extension,
                                 0);

  delete p_scene;

  if (err == aiReturn_SUCCESS)
    return true;
  else
  {
    hmap::log::error("failed to export asset");
    return false;
  }
}

bool export_asset(const std::string &fname,
                  const Array       &array,
                  const Array       &mask,
                  AssetExportFormat  export_format,
                  float              elevation_scaling,
                  const std::string &texture_fname,
                  const std::string &normal_map_fname)
{
  if (!validate_non_empty(array) || !validate_same_shape(array, mask))
    return false;

  hmap::log::trace("exporting asset, format [{}] aka [{}]",
                   asset_export_format_as_string.at(export_format)[0],
                   asset_export_format_as_string.at(export_format)[1]);

  aiScene *p_scene = new aiScene();
  aiMesh  *p_mesh = new aiMesh();

  // --- generate mesh

  std::vector<float> x = linspace(1.f, 0.f, array.shape.x);
  std::vector<float> y = linspace(0.f, 1.f, array.shape.y);

  aiMaterial *p_material = new aiMaterial();
  aiNode     *p_root = new aiNode();
  p_root->mNumMeshes = 1;
  p_root->mMeshes = new unsigned[]{0};

  p_scene->mNumMeshes = 1;
  p_scene->mMeshes = new aiMesh *[] { p_mesh };
  p_scene->mNumMaterials = 1;
  p_scene->mMaterials = new aiMaterial *[] { p_material };
  p_scene->mRootNode = p_root;
  p_scene->mMetaData = new aiMetadata();

  if (texture_fname != "")
  {
    aiString fname = aiString(texture_fname);
    p_material->AddProperty(&fname, AI_MATKEY_TEXTURE_DIFFUSE(0));
  }

  if (normal_map_fname != "")
  {
    aiString fname = aiString(normal_map_fname);
    p_material->AddProperty(&fname, AI_MATKEY_TEXTURE_NORMALS(0));
  }

  // --- build mesh

  helper_build_mesh_masked(p_mesh, array, mask, elevation_scaling);

  // --- export

  Assimp::Exporter exporter;

  std::string format_id = asset_export_format_as_string.at(export_format)[1];
  std::string format_extension = "." + asset_export_format_as_string.at(
                                           export_format)[2];

  aiReturn err = exporter.Export(p_scene,
                                 format_id,
                                 fname + format_extension,
                                 0);

  delete p_scene;

  if (err == aiReturn_SUCCESS)
    return true;
  else
  {
    hmap::log::error("failed to export asset");
    return false;
  }

  return true;
}

} // namespace hmap
