/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file asset.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for 3D asset export functionalities.
 *
 * This header declares functions and types related to exporting heightmaps as
 * 3D mesh assets in various formats (OBJ, GLTF, PLY, STL, FBX, RAW, etc.) via
 * Assimp and native writers.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Enumeration for different mesh types.
 *
 * This enum defines the various types of mesh representations available. Each
 * type corresponds to a different way of constructing and representing mesh
 * data.
 */
enum MeshType : int
{
  TRI_OPTIMIZED, ///< Triangles with optimized Delaunay triangulation
  TRI,           ///< Triangle elements
};

/**
 * @brief Mapping between `MeshType` enum values and their plain text
 * descriptions.
 *
 * This static map provides a conversion between `MeshType` enum values and
 * their corresponding human-readable descriptions. It is used for displaying or
 * logging the mesh type in a human-friendly format.
 */
static std::map<MeshType, std::string> mesh_type_as_string = {
    {TRI_OPTIMIZED, "triangles (optimized)"},
    {TRI, "triangles"}};

/**
 * @brief Enumeration for asset export formats supported by Assimp.
 *
 * This enum lists the various file formats supported for asset export, as
 * recognized by the Assimp library. Each format is associated with a specific
 * file extension and usage.
 */
enum AssetExportFormat : int
{
  _3DS,    ///< Autodesk 3DS (legacy) - *.3ds
  _3MF,    ///< The 3MF-File-Format - *.3mf
  ASSBIN,  ///< Assimp Binary - *.assbin
  ASSXML,  ///< Assxml Document - *.assxml
  FXBA,    ///< Autodesk FBX (ascii) - *.fbx
  FBX,     ///< Autodesk FBX (binary) - *.fbx
  COLLADA, ///< COLLADA - Digital Asset Exchange Schema - *.dae
  X3D,     ///< Extensible 3D - *.x3d
  GLTF,    ///< GL Transmission Format - *.gltf
  GLB,     ///< GL Transmission Format (binary) - *.glb
  GTLF2,   ///< GL Transmission Format v. 2 - *.gltf
  GLB2,    ///< GL Transmission Format v. 2 (binary) - *.glb
  PLY,     ///< Stanford Polygon Library - *.ply
  PLYB,    ///< Stanford Polygon Library (binary) - *.ply
  STP,     ///< Step Files - *.stp
  STL,     ///< Stereolithography - *.stl
  STLB,    ///< Stereolithography (binary) - *.stl
  OBJ,     ///< Wavefront OBJ format - *.obj
  OBJNOMTL ///< Wavefront OBJ format without material file - *.obj
};

/**
 * @brief Mapping between asset export formats and their plain text
 * representations.
 *
 * This static map provides a mapping between `AssetExportFormat` enumeration
 * values and their corresponding plain text descriptions. Each entry includes a
 * human-readable format description, the format ID used by the Assimp library,
 * and the associated file extension. This mapping is used for converting
 * between enum values and their string representations in various export
 * scenarios.
 *
 * The format is structured as follows:
 * - Human-readable description of the format.
 * - Format ID as recognized by the Assimp library.
 * - File extension commonly used for that format.
 *
 * **Note**: For more details on the Assimp library formats, refer to [Assimp
 * Issue #2481](https://github.com/assimp/assimp/issues/2481).
 */
// clang-format off
static std::map<AssetExportFormat, std::vector<std::string> >
asset_export_format_as_string = {
	{_3DS, {"Autodesk 3DS (legacy) - *.3ds", "3ds", "3ds"}},
	{_3MF, {"The 3MF-File-Format - *.3mf", "3mf", "3mf"}},
	{ASSBIN, {"Assimp Binary - *.assbin", "assbin", "assbin"}},
	{ASSXML, {"Assxml Document - *.assxml", "assxml", "assxml"}},
	{FXBA, {"Autodesk FBX (ascii) - *.fbx", "fbxa", "fbx"}},
	{FBX, {"Autodesk FBX (binary) - *.fbx", "fbx", "fbx"}},
	{COLLADA, {"COLLADA - Digital Asset Exchange Schema - *.dae", "collada", "dae"}},
	{X3D, {"Extensible 3D - *.x3d", "x3d", "x3d"}},
	{GLTF, {"GL Transmission Format - *.gltf", "gltf", "gltf"}},
	{GLB, {"GL Transmission Format (binary) - *.glb", "glb", "glb"}},
	{GTLF2, {"GL Transmission Format v. 2 - *.gltf", "gltf2", "gltf"}},
	{GLB2, {"GL Transmission Format v. 2 (binary) - *.glb", "glb2", "glb"}},
	{PLY, {"Stanford Polygon Library - *.ply", "ply", "ply"}},
	{PLYB, {"Stanford Polygon Library (binary) - *.ply", "plyb", "ply"}},
	{STP, {"Step Files - *.stp", "stp", "stp"}},
	{STL, {"Stereolithography - *.stl", "stl", "stl"}},
	{STLB, {"Stereolithography (binary) - *.stl", "stlb", "stl"}},
	{OBJ, {"Wavefront OBJ format - *.obj", "obj", "obj"}},
	{OBJNOMTL, {"Wavefront OBJ format without material file - *.obj", "objnomtl", "obj"}},
};
// clang-format on

/**
 * @brief Exports a heightmap to various 3D file formats.
 *
 * This function exports the input heightmap array as a 3D asset in the
 * specified format. The export can include different mesh types, elevation
 * scaling, and optional texture and normal maps. The function supports
 * optimized Delaunay triangulation for mesh generation, with a configurable
 * maximum error.
 *
 * @param  fname             The name of the file to which the 3D asset will be
 *                           exported.
 * @param  array             The input heightmap array to be converted into a 3D
 *                           asset.
 * @param  mesh_type         The type of mesh to generate (see {@link
 *                           MeshType}).
 * @param  export_format     The format in which to export the asset (see {@link
 *                           AssetExportFormat}).
 * @param  elevation_scaling A scaling factor applied to the elevation values of
 *                           the heightmap. Default is 0.2f.
 * @param  texture_fname     The name of the texture file to be applied to the
 *                           asset (optional).
 * @param  normal_map_fname  The name of the normal map file to be applied to
 *                           the asset (optional).
 * @param  max_error         The maximum allowable error for optimized Delaunay
 *                           triangulation. Default is 5e-4f.
 * @param  fit_boundaries    If true, domain coordinates are scaled to [0, 1] x
 *                           [0, 1]. Default is false (domain scaled to [0, 1 -
 *                           1/shape.x] x [0, 1 -
 * 1/shape.y]).
 * @return                   `true` if the export is successful, `false`
 *                           otherwise.
 */
bool export_asset(const std::string &fname,
                  const Array       &array,
                  MeshType           mesh_type = MeshType::TRI,
                  AssetExportFormat  export_format = AssetExportFormat::GLB2,
                  float              elevation_scaling = 0.2f,
                  const std::string &texture_fname = "",
                  const std::string &normal_map_fname = "",
                  float              max_error = 5e-4f,
                  bool               fit_boundaries = false);

bool export_asset(const std::string &fname,
                  const Array       &array,
                  const Array       &mask,
                  AssetExportFormat  export_format = AssetExportFormat::GLB2,
                  float              elevation_scaling = 0.2f,
                  const std::string &texture_fname = "",
                  const std::string &normal_map_fname = "",
                  bool               fit_boundaries = false);

/**
 * @brief Exports 3D points with optional custom fields to an ASCII PLY file.
 *
 * This function writes vertex data (x, y, z) and optional additional per-point
 * attributes to a PLY file in ASCII format. Each key in @p custom_fields
 * defines a new property in the PLY header, and its associated vector provides
 * per-point values for that property.
 *
 * @param fname         The output file name.
 * @param x             Vector of x coordinates.
 * @param y             Vector of y coordinates.
 * @param z             Vector of z coordinates.
 * @param custom_fields A map of custom field names to their per-point float
 *                      values.
 *
 * **Example**
 * @include ex_export_points_to_ply.cpp
 */
void export_points_to_ply(
    const std::string                               &fname,
    const std::vector<float>                        &x,
    const std::vector<float>                        &y,
    const std::vector<float>                        &z,
    const std::map<std::string, std::vector<float>> &custom_fields = {});

/**
 * @brief Exports an array to a 16-bit 'raw' file format, commonly used for
 * Unity terrain imports.
 *
 * This function saves the input array to a file in a 16-bit 'raw' format, which
 * is suitable for importing heightmaps into Unity or other applications that
 * support this format. The array values are converted and written to the file
 * specified by `fname`.
 *
 * @param fname The name of the file to which the array will be exported.
 * @param array The input array containing the data to be exported.
 */
void write_raw_16bit(const std::string &fname, const Array &array);

} // namespace hmap
