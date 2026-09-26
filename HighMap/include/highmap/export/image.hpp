/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file image.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for 2D image and raster export functionalities.
 *
 * This header declares functions for exporting heightmaps and arrays to image
 * formats (PNG, cubemaps, banners, splatmaps, tiled textures, ASCII art) and
 * reading raster images into arrays.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <string>
#include <vector>

#include <opencv2/imgcodecs.hpp>

#include "highmap/array.hpp"
#include "highmap/export/image_writer.hpp"

namespace hmap
{

enum Cmap : int; // highmap/colormaps.hpp

struct ComputeMode;  // highmap/virtual_array/virtual_array.hpp
struct VirtualArray; // highmap/virtual_array/virtual_array.hpp

/**
 * @brief Export a 2D array as an ASCII-art string representation.
 *
 * This function rescales the input array to a given export resolution, remaps
 * its values to the [0,1] range, and then converts each value to a character
 * from a user-provided character map. The characters are arranged row by row to
 * form a textual representation of the array, similar to rendering a heightmap
 * in ASCII.
 *
 * @param  array        The input 2D array to be exported.
 * @param  export_shape The desired shape (width, height) of the ASCII output.
 * @param  chars_map    A string containing characters ordered from "low"
 *                      to "high" intensity. For example: " .:-=+*#%@"
 *
 * @return              A string containing the ASCII-art representation of the
 *                      array.
 *
 * @note The function resamples the array using nearest-neighbor interpolation.
 * @note The y-axis is flipped so that higher indices appear lower in the text
 * output.
 */
std::string export_as_ascii(const Array      &array,
                            const glm::ivec2 &export_shape = {64, 64},
                            const std::string chars_map = " .:-=+*#%@");

/**
 * @brief Exports a 2D array as a cubemap texture with continuity enforcement
 * and overlapping regions.
 *
 * This function generates a cubemap texture from the input array `z`, resamples
 * the data to fit the cubemap resolution with optional overlapping regions, and
 * ensures seamless transitions between the six faces of the cubemap. The
 * cubemap can either be saved as a single texture or split into individual face
 * textures.
 *
 * @param fname              Output file name or base name for the cubemap
 *                           files.
 * @param z                  Input 2D array representing the data to be
 *                           converted into a cubemap.
 * @param cubemap_resolution Resolution (width and height) of each individual
 *                           face of the cubemap.
 * @param overlap            Fraction (0 to 1) of overlap between adjacent faces
 *                           to ensure smooth transitions.
 * @param ir                 Radius parameter for smoothing at triple corners.
 * @param cmap               Colormap to be applied when exporting the textures.
 * @param splitted           If true, exports each face of the cubemap as a
 *                           separate image; otherwise, exports the entire
 *                           cubemap as a single texture.
 * @param p_cubemap          Pointer to an optional output array where the final
 *                           cubemap data will be stored.
 *
 * The generated cubemap maintains continuity between faces, adjusting values at
 * overlapping regions and corners using smooth transitions. If the `splitted`
 * flag is set, six individual PNG images are generated for the cubemap faces
 * with appropriate suffixes appended to the base name. Otherwise, the entire
 * cubemap is exported as a single texture file.
 *
 * **Example**
 * @include ex_export_as_cubemap.cpp
 *
 * **Result**
 * @image html ex_export_as_cubemap.png
 */
void export_as_cubemap(const std::string &fname,
                       const Array       &z,
                       int                cubemap_resolution = 128,
                       float              overlap = 0.25f,
                       int                ir = 16,
                       Cmap               cmap = static_cast<Cmap>(0),
                       bool               splitted = false,
                       Array             *p_cubemap = nullptr);

/**
 * @brief Exports a set of arrays as a banner PNG image file.
 *
 * This function takes a vector of arrays and exports them as a single banner
 * PNG image. The arrays are displayed side by side in the image, using the
 * specified colormap `cmap`. Optionally, hillshading can be applied to enhance
 * the visual representation of the data.
 *
 * @param fname       The name of the file to which the banner image will be
 *                    exported.
 * @param arrays      A vector of arrays to be included in the banner image.
 * @param cmap        An integer representing the colormap to be applied to the
 *                    arrays.
 * @param hillshading A boolean flag to activate hillshading for enhanced visual
 *                    depth. Default is `false`.
 */
void export_banner_png(const std::string        &fname,
                       const std::vector<Array> &arrays,
                       int                       cmap,
                       bool                      hillshading = false,
                       bool                      normalize_arrays = false);

/**
 * @brief Exports an array using the generic ImageWriter pipeline.
 *
 * @param array  Input 2D array.
 * @param fname  Output file or directory path.
 * @param config Optional writer configuration.
 * @return       True if successful, false otherwise.
 */
bool export_image(const Array             &array,
                  const std::string       &fname,
                  const ImageWriterConfig &config = {});

/**
 * @brief Exports a VirtualArray incrementally using the generic ImageWriter
 * pipeline.
 *
 * @param va     Input VirtualArray.
 * @param fname  Output file or directory path.
 * @param config Optional writer configuration.
 * @return       True if successful, false otherwise.
 */
bool export_image(const VirtualArray      &va,
                  const std::string       &fname,
                  const ImageWriterConfig &config = {});

/**
 * @brief Exports a VirtualArray incrementally with custom compute mode.
 *
 * @param va     Input VirtualArray.
 * @param fname  Output file or directory path.
 * @param config Writer configuration.
 * @param cm     Compute mode for tile processing.
 * @return       True if successful, false otherwise.
 */
bool export_image(const VirtualArray      &va,
                  const std::string       &fname,
                  const ImageWriterConfig &config,
                  const ComputeMode       &cm);

/**
 * @brief Streams a VirtualArray tile-by-tile into an already opened
 * ImageWriter.
 *
 * @param va     Input VirtualArray.
 * @param writer Reference to target ImageWriter.
 * @return       True if successful, false otherwise.
 */
bool export_virtual_array(const VirtualArray &va, ImageWriter &writer);

/**
 * @brief Streams a VirtualArray tile-by-tile into an already opened ImageWriter
 * with custom compute mode.
 *
 * @param va     Input VirtualArray.
 * @param writer Reference to target ImageWriter.
 * @param cm     Compute mode for tile processing.
 * @return       True if successful, false otherwise.
 */
bool export_virtual_array(const VirtualArray &va,
                          ImageWriter        &writer,
                          const ComputeMode  &cm);

/**
 * @brief Exports the heightmap normal map as an 8-bit PNG file.
 *
 * This function generates a normal map from the input heightmap array and
 * exports it as an 8-bit PNG image. The normal map can be used in 3D rendering
 * engines to create realistic lighting and shading effects.
 *
 * @param fname The name of the file to which the normal map will be exported.
 * @param array The input heightmap array from which the normal map is derived.
 * @param depth The depth of the PNG image, e.g., `CV_8U` for 8-bit or `CV_16U`
 * for 16-bit. Default is `CV_8U`.
 *
 * **Example**
 * @include ex_export_normal_map.cpp
 */
void export_normal_map_png(const std::string &fname,
                           const Array       &array,
                           int                depth = CV_8U);

/**
 * @brief Exports four arrays as an RGBA PNG splatmap.
 *
 * This function combines four input arrays, representing the red (R), green
 * (G), blue (B), and alpha (A) channels, into a single RGBA PNG image. The
 * resulting splatmap can be used in applications like terrain texturing. The
 * PNG image is saved to the specified file name `fname`. Channels G, B, and A
 * are optional; if not provided, they will default to zero.
 *
 * @param fname The name of the file to which the RGBA splatmap will be
 *              exported.
 * @param p_r   Pointer to the array representing the red (R) channel.
 * @param p_g   Pointer to the array representing the green (G) channel. Default
 *              is `nullptr`.
 * @param p_b   Pointer to the array representing the blue (B) channel. Default
 *              is
 * `nullptr`.
 * @param p_a   Pointer to the array representing the alpha (A) channel. Default
 *              is `nullptr`.
 * @param depth The depth of the PNG image, e.g., `CV_8U` for 8-bit or `CV_16U`
 * for 16-bit. Default is `CV_8U`.
 *
 * **Example**
 * @include ex_export_splatmap_png_16bit.cpp
 *
 * **Result**
 * @image html ex_export_splatmap_png_16bit0.png
 * @image html ex_export_splatmap_png_16bit1.png
 */
void export_splatmap_png(const std::string &fname,
                         const Array       *p_r,
                         const Array       *p_g = nullptr,
                         const Array       *p_b = nullptr,
                         const Array       *p_a = nullptr,
                         int                depth = CV_8U);

/**
 * @brief Exports a 2D array as a set of grayscale PNG image tiles.
 *
 * This function divides a given 2D array into smaller rectangular tiles and
 * saves each tile as a grayscale PNG image file. Tiles are named using a
 * combination of the provided file name radical, tile indices, and file
 * extension.
 *
 * @param fname_radical           Base name (radical) for output image files.
 * @param fname_extension         File extension to use for exported images
 *                                (e.g.,
 * "png").
 * @param array                   The input 2D array to be tiled and exported.
 * @param tiling                  A 2D vector specifying the number of tiles in
 *                                the x and y directions.
 * @param leading_zeros           Number of digits used to pad the tile indices
 *                                in the filename.
 * @param depth                   Bit depth of the output PNG images (commonly 8
 *                                or 16).
 * @param overlapping_edges       If true, each tile includes an extra
 *                                row/column from neighboring tiles (for
 *                                overlap).
 * @param reverse_tile_y_indexing If true, Y tile indices are reversed (tile 0
 *                                is at the top).
 *
 * Each tile is extracted using slicing, adjusted for overlap if specified, and
 * then exported as an individual image file named with its tile indices. For
 * example, an output file might be named `radical_01_03.png`.
 *
 * **Example**
 * @include ex_export_tiled.cpp
 */
void export_tiled(const std::string &fname_radical,
                  const std::string &fname_extension,
                  const Array       &array,
                  const glm::ivec2  &tiling,
                  int                leading_zeros = 0,
                  int                depth = CV_8U,
                  bool               overlapping_edges = false,
                  bool               reverse_tile_y_indexing = false);

/**
 * @brief Reads an image file and converts it to a 2D array.
 *
 * This function uses the OpenCV `imread` function to load an image from the
 * specified file. The supported file formats are those recognized by OpenCV's
 * `imread` function, such as JPEG, PNG, BMP, and others. If the image is in
 * color, it is automatically converted to grayscale using the built-in OpenCV
 * codec converter. This conversion process may introduce artifacts depending on
 * the image's original format and content.
 *
 * @param  fname The name of the image file to be read.
 * @return       Array A 2D array containing the pixel values of the grayscale
 *               image.
 */
Array read_to_array(const std::string &fname,
                    bool               flip_j = false,
                    bool               remap = true);

} // namespace hmap
