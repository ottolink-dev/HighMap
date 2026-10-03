/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file features.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file defining a collection of functions for terrain analysis
 * and feature extraction from heightmaps.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <array>

#include "highmap/array.hpp"

/**
 * @brief Packs eight 2-bit values into a 16-bit integer.
 *
 * This macro takes eight 2-bit values (`a` through `h`) and packs them into a
 * single 16-bit integer. The bits are shifted into position according to their
 * specified offsets to form the final packed value.
 *
 * @param  a The first 2-bit value (will be shifted left by 15 bits).
 * @param  b The second 2-bit value (will be shifted left by 13 bits).
 * @param  c The third 2-bit value (will be shifted left by 11 bits).
 * @param  d The fourth 2-bit value (will be shifted left by 9 bits).
 * @param  e The fifth 2-bit value (will be shifted left by 7 bits).
 * @param  f The sixth 2-bit value (will be shifted left by 5 bits).
 * @param  g The seventh 2-bit value (will be shifted left by 3 bits).
 * @param  h The eighth 2-bit value (will be shifted left by 1 bit).
 *
 * @return   A 16-bit integer where the input values have been packed together
 *           according to their respective shifts.
 *
 * @note Each parameter (`a` through `h`) should only occupy 2 bits (values
 * between 0 and 3). If the values exceed this range, the result may be
 * incorrect.
 */
#define HMAP_PACK8(a, b, c, d, e, f, g, h)                                     \
  ((a << 15) + (b << 13) + (c << 11) + (d << 9) + (e << 7) + (f << 5) +        \
   (g << 3) + (h << 1))

namespace hmap
{

/**
 * @brief Identifies and labels connected components within a binary or labeled
 * array, with optional filtering by size.
 *
 * Connected-component labeling is a technique used to identify clusters of
 * connected pixels (or components) in an array. This function can be used in
 * image processing and spatial analysis to isolate regions of interest, such as
 * detecting distinct objects or areas within a heightmap.
 *
 * **Usage**
 * - Use this function to label distinct features within an array, such as
 * individual water bodies, forest patches, or elevation features.
 * - Apply the surface threshold to filter out small, insignificant components
 * that might be noise or irrelevant.
 *
 * @param  array             The input array where connected components are to
 *                           be identified.
 * @param  surface_threshold The minimum number of pixels a component must have
 *                           to be retained. Components smaller than this
 *                           threshold will be removed. The default value is 0
 *                           (no filtering).
 * @param  background_value  The value used to represent background pixels,
 *                           which are not part of any component. Default is 0.
 * @return                   Array An array with labeled connected components,
 *                           where each component is assigned a unique
 *                           identifier.
 *
 *
 * **Example**
 * @include ex_connected_components.cpp
 *
 * **Result**
 * @image html ex_connected_components0.png
 * @image html ex_connected_components1.png
 *
 * See unit tests: @ref test_connected_components.cpp
 */
Array connected_components(
    const Array                           &array,
    float                                  surface_threshold = 0.f,
    float                                  background_value = 0.f,
    std::map<float, float>                *p_surfaces = nullptr,
    std::map<float, std::array<float, 2>> *p_centroids = nullptr);

/**
 * @brief Classifies terrain into geomorphological features based on the
 * geomorphons method.
 *
 * The geomorphons method classifies each point in a terrain into one of several
 * landform types (e.g., ridge, valley, plain) based on the surrounding
 * topography. This classification is useful for automated landform mapping and
 * terrain analysis.
 *
 * **Usage**
 * - Use this function to automatically categorize terrain into different
 * geomorphological units.
 * - Useful in large-scale landform mapping and environmental modeling.
 *
 * @param  array   The input array representing the terrain elevation data.
 * @param  irmin   The minimum radius (in pixels) for considering the
 *                 surrounding area during classification.
 * @param  irmax   The maximum radius (in pixels) for considering the
 *                 surrounding area during classification.
 * @param  epsilon The slope tolerance that defines 'flat' regions, affecting
 *                 the classification.
 * @return         Array An output array with each pixel classified into a
 *                 geomorphological feature.
 *
 * **Example**
 * @include ex_geomorphons.cpp
 *
 * **Result**
 * @image html ex_geomorphons0.png
 * @image html ex_geomorphons1.png
 */
Array geomorphons(const Array &array, int irmin, int irmax, float epsilon);

/**
 * @brief Performs k-means clustering on two input arrays, grouping similar data
 * points into clusters.
 *
 * K-means clustering is a popular unsupervised learning algorithm used to
 * partition data into clusters. This function applies k-means clustering to two
 * arrays, which might represent different terrain attributes or environmental
 * variables.
 *
 * **Usage**:
 * - Use this function to identify regions with similar characteristics based on
 * multiple terrain features.
 * - Useful in ecological modeling, land cover classification, and resource
 * management.
 *
 * @param[in]  array1              The first input array for clustering,
 *                                 typically representing one attribute of the
 *                                 terrain.
 * @param[in]  array2              The second input array for clustering,
 *                                 representing another attribute.
 * @param[in]  nclusters           The number of clusters (k) to be formed.
 * @param[out] p_scoring           (optional) A pointer to a vector of arrays
 *                                 where the clustering scores will be stored.
 *                                 Pass nullptr if scoring is not required.
 * @param[out] p_aggregate_scoring (optional) A pointer to an array where the
 *                                 aggregate score across all clusters will be
 *                                 stored. Pass nullptr if not required.
 * @param[in]  weights             A vector of two floats representing the
 *                                 weight given to
 *                                 `array1` and `array2`. The default weights
 *                                 are {1.f, 1.f}.
 * @param[in]  seed                A seed value for random number generation,
 *                                 ensuring reproducibility of the clustering
 *                                 results. The default value is 1.
 * @return                         Array An array representing the clustered
 *                                 data, with each pixel assigned to a cluster.
 *
 * @return                         An Array containing the clustered data
 *                                 resulting from the k-means algorithm.
 *
 * **Example**
 * @include ex_kmeans_clustering.cpp
 *
 * **Result**
 * @image html ex_kmeans_clustering0.png
 * @image html ex_kmeans_clustering1.png
 * @image html ex_kmeans_clustering2.png
 * @image html ex_kmeans_clustering3.png
 * @image html ex_kmeans_clustering4.png
 * @image html ex_kmeans_clustering5.png
 */
Array kmeans_clustering2(const Array        &array1,
                         const Array        &array2,
                         int                 nclusters,
                         std::vector<Array> *p_scoring = nullptr,
                         Array              *p_aggregate_scoring = nullptr,
                         glm::vec2           weights = {1.f, 1.f},
                         std::uint32_t       seed = 1);

/**
 * @brief Performs k-means clustering on three input arrays, providing more
 * detailed cluster analysis by considering an additional dimension.
 *
 * This version of k-means clustering includes a third array, enabling more
 * complex clustering based on three terrain attributes or environmental
 * variables.
 *
 * **Usage**
 * - Apply this function to analyze the relationships between three different
 * terrain attributes, revealing complex patterns.
 * - Useful in multi-dimensional environmental modeling and resource management.
 *
 * @param[in]  array1              The first input array for clustering.
 * @param[in]  array2              The second input array for clustering.
 * @param[in]  array3              The third input array for clustering.
 * @param[in]  nclusters           The number of clusters (k) to be formed.
 * @param[out] p_scoring           (optional) A pointer to a vector of arrays
 *                                 where the clustering scores will be stored.
 *                                 Pass nullptr if scoring is not required.
 * @param[out] p_aggregate_scoring (optional) A pointer to an array where the
 *                                 aggregate score across all clusters will be
 *                                 stored. Pass nullptr if not required.
 * @param[in]  weights             A vector of three floats representing the
 *                                 weight given to
 *                                 `array1`, `array2`, and `array3`. The default
 *                                 weights are {1.f, 1.f, 1.f}.
 * @param[in]  seed                A seed value for random number generation,
 *                                 ensuring reproducibility of the clustering
 *                                 results. The default value is 1.
 * @return                         Array An array representing the clustered
 *                                 data, with each pixel assigned to a cluster.
 *
 * **Example**
 * @include ex_kmeans_clustering.cpp
 *
 * **Result**
 * @image html ex_kmeans_clustering0.png
 * @image html ex_kmeans_clustering1.png
 * @image html ex_kmeans_clustering2.png
 * @image html ex_kmeans_clustering3.png
 * @image html ex_kmeans_clustering4.png
 * @image html ex_kmeans_clustering5.png
 */
Array kmeans_clustering3(const Array        &array1,
                         const Array        &array2,
                         const Array        &array3,
                         int                 nclusters,
                         std::vector<Array> *p_scoring = nullptr,
                         Array              *p_aggregate_scoring = nullptr,
                         glm::vec3           weights = {1.f, 1.f, 1.f},
                         std::uint32_t       seed = 1);

} // namespace hmap
