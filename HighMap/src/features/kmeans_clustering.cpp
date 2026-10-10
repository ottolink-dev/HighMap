/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"

#include "dkm.hpp"

#include <format>

namespace hmap
{

Array kmeans_clustering2(const Array        &array1,
                         const Array        &array2,
                         int                 nclusters,
                         std::vector<Array> *p_scoring,
                         Array              *p_aggregate_scoring,
                         glm::vec2           weights,
                         std::uint32_t       seed)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  glm::ivec2 shape = array1.shape;
  Array      kmeans = Array(shape); // output

  // recast data
  std::vector<std::array<float, 2>> data = {};
  data.resize(shape.x * shape.y);

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      int k = i + j * shape.x;
      data[k][0] = weights.x * array1(i, j);
      data[k][1] = weights.y * array2(i, j);
    }

  dkm::clustering_parameters<float> parameters =
      dkm::clustering_parameters<float>(nclusters);
  parameters.set_random_seed(seed);
  auto dkm = dkm::kmeans_lloyd(data, parameters);

  // modify labelling to ensure it remains fairly consistent when the
  // data are modified (centroid are sorted by their coordinates)
  std::vector<int>   isort_rev(nclusters);
  std::vector<Point> centroids = {};

  for (auto &p : std::get<0>(dkm))
    centroids.push_back(Point(p[0], p[1]));

  sort_points(centroids);

  // TODO dirty
  for (int i = 0; i < nclusters; i++)
    for (int j = 0; j < nclusters; j++)
      if ((centroids[i].x == std::get<0>(dkm)[j][0]) &
          (centroids[i].y == std::get<0>(dkm)[j][1]))
        isort_rev[j] = i;

  for (size_t k = 0; k < std::get<1>(dkm).size(); k++)
  {
    int j = int(k / shape.x);
    int i = k - j * shape.x;
    kmeans(i, j) = isort_rev[std::get<1>(dkm)[k]];
  }

  // --- compute a score of belonging to a given cluster (see
  // https://datascience.stackexchange.com/questions/14435)

  // those scores are necessary for the aggregate score. If only the
  // aggregate score is requested, then work on a local vector, if
  // not, work on the input score vector of arrays
  std::vector<Array>  scores = {};
  std::vector<Array> *p_working_scores;

  if (p_scoring)
    p_working_scores = p_scoring;
  else
    p_working_scores = &scores;

  if (p_scoring || p_aggregate_scoring)
  {
    p_working_scores->clear();
    p_working_scores->reserve(nclusters);
    for (int r = 0; r < nclusters; r++)
      p_working_scores->push_back(Array(shape));

    for (int j = 0; j < shape.y; j++)
      for (int i = 0; i < shape.x; i++)
      {
        int k = i + j * shape.x; // linear index

        // normalization factor
        float sum = 0.f;
        for (int r = 0; r < nclusters; r++)
          sum += 1.f / std::hypot(data[k][0] - centroids[r].x,
                                  data[k][1] - centroids[r].y);

        // compute score for each cluster
        for (int r = 0; r < nclusters; r++)
        {
          float score = 1.f / std::hypot(data[k][0] - centroids[r].x,
                                         data[k][1] - centroids[r].y);
          score /= sum;
          p_working_scores->at(r)(i, j) = score;
        }
      }
  }

  // --- compute an aggregate score

  if (p_aggregate_scoring)
  {
    *p_aggregate_scoring = Array(shape);

    for (int j = 0; j < shape.y; j++)
      for (int i = 0; i < shape.x; i++)
      {
        float max = 0.f;
        int   rmax = 0;
        for (int r = 0; r < nclusters; r++)
          if (p_working_scores->at(r)(i, j) > max)
          {
            max = p_working_scores->at(r)(i, j);
            rmax = r;
          }

        (*p_aggregate_scoring)(i, j) = ((float)rmax + max) / (float)nclusters;
      }
  }

  return kmeans;
}

Array kmeans_clustering3(const Array        &array1,
                         const Array        &array2,
                         const Array        &array3,
                         int                 nclusters,
                         std::vector<Array> *p_scoring,
                         Array              *p_aggregate_scoring,
                         glm::vec3           weights,
                         std::uint32_t       seed)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2) ||
      !validate_same_shape(array1, array3))
    return Array();

  glm::ivec2 shape = array1.shape;
  Array      kmeans = Array(shape); // output

  // recast data
  std::vector<std::array<float, 3>> data = {};
  data.resize(shape.x * shape.y);

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      int k = i + j * shape.x;
      data[k][0] = weights.x * array1(i, j);
      data[k][1] = weights.y * array2(i, j);
      data[k][2] = weights.z * array3(i, j);
    }

  dkm::clustering_parameters<float> parameters =
      dkm::clustering_parameters<float>(nclusters);
  parameters.set_random_seed(seed);
  auto dkm = dkm::kmeans_lloyd(data, parameters);

  // modify labelling to ensure it remains fairly consistent when the
  // data are modified (centroid are sorted by their coordinates)
  std::vector<int>   isort_rev(nclusters);
  std::vector<Point> centroids = {};

  for (auto &p : std::get<0>(dkm))
    centroids.push_back(Point(p[0], p[1], p[2]));

  sort_points(centroids);

  // TODO dirty
  for (int i = 0; i < nclusters; i++)
    for (int j = 0; j < nclusters; j++)
      if ((centroids[i].x == std::get<0>(dkm)[j][0]) &
          (centroids[i].y == std::get<0>(dkm)[j][1]) &
          (centroids[i].v == std::get<0>(dkm)[j][2]))
        isort_rev[j] = i;

  for (size_t k = 0; k < std::get<1>(dkm).size(); k++)
  {
    int j = int(k / shape.x);
    int i = k - j * shape.x;
    kmeans(i, j) = isort_rev[std::get<1>(dkm)[k]];
  }

  // --- compute a score of belonging to a given cluster (see
  // https://datascience.stackexchange.com/questions/14435)

  std::vector<Array>  scores = {};
  std::vector<Array> *p_working_scores;

  if (p_scoring)
    p_working_scores = p_scoring;
  else
    p_working_scores = &scores;

  if (p_working_scores)
  {
    p_working_scores->reserve(nclusters);
    for (int r = 0; r < nclusters; r++)
      p_working_scores->push_back(Array(shape));

    for (int i = 0; i < shape.x; i++)
      for (int j = 0; j < shape.y; j++)
      {
        int k = i + j * shape.x; // linear index

        // normalization factor
        float sum = 0.f;
        for (int r = 0; r < nclusters; r++)
          sum += 1.f / std::pow(std::hypot(data[k][0] - centroids[r].x,
                                           data[k][1] - centroids[r].y,
                                           data[k][2] - centroids[r].v),
                                2);

        // compute score for each cluster
        for (int r = 0; r < nclusters; r++)
        {
          float score = 1.f / std::pow(std::hypot(data[k][0] - centroids[r].x,
                                                  data[k][1] - centroids[r].y,
                                                  data[k][2] - centroids[r].v),
                                       2);
          score /= sum;
          p_working_scores->at(r)(i, j) = score;
        }
      }
  }

  // --- compute an aggregate score

  if (p_aggregate_scoring)
  {
    *p_aggregate_scoring = Array(shape);

    for (int i = 0; i < shape.x; i++)
      for (int j = 0; j < shape.y; j++)
      {
        float max = 0.f;
        int   rmax = 0;
        for (int r = 0; r < nclusters; r++)
          if (p_working_scores->at(r)(i, j) > max)
          {
            max = p_working_scores->at(r)(i, j);
            rmax = r;
          }

        (*p_aggregate_scoring)(i, j) = ((float)rmax + max) / (float)nclusters;
      }
  }

  return kmeans;
}

} // namespace hmap
