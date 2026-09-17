/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/kd_tree.hpp"
#include "highmap/geometry/point.hpp"

#include "nanoflann.hpp"

namespace hmap
{

// ==========================================================================
//  Internal Adaptors & Implementation
// ==========================================================================

namespace
{

struct PointVectorAdaptor
{
  const Point *points{nullptr};
  size_t       count{0};

  PointVectorAdaptor(const Point *p, size_t c) : points(p), count(c)
  {
  }

  size_t kdtree_get_point_count() const
  {
    return count;
  }

  float kdtree_get_pt(size_t idx, int dim) const
  {
    return (dim == 0) ? points[idx].x : points[idx].y;
  }

  template <class BBOX> bool kdtree_get_bbox(BBOX &) const
  {
    return false;
  }
};

struct XYVectorAdaptor
{
  const std::vector<float> &x;
  const std::vector<float> &y;

  XYVectorAdaptor(const std::vector<float> &x_, const std::vector<float> &y_)
      : x(x_), y(y_)
  {
  }

  size_t kdtree_get_point_count() const
  {
    return std::min(x.size(), y.size());
  }

  float kdtree_get_pt(size_t idx, int dim) const
  {
    return (dim == 0) ? x[idx] : y[idx];
  }

  template <class BBOX> bool kdtree_get_bbox(BBOX &) const
  {
    return false;
  }
};

} // namespace

struct KDTree::Impl
{
  virtual ~Impl() = default;
  virtual size_t size() const = 0;
  virtual float  get_x(size_t i) const = 0;
  virtual float  get_y(size_t i) const = 0;

  virtual std::pair<size_t, float> nearest_with_distance_squared(
      float x,
      float y) const = 0;

  virtual void neighbor_search(float                x,
                               float                y,
                               size_t               k,
                               std::vector<size_t> &indices,
                               std::vector<float>  &distances) const = 0;

  virtual std::vector<std::pair<size_t, float>> radius_search(
      float x,
      float y,
      float radius) const = 0;
};

namespace
{

template <typename AdaptorT> struct NanoImpl : public KDTree::Impl
{
  AdaptorT adaptor;

  using NanoKDTree = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<float, AdaptorT>,
      AdaptorT,
      2>;

  std::unique_ptr<NanoKDTree> index;

  template <typename... Args>
  explicit NanoImpl(Args &&...args) : adaptor(std::forward<Args>(args)...)
  {
    if (adaptor.kdtree_get_point_count() > 0)
    {
      index = std::make_unique<NanoKDTree>(
          2,
          adaptor,
          nanoflann::KDTreeSingleIndexAdaptorParams(10));
      index->buildIndex();
    }
  }

  size_t size() const override
  {
    return adaptor.kdtree_get_point_count();
  }

  float get_x(size_t i) const override
  {
    return adaptor.kdtree_get_pt(i, 0);
  }

  float get_y(size_t i) const override
  {
    return adaptor.kdtree_get_pt(i, 1);
  }

  std::pair<size_t, float> nearest_with_distance_squared(float x,
                                                         float y) const override
  {
    if (!index || size() == 0) return {0, 0.f};

    size_t out_index = 0;
    float  out_dist_sq = 0.f;

    nanoflann::KNNResultSet<float> result_set(1);
    result_set.init(&out_index, &out_dist_sq);

    float query_pt[2] = {x, y};
    index->findNeighbors(result_set, query_pt, nanoflann::SearchParameters());

    return {out_index, out_dist_sq};
  }

  void neighbor_search(float                x,
                       float                y,
                       size_t               k,
                       std::vector<size_t> &indices,
                       std::vector<float>  &distances) const override
  {
    const size_t n = size();
    if (!index || n == 0 || k == 0)
    {
      indices.clear();
      distances.clear();
      return;
    }

    size_t k_clamped = std::min(k, n);
    indices.resize(k_clamped);
    distances.resize(k_clamped);

    nanoflann::KNNResultSet<float> result_set(k_clamped);
    result_set.init(indices.data(), distances.data());

    float query_pt[2] = {x, y};
    index->findNeighbors(result_set, query_pt, nanoflann::SearchParameters());
  }

  std::vector<std::pair<size_t, float>> radius_search(
      float x,
      float y,
      float radius) const override
  {
    if (!index || radius <= 0.f || size() == 0) return {};

    std::vector<nanoflann::ResultItem<uint32_t, float>> matches;
    float                                               query_pt[2] = {x, y};
    index->radiusSearch(query_pt,
                        radius * radius,
                        matches,
                        nanoflann::SearchParameters());

    std::vector<std::pair<size_t, float>> result;
    result.reserve(matches.size());
    for (const auto &m : matches)
    {
      result.emplace_back(static_cast<size_t>(m.first), m.second);
    }
    return result;
  }
};

} // namespace

// ==========================================================================
//  KDTree Constructors & Destructor
// ==========================================================================

KDTree::KDTree() = default;

KDTree::~KDTree() = default;

KDTree::KDTree(const std::vector<float> &x, const std::vector<float> &y)
    : p_impl(std::make_unique<NanoImpl<XYVectorAdaptor>>(x, y))
{
}

KDTree::KDTree(const std::vector<Point> &points)
    : p_impl(std::make_unique<NanoImpl<PointVectorAdaptor>>(points.data(),
                                                            points.size()))
{
}

KDTree::KDTree(const Cloud &cloud)
    : p_impl(std::make_unique<NanoImpl<PointVectorAdaptor>>(cloud.data(),
                                                            cloud.size()))
{
}

KDTree::KDTree(KDTree &&other) noexcept = default;

KDTree &KDTree::operator=(KDTree &&other) noexcept = default;

// ==========================================================================
//  KDTree Capacity & Status
// ==========================================================================

bool KDTree::empty() const noexcept
{
  return !p_impl || p_impl->size() == 0;
}

size_t KDTree::size() const noexcept
{
  return p_impl ? p_impl->size() : 0;
}

// ==========================================================================
//  KDTree Queries
// ==========================================================================

glm::vec2 KDTree::compute_neighbor_distance_range(size_t k_neighbors) const
{
  if (!p_impl || p_impl->size() == 0) return {0.f, 0.f};

  const size_t n = p_impl->size();
  size_t       k_clamped = std::clamp(k_neighbors, size_t{1}, n);

  float dmin_sq = std::numeric_limits<float>::max();
  float dmax_sq = 0.f;

  std::vector<size_t> indices;
  std::vector<float>  distances;

  for (size_t i = 0; i < n; ++i)
  {
    float px = p_impl->get_x(i);
    float py = p_impl->get_y(i);
    p_impl->neighbor_search(px, py, k_clamped, indices, distances);

    for (const auto &d_sq : distances)
    {
      dmax_sq = std::max(dmax_sq, d_sq);
      if (d_sq > 1e-12f)
      {
        dmin_sq = std::min(dmin_sq, d_sq);
      }
    }
  }

  if (dmin_sq == std::numeric_limits<float>::max())
  {
    dmin_sq = dmax_sq;
  }

  return {std::sqrt(dmin_sq), std::sqrt(dmax_sq)};
}

size_t KDTree::nearest(float x_query, float y_query) const
{
  return this->nearest_with_distance_squared(x_query, y_query).first;
}

size_t KDTree::nearest(const Point &p) const
{
  return this->nearest(p.x, p.y);
}

size_t KDTree::nearest(const glm::vec2 &xy) const
{
  return this->nearest(xy.x, xy.y);
}

std::pair<size_t, float> KDTree::nearest_with_distance_squared(
    float x_query,
    float y_query) const
{
  if (!p_impl) return {0, 0.f};
  return p_impl->nearest_with_distance_squared(x_query, y_query);
}

std::pair<size_t, float> KDTree::nearest_with_distance_squared(
    const Point &p) const
{
  return this->nearest_with_distance_squared(p.x, p.y);
}

std::pair<size_t, float> KDTree::nearest_with_distance_squared(
    const glm::vec2 &xy) const
{
  return this->nearest_with_distance_squared(xy.x, xy.y);
}

void KDTree::neighbor_search(float                x_query,
                             float                y_query,
                             size_t               k_neighbors,
                             std::vector<size_t> &indices,
                             std::vector<float>  &distances) const
{
  if (!p_impl)
  {
    indices.clear();
    distances.clear();
    return;
  }
  p_impl->neighbor_search(x_query, y_query, k_neighbors, indices, distances);
}

void KDTree::neighbor_search(const Point         &p,
                             size_t               k_neighbors,
                             std::vector<size_t> &indices,
                             std::vector<float>  &distances) const
{
  this->neighbor_search(p.x, p.y, k_neighbors, indices, distances);
}

void KDTree::neighbor_search(const glm::vec2     &xy,
                             size_t               k_neighbors,
                             std::vector<size_t> &indices,
                             std::vector<float>  &distances) const
{
  this->neighbor_search(xy.x, xy.y, k_neighbors, indices, distances);
}

std::vector<std::pair<size_t, float>> KDTree::radius_search(float x_query,
                                                            float y_query,
                                                            float radius) const
{
  if (!p_impl) return {};
  return p_impl->radius_search(x_query, y_query, radius);
}

std::vector<std::pair<size_t, float>> KDTree::radius_search(const Point &p,
                                                            float radius) const
{
  return this->radius_search(p.x, p.y, radius);
}

std::vector<std::pair<size_t, float>> KDTree::radius_search(const glm::vec2 &xy,
                                                            float radius) const
{
  return this->radius_search(xy.x, xy.y, radius);
}

} // namespace hmap
