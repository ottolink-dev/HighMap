/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/curvature.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"

#include "appMetrics/core.h"
#include "appMetrics/heightfield.h"
#include "appMetrics/scalarfield.h"

namespace hmap
{

Array curvature_quadric(const Array &z, int ir, CurvatureType curvature_type)
{
  if (!validate_non_empty(z)) return Array();

  const glm::ivec2 &shape = z.shape;
  const int         w = 2 * ir + 1; // window size

  // --- Adaptor to terrain-descriptors classes

  std::vector<double> storage;
  storage.reserve(shape.x * shape.y);

  for (const auto &v : z.vector)
    storage.push_back(double(v));

  ScalarField2 field = ScalarField2(
      Box2(double(shape.x - 1), double(shape.y - 1)),
      shape.x,
      shape.y,
      storage);
  HeightField  h(field);
  ScalarField2 out;

  // --- Compute curvatures

  switch (curvature_type)
  {
  case hmap::CurvatureType::CT_MIN:
    out = h.Curvature(HeightField::CurvatureType::MIN, w);
    break;
    //
  case hmap::CurvatureType::CT_MAX:
    out = h.Curvature(HeightField::CurvatureType::MAX, w);
    break;
    //
  case hmap::CurvatureType::CT_MEAN:
    out = h.Curvature(HeightField::CurvatureType::MEAN, w);
    break;
    //
  case hmap::CurvatureType::CT_GAUSSIAN:
    out = h.Curvature(HeightField::CurvatureType::GAUSSIAN, w);
    break;
    //
  case hmap::CurvatureType::CT_PROFILE:
    out = h.Curvature(HeightField::CurvatureType::PROFILE, w);
    break;
    //
  case hmap::CurvatureType::CT_CONTOUR:
    out = h.Curvature(HeightField::CurvatureType::CONTOUR, w);
    break;
    //
  case hmap::CurvatureType::CT_TANGENTIAL:
    out = h.Curvature(HeightField::CurvatureType::TANGENTIAL, w);
    break;
    //
    break;
  default:
  {
    hmap::log::trace("unknown CurvatureType");
    return Array(shape);
  }
  }

  // --- Output

  Array                      c(shape);
  const std::vector<double> &values = out.values();

  for (size_t k = 0; k < values.size(); ++k)
    c.vector[k] = float(values[k]);

  return c;
}

} // namespace hmap

namespace hmap::gpu
{

Array curvature_quadric(const Array  &z,
                        int           ir,
                        CurvatureType curvature_type,
                        bool          approx_algo)
{
  if (!validate_non_empty(z)) return Array();

  const glm::ivec2 shape = z.shape;
  int              type_id = static_cast<int>(curvature_type);

  Array out(shape);

  auto run = clwrapper::Run("curvature_quadric");

  Array zf;

  if (approx_algo)
  {
    zf = z;
    gpu::smooth_cpulse(zf, ir);
    ir = 0;
    run.bind_imagef("z", zf.vector, shape.x, shape.y);
  }
  else
  {
    run.bind_imagef("z", z.vector, shape.x, shape.y);
  }

  run.bind_imagef("out", out.vector, shape.x, shape.y, true);

  run.bind_arguments(shape.x, shape.y, ir, type_id);

  run.execute({shape.x, shape.y});

  run.read_imagef("out");

  if (ir == 0) extrapolate_borders(out);

  return out;
}

} // namespace hmap::gpu
