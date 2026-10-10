/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/curvature.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/morphology.hpp"

namespace hmap
{

Array valley_width(const Array &z, int ir, bool ridge_select)
{
  if (!validate_non_empty(z)) return Array();

  Array vw = z;

  if (ridge_select) vw *= -1.f;

  vw = -curvature_quadric(vw, ir, CurvatureType::CT_MEAN);
  vw = distance_transform_approx(vw);

  return vw;
}

} // namespace hmap

namespace hmap::gpu
{

Array valley_width(const Array &z, int ir, bool ridge_select)
{
  if (!validate_non_empty(z)) return Array();

  Array vw = z;

  if (ridge_select) vw *= -1.f;

  vw = -gpu::curvature_quadric(vw, ir, CurvatureType::CT_MEAN);
  vw = distance_transform_approx(vw);

  return vw;
}

} // namespace hmap::gpu
