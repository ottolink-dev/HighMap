/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/curvature.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/morphology.hpp"
#include "highmap/range.hpp"
#include "highmap/selector.hpp"

namespace hmap::gpu
{

Array select_cavities(const Array   &array,
                      int            ir,
                      bool           concave,
                      bool           normalize,
                      ICommunicator &comm)
{
  if (!validate_non_empty(array)) return Array();

  Array array_smooth = array;
  gpu::smooth_cpulse(array_smooth, ir);
  Array c = -gpu::curvature_quadric(array_smooth, 0, CurvatureType::CT_MEAN);

  if (!concave) c *= -1.f;

  clamp_min(c, 0.f);

  if (normalize)
  {
    normalize_and_clip(c, 2.5f, true, true, comm);
  }

  return c;
}

Array select_slope(const Array      &z,
                   int               ir,
                   bool              normalize,
                   MinMaxKernel      kernel_type,
                   const TileRegion &region,
                   ICommunicator    &comm)
{
  if (!validate_non_empty(z)) return Array();

  Array out;

  if (ir == 0)
    out = hmap::gradient_norm(z);
  else
    out = gpu::morphological_gradient(z, ir, kernel_type);

  comm.sync_overlap_buffers(out, region, SyncOperation::SmoothBlend);

  if (normalize) normalize_and_clip(out, 2.5f, true, true, comm);

  return out;
}

Array select_valley(const Array   &z,
                    int            ir,
                    bool           ridge_select,
                    bool           normalize,
                    ICommunicator &comm)
{
  if (!validate_non_empty(z)) return Array();

  Array out = ridge_select ? gpu::morphological_top_hat(z, ir)
                           : gpu::morphological_black_hat(z, ir);

  if (normalize)
  {
    normalize_and_clip(out, 2.5f, true, true, comm);
  }

  return out;
}

} // namespace hmap::gpu
