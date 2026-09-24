/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/morphology.hpp"

namespace hmap
{

Array morphological_operators(const Array        &array,
                              int                 ir,
                              MorphologyOperation operation,
                              MinMaxKernel        kernel_type)
{
  if (!validate_non_empty(array)) return Array();

  (void)kernel_type;

  switch (operation)
  {
  case MorphologyOperation::MO_BORDER:
    return border(array, ir);
    //
  case MorphologyOperation::MO_CLOSING:
    return closing(array, ir);
    //
  case MorphologyOperation::MO_DILATION:
    return dilation(array, ir);
    //
  case MorphologyOperation::MO_EROSION:
    return erosion(array, ir);
    //
  case MorphologyOperation::MO_OPENING:
    return opening(array, ir);
    //
  case MorphologyOperation::MO_BLACK_HAT:
    return morphological_black_hat(array, ir);
    //
  case MorphologyOperation::MO_TOP_HAT:
    return morphological_top_hat(array, ir);
    //
  case MorphologyOperation::MO_GRADIENT:
    return morphological_gradient(array, ir);
    //
  case MorphologyOperation::MO_LAPLACIAN:
    return morphological_laplacian(array, ir);
    //
  case MorphologyOperation::MO_CLOSING_BY_RECONSTRUCTION:
    return closing_by_reconstruction(array, ir);
    //
  case MorphologyOperation::MO_OPENING_BY_RECONSTRUCTION:
    return opening_by_reconstruction(array, ir);
    //
  default: return Array(array.shape);
  }
}

} // namespace hmap

namespace hmap::gpu
{

Array morphological_operators(const Array        &array,
                              int                 ir,
                              MorphologyOperation operation,
                              MinMaxKernel        kernel_type)
{
  if (!validate_non_empty(array)) return Array();

  switch (operation)
  {
  case MorphologyOperation::MO_BORDER:
    return gpu::border(array, ir, kernel_type);
    //
  case MorphologyOperation::MO_CLOSING:
    return gpu::closing(array, ir, kernel_type);
    //
  case MorphologyOperation::MO_DILATION:
    return gpu::dilation(array, ir, kernel_type);
    //
  case MorphologyOperation::MO_EROSION:
    return gpu::erosion(array, ir, kernel_type);
    //
  case MorphologyOperation::MO_OPENING:
    return gpu::opening(array, ir, kernel_type);
    //
  case MorphologyOperation::MO_BLACK_HAT:
    return gpu::morphological_black_hat(array, ir, kernel_type);
    //
  case MorphologyOperation::MO_TOP_HAT:
    return gpu::morphological_top_hat(array, ir, kernel_type);
    //
  case MorphologyOperation::MO_GRADIENT:
    return gpu::morphological_gradient(array, ir, kernel_type);
    //
  case MorphologyOperation::MO_LAPLACIAN:
    return gpu::morphological_laplacian(array, ir, kernel_type);
    //
  case MorphologyOperation::MO_CLOSING_BY_RECONSTRUCTION:
    return gpu::closing_by_reconstruction(array, ir, 0.f, kernel_type);
    //
  case MorphologyOperation::MO_OPENING_BY_RECONSTRUCTION:
    return gpu::opening_by_reconstruction(array, ir, 0.f, kernel_type);
    //
  default: return Array(array.shape);
  }
}

} // namespace hmap::gpu
