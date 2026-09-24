#include <stdexcept>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"

namespace hmap::gpu
{

Array local_metrics(const Array &array,
                    int          ir,
                    LocalMetrics metric,
                    MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();

  // clang-format off
  switch (metric)
  {
  case LocalMetrics::LM_LOCAL_ASPECT_VARIANCE:
    return gpu::local_aspect_variance(array, ir);
    //
  case LocalMetrics::LM_LOCAL_MAX:
    return gpu::local_max(array, ir, kernel_type);
    //
  case LocalMetrics::LM_LOCAL_MEDIAN_DEVIATION:
    return gpu::local_median_deviation(array, ir);
    //
  case LocalMetrics::LM_LOCAL_MIN:
    return gpu::local_min(array, ir, kernel_type);
    //
  case LocalMetrics::LM_LOCAL_RELIEF:
    return gpu::local_relief(array, ir, kernel_type);
    //
  case LocalMetrics::LM_LOCAL_VARIANCE:
    return gpu::local_variance(array, ir);
    //
  case LocalMetrics::LM_LOCAL_MEAN:
    return gpu::local_mean(array, ir);
    //
  case LocalMetrics::LM_LOCAL_SKEWNESS:
    return gpu::local_skewness(array, ir);
    //
  case LocalMetrics::LM_LOCAL_Z_SCORE:
    return gpu::local_z_score(array, ir);
    //
  case LocalMetrics::LM_TOPOGRAPHIC_POSITION_INDEX:
    return gpu::topographic_position_index(array, ir);
    //
  case LocalMetrics::LM_RELATIVE_ELEVATION:
    return gpu::relative_elevation(array, ir, kernel_type);
    //
  case LocalMetrics::LM_RUGGEDNESS:
    return gpu::ruggedness(array, ir);
    //
  case LocalMetrics::LM_RUGOSITY_CONCAVE:
    return gpu::rugosity(array, ir, /* convex */ false);
    //
  case LocalMetrics::LM_RUGOSITY_CONVEX:
    return gpu::rugosity(array, ir, /* convex */ true);
    //
  default: throw std::runtime_error("unknown LocalMetrics in local_metrics");
  }
  // clang-format on
}

} // namespace hmap::gpu
