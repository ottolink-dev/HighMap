/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>
#include <vector>

#include "cl_wrapper/device_manager.hpp"
#include "cl_wrapper/kernel_manager.hpp"
#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/opencl/gpu_opencl.hpp"

namespace hmap::gpu
{

void helper_bind_optional_buffer(clwrapper::Run    &run,
                                 const std::string &id,
                                 const Array       *p_array)
{
  std::vector<float> dummy_vector(1);

  if (p_array)
  {
    run.bind_buffer<float>(id, p_array->vector);
    run.write_buffer(id);
  }
  else
    run.bind_buffer<float>(id, dummy_vector);
}

bool init_opencl()
{
  if (!clwrapper::DeviceManager::get_instance().is_ready()) return false;

  auto &km = clwrapper::KernelManager::get_instance();
  km.clear_sources();

  std::string opencl_build_options = "-cl-fast-relaxed-math "
                                     "-cl-mad-enable "
                                     "-cl-no-signed-zeros "
                                     "-cl-denorms-are-zero "
                                     "-cl-finite-math-only ";

  km.set_build_options(opencl_build_options);

  // load and build kernels
  auto add = [&](const std::string &src) { km.add_kernel(src, false, false); };

  // --- Common utilities

  add(
#include "kernels/common/_common_index.cl"
  );
  add(
#include "kernels/common/_common_math.cl"
  );
  add(
#include "kernels/common/_common_rand.cl"
  );
  add(
#include "kernels/common/_common_sort.cl"
  );
  add(
#include "kernels/common/_common_d8.cl"
  );

  // --- Noise and procedural primitives

  add(
#include "kernels/noise/noise_a.cl"
  );
  add(
#include "kernels/noise/noise_b.cl"
  );
  add(
#include "kernels/noise/gabor_wave.cl"
  );
  add(
#include "kernels/noise/gavoronoise.cl"
  );
  add(
#include "kernels/noise/hemisphere_field.cl"
  );
  add(
#include "kernels/noise/laplacian_fract.cl"
  );
  add(
#include "kernels/noise/mountain_range_radial.cl"
  );
  add(
#include "kernels/noise/phase_averaging.cl"
  );
  add(
#include "kernels/noise/phase_field.cl"
  );
  add(
#include "kernels/noise/polygon_field.cl"
  );
  add(
#include "kernels/noise/wavelet_noise.cl"
  );

  // --- Voronoi diagrams and patterns

  add(
#include "kernels/voronoi/voronoi_base.cl"
  );
  add(
#include "kernels/voronoi/vorolines.cl"
  );
  add(
#include "kernels/voronoi/voronoi_edge_distance.cl"
  );
  add(
#include "kernels/voronoi/voronoi_fbm.cl"
  );
  add(
#include "kernels/voronoi/voronoi_main.cl"
  );
  add(
#include "kernels/voronoi/voronoise.cl"
  );
  add(
#include "kernels/voronoi/vororand_main.cl"
  );

  // --- Erosion and geomorphology

  add(
#include "kernels/erosion/hydraulic_mcdonald.cl"
  );
  add(
#include "kernels/erosion/hydraulic_musgrave.cl"
  );
  add(
#include "kernels/erosion/hydraulic_particle.cl"
  );
  add(
#include "kernels/erosion/hydraulic_schott.cl"
  );
  add(
#include "kernels/erosion/hydraulic_vpipes.cl"
  );
  add(
#include "kernels/erosion/rifts.cl"
  );
  add(
#include "kernels/erosion/strata.cl"
  );
  add(
#include "kernels/erosion/strata_cells.cl"
  );
  add(
#include "kernels/erosion/strata_terrace.cl"
  );
  add(
#include "kernels/erosion/thermal.cl"
  );
  add(
#include "kernels/erosion/thermal_flatten.cl"
  );
  add(
#include "kernels/erosion/thermal_inflate.cl"
  );
  add(
#include "kernels/erosion/thermal_olsen.cl"
  );
  add(
#include "kernels/erosion/thermal_rib.cl"
  );
  add(
#include "kernels/erosion/thermal_ridge.cl"
  );
  add(
#include "kernels/erosion/thermal_schott.cl"
  );
  add(
#include "kernels/erosion/thermal_scree.cl"
  );

  // --- Filters

  add(
#include "kernels/filters/bilateral_filter.cl"
  );
  add(
#include "kernels/filters/directional_blur.cl"
  );
  add(
#include "kernels/filters/expand.cl"
  );
  add(
#include "kernels/filters/laplace.cl"
  );
  add(
#include "kernels/filters/mean_shift.cl"
  );
  add(
#include "kernels/filters/median_3x3.cl"
  );
  add(
#include "kernels/filters/normal_displacement.cl"
  );
  add(
#include "kernels/filters/plateau.cl"
  );
  add(
#include "kernels/filters/smooth_cpulse.cl"
  );
  add(
#include "kernels/filters/sparse_max_convolution.cl"
  );
  add(
#include "kernels/filters/water_depth_filter.cl"
  );

  // --- Flow and hydrology

  add(
#include "kernels/flow/coastal_fetch.cl"
  );
  add(
#include "kernels/flow/flow_accum_stochastic.cl"
  );
  add(
#include "kernels/flow/flow_direction_d8.cl"
  );
  add(
#include "kernels/flow/generate_riverbed.cl"
  );
  add(
#include "kernels/flow/shallow_viscous_flow.cl"
  );
  add(
#include "kernels/flow/snow_simulation.cl"
  );

  // --- Local metrics

  add(
#include "kernels/local_metrics/curvature_quadric.cl"
  );
  add(
#include "kernels/local_metrics/local_max.cl"
  );
  add(
#include "kernels/local_metrics/local_max_octagon.cl"
  );
  add(
#include "kernels/local_metrics/local_max_square.cl"
  );
  add(
#include "kernels/local_metrics/local_mean.cl"
  );
  add(
#include "kernels/local_metrics/local_min.cl"
  );
  add(
#include "kernels/local_metrics/local_min_octagon.cl"
  );
  add(
#include "kernels/local_metrics/local_min_square.cl"
  );
  add(
#include "kernels/local_metrics/local_relief.cl"
  );
  add(
#include "kernels/local_metrics/local_skewness.cl"
  );
  add(
#include "kernels/local_metrics/local_variance.cl"
  );
  add(
#include "kernels/local_metrics/local_z_score.cl"
  );
  add(
#include "kernels/local_metrics/ridge_accentuate.cl"
  );
  add(
#include "kernels/local_metrics/ruggedness.cl"
  );
  add(
#include "kernels/local_metrics/rugosity.cl"
  );
  add(
#include "kernels/local_metrics/topographic_position_index.cl"
  );

  // --- Transformations, transport, and interpolation

  add(
#include "kernels/transform/advection_particle.cl"
  );
  add(
#include "kernels/transform/advection_warp.cl"
  );
  add(
#include "kernels/transform/blend_poisson_bf.cl"
  );
  add(
#include "kernels/transform/eulerian_transport.cl"
  );
  add(
#include "kernels/transform/harmonic_interpolation.cl"
  );
  add(
#include "kernels/transform/interpolate_array.cl"
  );
  add(
#include "kernels/transform/jagged.cl"
  );
  add(
#include "kernels/transform/jump_flooding.cl"
  );
  add(
#include "kernels/transform/project_slope_along_direction.cl"
  );
  add(
#include "kernels/transform/rotate.cl"
  );
  add(
#include "kernels/transform/sdf_2d_polyline.cl"
  );
  add(
#include "kernels/transform/skeleton.cl"
  );
  add(
#include "kernels/transform/voronoi_shrink.cl"
  );
  add(
#include "kernels/transform/warp.cl"
  );
  add(
#include "kernels/rocks/rock_simulation.cl"
  );

  km.build_program();

  return true;
}

} // namespace hmap::gpu
