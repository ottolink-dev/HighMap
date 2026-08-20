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

  add(
#include "kernels/_common_index.cl"
  );
  add(
#include "kernels/_common_math.cl"
  );
  add(
#include "kernels/_common_rand.cl"
  );
  add(
#include "kernels/_common_sort.cl"
  );
  //
  add(
#include "kernels/advection_particle.cl"
  );
  add(
#include "kernels/advection_warp.cl"
  );
  add(
#include "kernels/bilateral_filter.cl"
  );
  add(
#include "kernels/blend_poisson_bf.cl"
  );
  add(
#include "kernels/coastal_fetch.cl"
  );
  add(
#include "kernels/curvature_quadric.cl"
  );
  add(
#include "kernels/directional_blur.cl"
  );
  add(
#include "kernels/eulerian_transport.cl"
  );
  add(
#include "kernels/expand.cl"
  );
  add(
#include "kernels/flow_accum_stochastic.cl"
  );
  add(
#include "kernels/flow_direction_d8.cl"
  );
  add(
#include "kernels/gabor_wave.cl"
  );
  add(
#include "kernels/gavoronoise.cl"
  );
  add(
#include "kernels/generate_riverbed.cl"
  );
  add(
#include "kernels/gradient_norm.cl"
  );
  add(
#include "kernels/harmonic_interpolation.cl"
  );
  add(
#include "kernels/hemisphere_field.cl"
  );
  add(
#include "kernels/hydraulic_mcdonald.cl"
  );
  add(
#include "kernels/hydraulic_particle.cl"
  );
  add(
#include "kernels/hydraulic_schott.cl"
  );
  add(
#include "kernels/hydraulic_vpipes.cl"
  );
  add(
#include "kernels/interpolate_array.cl"
  );
  add(
#include "kernels/jump_flooding.cl"
  );
  add(
#include "kernels/laplace.cl"
  );
  add(
#include "kernels/laplacian_fract.cl"
  );
  add(
#include "kernels/local_max.cl"
  );
  add(
#include "kernels/local_mean.cl"
  );
  add(
#include "kernels/local_min.cl"
  );
  add(
#include "kernels/local_relief.cl"
  );
  add(
#include "kernels/local_skewness.cl"
  );
  add(
#include "kernels/local_variance.cl"
  );
  add(
#include "kernels/local_z_score.cl"
  );
  add(
#include "kernels/maximum_smooth.cl"
  );
  add(
#include "kernels/mean_shift.cl"
  );
  add(
#include "kernels/median_3x3.cl"
  );
  add(
#include "kernels/minimum_smooth.cl"
  );
  add(
#include "kernels/mountain_range_radial.cl"
  );
  add(
#include "kernels/noise_a.cl"
  );
  add(
#include "kernels/noise_b.cl"
  );
  add(
#include "kernels/normal_displacement.cl"
  );
  add(
#include "kernels/phase_averaging.cl"
  );
  add(
#include "kernels/phase_field.cl"
  );
  add(
#include "kernels/plateau.cl"
  );
  add(
#include "kernels/polygon_field.cl"
  );
  add(
#include "kernels/project_slope_along_direction.cl"
  );
  add(
#include "kernels/rotate.cl"
  );
  add(
#include "kernels/ruggedness.cl"
  );
  add(
#include "kernels/rugosity.cl"
  );
  add(
#include "kernels/sdf_2d_polyline.cl"
  );
  add(
#include "kernels/shallow_viscous_flow.cl"
  );
  add(
#include "kernels/skeleton.cl"
  );
  add(
#include "kernels/smooth_cpulse.cl"
  );
  add(
#include "kernels/snow_simulation.cl"
  );
  add(
#include "kernels/sparse_max_convolution.cl"
  );
  add(
#include "kernels/thermal.cl"
  );
  add(
#include "kernels/thermal_flatten.cl"
  );
  add(
#include "kernels/thermal_inflate.cl"
  );
  add(
#include "kernels/thermal_olsen.cl"
  );
  add(
#include "kernels/thermal_rib.cl"
  );
  add(
#include "kernels/thermal_ridge.cl"
  );
  add(
#include "kernels/thermal_schott.cl"
  );
  add(
#include "kernels/thermal_scree.cl"
  );
  add(
#include "kernels/topographic_position_index.cl"
  );
  add(
#include "kernels/vorolines.cl"
  );
  add(
#include "kernels/voronoi_base.cl"
  );
  add(
#include "kernels/voronoi_edge_distance.cl"
  );
  add(
#include "kernels/voronoi_fbm.cl"
  );
  add(
#include "kernels/voronoi_main.cl"
  );
  add(
#include "kernels/voronoise.cl"
  );
  add(
#include "kernels/vororand_main.cl"
  );
  add(
#include "kernels/warp.cl"
  );
  add(
#include "kernels/water_depth_filter.cl"
  );
  add(
#include "kernels/wavelet_noise.cl"
  );
  //
  add(
#include "kernels/rifts.cl"
  );
  add(
#include "kernels/strata.cl"
  );
  add(
#include "kernels/strata_cells.cl"
  );
  add(
#include "kernels/strata_terrace.cl"
  );

  km.build_program();

  return true;
}

} // namespace hmap::gpu
