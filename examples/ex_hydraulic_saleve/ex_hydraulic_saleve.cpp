#include "highmap.hpp"
#include "highmap/dbg/timer.hpp"

#include <omp.h>

int main(void)
{
  hmap::gpu::init_opencl();
  hmap::init_openmp();

  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 0;

  // --- erosion of an heightmap

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   kw,
                                   seed,
                                   8,
                                   0.7f);
  z0 = hmap::bulkify(z0, hmap::PrimitiveType::PRIM_CONE_SMOOTH, 1.f);

  size_t control_points_count = 10000;

  float strength = 0.5f;

  float m_exp = 0.8f;
  float uplift_rate = 1.f;
  float tolerance = 1e-3f;
  int   max_iterations = 1000;

  float smin = 0.f;
  float smax = 6.f;
  bool  scale_erodibility_with_z = true;
  float erodibility_distrib_exp = 1.f;
  float noise_strength = 0.25f; // stream tree structural noise

  hmap::Array dx = 0.03f * hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                           shape,
                                           kw,
                                           ++seed,
                                           8,
                                           0.f);
  hmap::Array dy = 0.03f * hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                           shape,
                                           kw,
                                           ++seed,
                                           8,
                                           0.f);

  hmap::Timer::Start("hydraulic_saleve");
  hmap::Array z1 = hmap::hydraulic_saleve(
      z0,
      ++seed,
      control_points_count,
      m_exp,
      uplift_rate,
      tolerance,
      max_iterations,
      smin,
      smax,
      strength,
      scale_erodibility_with_z,
      erodibility_distrib_exp,
      noise_strength,
      /* enable_post_slope_limiter */ false,
      /* post_slope_limit */ 0.f,
      /* enable_post_smoothing */ false,
      hmap::InterpolationMethod2D::ITP2D_DELAUNAY_GRADIENT,
      &dx,
      &dy);
  hmap::Timer::Stop("hydraulic_saleve");

  // mimic deposition
  int   deposition_ir = int(0.1f * shape.x);
  float deposition_strength = 0.5f;
  int   iterations = 2;
  hmap::gpu::deposition_fill_holes(z1,
                                   deposition_ir,
                                   deposition_strength,
                                   iterations);

  hmap::Timer::Dump();

  z1.dump();

  // --- export

  hmap::export_banner_png("ex_hydraulic_saleve.png",
                          {z0, z1},
                          hmap::Cmap::TERRAIN,
                          true);
}
