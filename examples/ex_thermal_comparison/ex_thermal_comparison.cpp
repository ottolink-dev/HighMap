#include <iostream>

#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();
  hmap::init_openmp();

  const glm::ivec2 shape = {512, 512};
  const glm::vec2  kw = {3.f, 3.f};
  const int        seed = 42;

  // --- Initial terrain setup: rough mountain cone with sharp features

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   kw,
                                   seed,
                                   8,
                                   0.f);
  z0 = hmap::bulkify(z0, hmap::PrimitiveType::PRIM_CONE_SMOOTH, 0.9f);
  hmap::remap(z0);

  const float talus_val = 2.5f / static_cast<float>(shape.x);
  const hmap::Array talus = hmap::constant(shape, talus_val);
  const int   iterations = 30;

  // 1. Standard Thermal (Musgrave / Shadertoy exchange)
  auto z_std = z0;
  hmap::gpu::thermal(z_std, talus, iterations);

  // 2. Mass-Conserving Thermal
  auto z_conserve = z0;
  hmap::gpu::thermal_conserve(z_conserve, talus, iterations, 0.5f);

  // 3. Schott Model (Uniform threshold & neighbor count balance)
  auto z_schott = z0;
  hmap::gpu::thermal_schott(z_schott, talus, iterations, 0.2f);

  // 4. Olsen Model (Proportional redistribution with strict talus limits)
  auto z_olsen = z0;
  hmap::gpu::thermal_olsen(z_olsen, talus, iterations);

  // 5. Thermal Ridge (Sharpened crests and talus apron)
  auto z_ridge = z0;
  hmap::gpu::thermal_ridge(z_ridge, talus, iterations);

  // 6. Thermal Flatten (Slope-selective asymmetric relaxation)
  auto z_flatten = z0;
  hmap::gpu::thermal_flatten(z_flatten, talus, iterations, 0.5f, 0.0f);

  // --- Export individual images with hillshading
  z0.to_png("docs/images/thermal_comparison_original.png",
            hmap::Cmap::TERRAIN,
            true);
  z_std.to_png("docs/images/thermal_comparison_standard.png",
               hmap::Cmap::TERRAIN,
               true);
  z_conserve.to_png("docs/images/thermal_comparison_conserve.png",
                    hmap::Cmap::TERRAIN,
                    true);
  z_schott.to_png("docs/images/thermal_comparison_schott.png",
                  hmap::Cmap::TERRAIN,
                  true);
  z_olsen.to_png("docs/images/thermal_comparison_olsen.png",
                 hmap::Cmap::TERRAIN,
                 true);
  z_ridge.to_png("docs/images/thermal_comparison_ridge.png",
                 hmap::Cmap::TERRAIN,
                 true);
  z_flatten.to_png("docs/images/thermal_comparison_flatten.png",
                   hmap::Cmap::TERRAIN,
                   true);

  // --- Export composite comparison banners
  hmap::export_banner_png("docs/images/thermal_comparison_banner1.png",
                          {z0, z_std, z_conserve, z_schott},
                          hmap::Cmap::TERRAIN,
                          true);

  hmap::export_banner_png("docs/images/thermal_comparison_banner2.png",
                          {z0, z_olsen, z_ridge, z_flatten},
                          hmap::Cmap::TERRAIN,
                          true);

  std::cout << "Thermal comparison images successfully generated.\n";
  return 0;
}
