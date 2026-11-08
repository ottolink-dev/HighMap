#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  const glm::ivec2 shape = {256, 256};
  const glm::vec2  kw = {4.f, 4.f};
  const int        seed = 2;

  // Generate synthetic mountain terrain
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_CUBIC_PULSE, 1.f);
  hmap::remap(z);

  const float strength = 0.05f;
  const float talus = 1.f / float(shape.x);

  // Base snow mantle (altitude, slope exclusion, peaks)
  hmap::Array snow_base = hmap::gpu::snow_mantle(z,
                                                 strength, // strength
                                                 0.45f,    // snowline
                                                 0.10f,    // falloff
                                                 0.80f,    // talus_coverage
                                                 8,
                                                 talus, // talus_threshold
                                                 0.00f, // directional
                                                 0.00f, // wind_direction
                                                 0.30f, // peaks
                                                 0.00f, // melt
                                                 180.f, // sun_azimuth
                                                 45.f,  // sun_zenith
                                                 100,   // thermal_iterations
                                                 talus, // thermal_talus
                                                 0.50f  // coverage
  );

  // Windward vs Leeward directional bias (wind from West -> azimuth 270 deg)
  hmap::Array snow_wind = hmap::gpu::snow_mantle(z,
                                                 strength,
                                                 0.45f,
                                                 0.10f,
                                                 0.80f,
                                                 8,
                                                 talus,
                                                 1.00f, // directional bias
                                                 270.f, // wind azimuth
                                                 0.30f,
                                                 0.00f,
                                                 180.f,
                                                 45.f,
                                                 0,
                                                 talus,
                                                 0.50f);

  // Optional high-frequency texture noise
  hmap::Array noise = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      {8.f, 8.f},
                                      seed);
  hmap::remap(noise, 0.2f, 1.f);

  // Full mantle with solar melt, texture noise, and downhill GPU thermal flow
  hmap::Array snow_full = hmap::gpu::snow_mantle(
      z,
      strength,
      0.45f,
      0.10f,
      0.80f,
      8,
      talus,
      0.40f,
      270.f,
      0.30f,
      0.50f, // melt
      180.f,
      45.f,
      100,   // thermal_iterations (GPU thermal_conserve)
      talus, // thermal_talus
      0.50f,
      nullptr,
      nullptr,
      &noise);

  z.dump("z.png");
  snow_full.dump("snow_depth.png");
  (z + snow_full).dump("z_snow.png");

  hmap::export_banner_png("ex_snow_mantle.png",
                          {z, snow_base, snow_wind, z + snow_full},
                          hmap::Cmap::JET);

  return 0;
}
