#include <iostream>

#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  // shape = {1024, 1024};
  // shape = {512, 512};
  std::uint32_t seed = 0;

  constexpr int       steps = 500;
  constexpr float     alpha = 0.f / 180.f * M_PI;
  constexpr glm::vec2 wind_dir = {std::cos(alpha), std::sin(alpha)};

  // --- 1. Transverse sand dunes (unimodal wind with abundant sand supply)

  hmap::Array z_transverse(shape, 0.05f);

  hmap::aeolian_erosion(z_transverse,
                        steps,
                        /* seed */ seed,
                        3.f * wind_dir,
                        /* p_bedrock */ nullptr,
                        /* p_vegetation */ nullptr,
                        /* p_abrasion_mask */ nullptr,
                        /* p_sediment_map */ nullptr,
                        /* p_bedrock_map */ nullptr,
                        /* cell_size */ -1.f,
                        /* matter_to_move */ 0.001f,
                        /* samples_per_step */ shape.x * shape.y,
                        /* sand_repose_angle */ 33.f,
                        /* shadow_angle_min */ 5.f,
                        /* shadow_angle_max */ 15.f,
                        /* max_shadow_dist */ 12.f,
                        /* abrasion_rate */ 0.f,
                        /* repose_rate */ 0.35f);

  z_transverse.dump("zt.png");

  // --- 2. Transverse dunes with bedrock bumps (sand flowing around and piling)

  hmap::Array bed_bumps = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                          shape,
                                          {3.f, 3.f},
                                          seed + 10u,
                                          8,
                                          0.f);
  hmap::remap(bed_bumps, 0.f, 0.5f);
  hmap::clamp_min(bed_bumps, 0.5f * bed_bumps.max());
  bed_bumps -= 0.5f * bed_bumps.max();

  // Initial terrain: bedrock bumps + abundant sand layer (0.05 units)
  const float smax = 0.01f;

  hmap::Array z_bumps = bed_bumps + hmap::Array(shape, smax);

  hmap::aeolian_erosion(z_bumps,
                        steps,
                        /* seed */ seed,
                        3.f * wind_dir,
                        /* p_bedrock */ &bed_bumps,
                        /* p_vegetation */ nullptr,
                        /* p_abrasion_mask */ nullptr,
                        /* p_sediment_map */ nullptr,
                        /* p_bedrock_map */ nullptr,
                        /* cell_size */ -1.f,
                        /* matter_to_move */ smax / 50.f,
                        /* samples_per_step */ shape.x * shape.y,
                        /* sand_repose_angle */ 33.f,
                        /* shadow_angle_min */ 5.f,
                        /* shadow_angle_max */ 15.f,
                        /* max_shadow_dist */ 12.f,
                        /* abrasion_rate */ 0.f,
                        /* repose_rate */ 0.35f);

  bed_bumps.dump("z.png");
  (z_bumps - bed_bumps).dump("s.png");

  // --- 3. Barchan dunes (higher wind speed with lower sand supply)

  hmap::Array z_barchan(shape, 0.0125f);

  hmap::aeolian_erosion(z_barchan,
                        steps,
                        /* seed */ seed,
                        5.f * wind_dir,
                        /* p_bedrock */ nullptr,
                        /* p_vegetation */ nullptr,
                        /* p_abrasion_mask */ nullptr,
                        /* p_sediment_map */ nullptr,
                        /* p_bedrock_map */ nullptr,
                        /* cell_size */ -1.f,
                        /* matter_to_move */ 0.001f,
                        /* samples_per_step */ shape.x * shape.y,
                        /* sand_repose_angle */ 33.f,
                        /* shadow_angle_min */ 5.f,
                        /* shadow_angle_max */ 15.f,
                        /* max_shadow_dist */ 12.f,
                        /* abrasion_rate */ 0.f,
                        /* repose_rate */ 0.35f);

  // --- 4. Yardangs (wind abrasion carving corrugated bedrock corridors)

  hmap::Array bed_yardang(shape, 0.05f);
  hmap::Array z_yardang = bed_yardang + hmap::Array(shape, 0.005f);

  hmap::Array abrasion_mask = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                              shape,
                                              {4.f, 16.f},
                                              seed);
  hmap::remap(abrasion_mask, 0.2f, 1.f);

  hmap::aeolian_erosion(z_yardang,
                        steps,
                        /* seed */ seed,
                        3.f * wind_dir,
                        /* p_bedrock */ &bed_yardang,
                        /* p_vegetation */ nullptr,
                        /* p_abrasion_mask */ &abrasion_mask,
                        /* p_sediment_map */ nullptr,
                        /* p_bedrock_map */ nullptr,
                        /* cell_size */ -1.f,
                        /* matter_to_move */ 0.001f,
                        /* samples_per_step */ shape.x * shape.y,
                        /* sand_repose_angle */ 33.f,
                        /* shadow_angle_min */ 5.f,
                        /* shadow_angle_max */ 15.f,
                        /* max_shadow_dist */ 12.f,
                        /* abrasion_rate */ 0.0005f,
                        /* repose_rate */ 0.35f);

  // --- 5. Export banner

  hmap::export_banner_png("ex_aeolian_erosion.png",
                          {z_transverse, z_bumps, z_barchan, z_yardang},
                          hmap::Cmap::TERRAIN,
                          true,
                          true);

  std::cout << "ex_aeolian_erosion successfully executed." << std::endl;
  return 0;
}
