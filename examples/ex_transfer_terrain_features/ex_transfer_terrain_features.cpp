#include "highmap.hpp"

int main(void)
{
  // --- Parameters
  const glm::ivec2 shape = {256, 256};
  const glm::ivec2 tile_shape = {128, 128};
  const int        halo = 32;
  const int        seed = 1;

  const int           ir = 32;
  const float         amplitude = 2.f;
  const bool          target_prefiltering = false;
  const int           factor = 8;
  const int           patch_size = 8;
  const int           analysis_stride = 4;
  const int           synthesis_stride = 4;
  const int           sparsity = 3;
  const int           max_atoms = 1024;
  const std::uint32_t dict_seed = 0;

  const std::vector<float> source_weights = {1.f, 0.6f};

  // Coarse target terrain (large-scale structure to preserve)
  hmap::Array target_base = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                            shape,
                                            {2.f, 2.f},
                                            seed,
                                            2);
  hmap::remap(target_base, 0.f, 0.8f);

  // Two detailed source exemplars with different physical scales
  hmap::Array source1 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                        shape,
                                        {8.f, 8.f},
                                        seed + 1,
                                        6);
  hmap::remap(source1, 0.f, 1.f);

  hmap::Array source2 = hmap::noise_ridged(hmap::NoiseType::WORLEY,
                                           {128, 128},
                                           {8.f, 8.f},
                                           seed + 2,
                                           6);

  hmap::remap(source2, 0.f, 1.f);
  hmap::CoordFrame frame_target({0.f, 0.f}, {1.f, 1.f}, 0.f);
  hmap::CoordFrame frame_src1({2.f, 2.f}, {1.f, 1.f}, 0.f);
  hmap::CoordFrame frame_src2({-3.f, 1.f}, {0.5f, 0.5f}, 0.f);

  // Single-source transfer
  hmap::Array out_single = target_base;
  hmap::transfer_terrain_features(out_single,
                                  frame_target,
                                  source1,
                                  frame_src1,
                                  ir,
                                  amplitude,
                                  target_prefiltering,
                                  factor,
                                  patch_size,
                                  analysis_stride,
                                  synthesis_stride,
                                  sparsity,
                                  max_atoms,
                                  dict_seed);

  // Multi-source transfer with custom weights
  hmap::Array out_multi = target_base;
  hmap::transfer_terrain_features(out_multi,
                                  frame_target,
                                  {&source1, &source2},
                                  {&frame_src1, &frame_src2},
                                  ir,
                                  amplitude,
                                  target_prefiltering,
                                  source_weights,
                                  factor,
                                  patch_size,
                                  analysis_stride,
                                  synthesis_stride,
                                  sparsity,
                                  max_atoms,
                                  dict_seed);

  // Tiled VirtualArray transfer
  hmap::VirtualArray va_target(shape,
                               tile_shape,
                               halo,
                               hmap::StorageMode::VA_RAM);
  va_target.from_array(target_base, hmap::ComputeMode{});

  hmap::VirtualArray va_src1(shape, tile_shape, 0, hmap::StorageMode::VA_RAM);
  va_src1.from_array(source1, hmap::ComputeMode{});

  hmap::transfer_terrain_features(
      va_target,
      frame_target,
      va_src1,
      frame_src1,
      ir,
      amplitude,
      target_prefiltering,
      factor,
      patch_size,
      analysis_stride,
      synthesis_stride,
      sparsity,
      max_atoms,
      dict_seed,
      hmap::ComputeMode{.mode = hmap::ForEachMode::VA_DISTRIBUTED});

  va_target.sync_overlap_buffers();

  hmap::Array out_va = va_target.to_array(hmap::ComputeMode{});

  hmap::export_banner_png("ex_transfer_terrain_features.png",
                          {target_base, source1, out_single, out_multi, out_va},
                          hmap::Cmap::TERRAIN);
}
