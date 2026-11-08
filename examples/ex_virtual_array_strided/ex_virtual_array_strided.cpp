#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "highmap.hpp"

int main()
{
  hmap::gpu::init_opencl();

  // ===========================================================================
  // Configuration
  // ===========================================================================

  const glm::vec4 bbox{0.f, 1.f, 0.f, 1.f};

  const glm::ivec2 shape{1024, 1024};
  const glm::ivec2 tile_shape{220, 312};
  const size_t     halo = 32;

  // const auto storage_mode = hmap::StorageMode::VA_RAM;
  const auto storage_mode = hmap::StorageMode::VA_DISK_LRU_FULL;
  // const auto storage_mode = hmap::StorageMode::VA_DISK_LRU_STREAM;
  // const auto storage_mode = hmap::StorageMode::VA_DISK_LRU_STREAM;

  // const auto foreach_mode = hmap::ForEachMode::VA_DISTRIBUTED;
  // const auto foreach_mode = hmap::ForEachMode::VA_SEQUENTIAL;
  const auto foreach_mode = hmap::ForEachMode::VA_SINGLE_ARRAY;

  const hmap::ComputeMode compute_mode{.mode = foreach_mode,
                                       .trim_storage = true,
                                       .stride = 1};

  // ===========================================================================
  // Virtual arrays
  // ===========================================================================

  hmap::VirtualArray varray(shape, bbox, tile_shape, halo, storage_mode);

  // ===========================================================================
  // Tile operations
  // ===========================================================================

  const auto generate_noise =
      [](hmap::Array &tile, const hmap::TileRegion &region)
  {
    tile = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                           region.shape,
                           {2.f, 2.f},
                           0,
                           8,
                           0.7f,
                           0.5f,
                           2.f,
                           nullptr,
                           nullptr,
                           nullptr,
                           region.bbox);
  };

  const auto apply_gamma = [](hmap::Array &tile, const hmap::TileRegion &)
  { hmap::gamma_correction(tile, 2.f); };

  hmap::for_each_tile(varray, generate_noise, compute_mode);
  hmap::for_each_tile(varray, apply_gamma, compute_mode);

  // specific
  {
    const auto erosion = [](hmap::Array &tile, const hmap::TileRegion &)
    { hmap::gpu::hydraulic_stream_log(tile, 0.2f, 1e0f); };

    hmap::ComputeMode cm_dwn = {.mode = hmap::ForEachMode::VA_SINGLE_ARRAY,
                                .stride = 8};

    hmap::for_each_tile(varray, erosion, cm_dwn);
  }

  // visual output
  {
    const hmap::ComputeMode cm_local = {
        .mode = hmap::ForEachMode::VA_SEQUENTIAL,
        .trim_storage = false};
    varray.to_array(cm_local).to_png("out_noise.png", hmap::Cmap::JET);
  }

  return 0;
}
