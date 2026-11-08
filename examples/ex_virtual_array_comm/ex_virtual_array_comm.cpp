/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file ex_virtual_array_comm.cpp
 * @brief Demonstrator: comparing smooth_cpulse with and without communicator
 * boundary synchronization across tiles in a VirtualArray against a monolithic
 * reference Array.
 */

#include "highmap.hpp"

int main()
{
  const glm::ivec2 shape{512, 512};
  const glm::ivec2 tile_shape{64, 64};
  const int        halo = 16;
  const int        ir = 64;

  const hmap::StorageMode storage = hmap::StorageMode::VA_RAM;
  const hmap::ComputeMode cm{.mode = hmap::ForEachMode::VA_DISTRIBUTED};

  // 1. Generate base elevation noise on monolithic Array
  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   {4.f, 4.f},
                                   /*seed=*/42);

  // 2. Monolithic reference (single Array, standalone)
  hmap::Array ref_smooth = z0;
  hmap::smooth_cpulse(ref_smooth, ir);

  // 3. Tiled VirtualArray without any boundary synchronization
  hmap::VirtualArray va_raw(shape, tile_shape, halo, storage);
  va_raw.from_array(z0, cm);

  hmap::for_each_tile(
      va_raw,
      [&](hmap::Array &tile, const hmap::TileRegion &)
      { hmap::smooth_cpulse(tile, ir); },
      cm);

  hmap::Array img_raw = va_raw.to_array(cm);

  // 4. Tiled VirtualArray with post-processing
  // VirtualArray::sync_overlap_buffers()
  hmap::VirtualArray va_post_sync(shape, tile_shape, halo, storage);
  va_post_sync.from_array(z0, cm);

  hmap::for_each_tile(
      va_post_sync,
      [&](hmap::Array &tile, const hmap::TileRegion &)
      { hmap::smooth_cpulse(tile, ir); },
      cm);

  va_post_sync.sync_overlap_buffers();
  hmap::Array img_post_sync = va_post_sync.to_array(cm);

  // 5. Tiled VirtualArray with communicator boundary synchronization inside
  // loop
  hmap::VirtualArray va_comm_sync(shape, tile_shape, halo, storage);
  va_comm_sync.from_array(z0, cm);

  hmap::for_each_tile(
      va_comm_sync,
      [&](hmap::Array            &tile,
          const hmap::TileRegion &region,
          hmap::ICommunicator    &comm)
      {
        // Demonstrator: smooth_cpulse using communicator to sync overlap
        // buffers
        hmap::smooth_cpulse(tile, ir, nullptr, region, comm);
      },
      cm);

  hmap::Array img_comm_sync = va_comm_sync.to_array(cm);

  // 6. Compute difference maps against the monolithic reference
  hmap::Array diff_raw = hmap::abs(img_raw - ref_smooth);
  hmap::Array diff_post_sync = hmap::abs(img_post_sync - ref_smooth);
  hmap::Array diff_comm_sync = hmap::abs(img_comm_sync - ref_smooth);
  hmap::Array diff_post_vs_comm = hmap::abs(img_post_sync - img_comm_sync);

  hmap::log::info("Max diff raw (no sync) vs ref: {}", diff_raw.max());
  hmap::log::info("Max diff VA post-sync vs ref:  {}", diff_post_sync.max());
  hmap::log::info("Max diff comm sync vs ref:     {}", diff_comm_sync.max());
  hmap::log::info("Max diff post-sync vs comm:    {}", diff_post_vs_comm.max());

  // 7. Export comparative banner
  hmap::export_banner_png("ex_virtual_array_comm.png",
                          {ref_smooth, img_raw, img_post_sync, img_comm_sync},
                          hmap::Cmap::TERRAIN,
                          /*draw_contour=*/false,
                          /*normalize=*/true);

  hmap::export_banner_png(
      "ex_virtual_array_comm_diff.png",
      {diff_raw, diff_post_sync, diff_comm_sync, diff_post_vs_comm},
      hmap::Cmap::VIRIDIS,
      /*draw_contour=*/false,
      /*normalize=*/true);

  return 0;
}
