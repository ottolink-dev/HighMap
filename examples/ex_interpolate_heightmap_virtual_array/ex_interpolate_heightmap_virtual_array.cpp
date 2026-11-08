#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 256};
  glm::ivec2 tile_shape = {64, 64};
  int        halo = 16;

  // --- virtual array params

  auto              storage_mode = hmap::StorageMode::VA_DISK_LRU_FULL;
  hmap::ComputeMode cm = {.mode = hmap::ForEachMode::VA_DISTRIBUTED,
                          .trim_storage = false};
  glm::vec4 bbox(0.f, 1.0f, 0.f, 1.f); // local, within the virtual array...

  // --- data

  hmap::VirtualArray varray1(shape, bbox, tile_shape, halo, storage_mode);
  hmap::VirtualArray varray2(shape, bbox, tile_shape, halo, storage_mode);

  auto lambda_noise = [](hmap::Array &tile, const hmap::TileRegion &region)
  {
    tile = hmap::noise(hmap::NoiseType::PERLIN,
                       region.shape,
                       {2.f, 2.f},
                       0,
                       nullptr,
                       nullptr,
                       region.bbox);
  };

  hmap::for_each_tile(varray1, lambda_noise, cm);

  // --- coord frames (global, not related to virtual array bbox)

  hmap::CoordFrame frame1 = hmap::CoordFrame(glm::vec2(10.f, 20.f),
                                             glm::vec2(50.f, 100.f),
                                             30.f);

  hmap::CoordFrame frame2 = hmap::CoordFrame(glm::vec2(-20.f, 50.f),
                                             glm::vec2(100.f, 70.f),
                                             -30.f);

  // --- interpolate

  hmap::interpolate_heightmap(varray1, varray2, frame1, frame2, cm);

  varray2.to_array(cm).to_png("out.png", hmap::Cmap::JET);
}
