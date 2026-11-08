/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <memory>
#include <string>

#include "highmap/logger.hpp"
#include "highmap/virtual_array/tile_storage.hpp"

namespace hmap
{

void TileStorage::trim()
{
}

std::unique_ptr<TileStorage> make_storage(glm::ivec2 /* shape */,
                                          glm::ivec2 /* tile_shape */,
                                          StorageMode storage_mode)
{
  switch (storage_mode)
  {
  case StorageMode::VA_RAM: return std::make_unique<RamTileStorage>();

  case StorageMode::VA_DISK_LRU_STREAM:
  case StorageMode::VA_DISK_LRU_FULL:
    return std::make_unique<DiskLruTileStorage>(storage_mode);
  }

  // fallback
  hmap::log::error("Unknown StorageMode, defaulting to RAM");
  return std::make_unique<RamTileStorage>();
}

std::string to_string(StorageMode m)
{
  switch (m)
  {
  case StorageMode::VA_RAM: return "VA_RAM";
  case StorageMode::VA_DISK_LRU_STREAM: return "VA_DISK_LRU_STREAM";
  case StorageMode::VA_DISK_LRU_FULL: return "VA_DISK_LRU_FULL";
  }
  return "UNKNOWN";
}

} // namespace hmap
