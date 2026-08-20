/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <limits>
#include <memory>
#include <mutex>
#include <utility>

#include "highmap/array.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/tile_storage.hpp"

#include <unordered_map>

namespace hmap
{

std::unique_ptr<TileStorage> RamTileStorage::clone() const
{
  // std::mutex is not copyable: copy the tiles manually under the lock
  auto clone = std::make_unique<RamTileStorage>();

  const std::lock_guard<std::mutex> lock(this->mutex);
  clone->tiles = this->tiles;

  return clone;
}

Array &RamTileStorage::get_tile(const TileRegion &region)
{
  // concurrent workers (VA_DISTRIBUTED mode) may request/insert tiles at
  // the same time: guard the lazy creation path
  const std::lock_guard<std::mutex> lock(this->mutex);

  auto it = tiles.find(region.key);
  if (it != tiles.end()) return it->second;

  glm::ivec2 total = region.shape; // include halo
  Array      tile(total);
  auto [inserted_it, _] = tiles.emplace(region.key, std::move(tile));
  return inserted_it->second;
}

size_t RamTileStorage::max_live_tiles() const
{
  return std::numeric_limits<size_t>::max();
}

void RamTileStorage::release_tile(const TileRegion & /* region */)
{
  // nothing for RAM
}

} // namespace hmap
