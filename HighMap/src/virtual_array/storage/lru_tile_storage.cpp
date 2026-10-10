/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <list>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/tile_storage.hpp"

#include <unordered_map>

namespace hmap
{

LruTileStorage::LruTileStorage(size_t max_tiles) : max_tiles(max_tiles)
{
}

std::unique_ptr<TileStorage> LruTileStorage::clone() const
{
  return std::make_unique<LruTileStorage>(this->max_tiles);
}

Array &LruTileStorage::get_tile(const TileRegion &region)
{
  std::lock_guard<std::mutex> lock(mutex);
  return this->get_tile_no_mutex_lock(region);
}

Array &LruTileStorage::get_tile_no_mutex_lock(const TileRegion &region)
{
  const TileKey &key = region.key;

  // ---- hit

  auto it = tiles.find(key);
  if (it != tiles.end())
  {
    this->lru.splice(this->lru.begin(), lru, it->second.lru_it);
    return it->second.value;
  }

  // ---- miss → allocate

  glm::ivec2 total = region.shape;
  Array      tile(total);

  // ---- eviction if needed

  if (this->tiles.size() >= this->max_tiles)
  {
    const TileKey &evict_key = this->lru.back();
    auto           evict_it = this->tiles.find(evict_key);

    this->on_evict(evict_key, evict_it->second.value);

    this->tiles.erase(evict_it);
    this->lru.pop_back();
  }

  // ---- insert new tile

  this->lru.push_front(key);
  auto [inserted_it, _] = this->tiles.emplace(
      key,
      LruTileEntry{std::move(tile), this->lru.begin()});

  return inserted_it->second.value;
}

std::string LruTileStorage::info_string() const
{
  return "LRU (capacity=" + std::to_string(this->max_live_tiles()) + ")";
}

size_t LruTileStorage::max_live_tiles() const
{
  return this->max_tiles;
}

size_t LruTileStorage::live_tile_count() const
{
  std::lock_guard<std::mutex> lock(mutex);
  return tiles.size();
}

size_t LruTileStorage::live_memory_bytes() const
{
  std::lock_guard<std::mutex> lock(mutex);
  size_t                      bytes = 0;
  for (const auto &[key, entry] : tiles)
  {
    bytes += entry.value.vector.capacity() * sizeof(float) +
             sizeof(LruTileEntry) + sizeof(key);
  }
  return bytes;
}

void LruTileStorage::release_tile(const TileRegion & /* region */)
{
  // no-op
}

void LruTileStorage::on_evict(const TileKey & /* key */, Array & /* tile */)
{
  // default: do nothing
}

} // namespace hmap
