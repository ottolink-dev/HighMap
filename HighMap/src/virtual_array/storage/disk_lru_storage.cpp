/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <list>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/string_utils.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/tile_storage.hpp"

#include <unordered_map>

namespace hmap
{

DiskLruTileStorage::DiskLruTileStorage(size_t max_tiles)
    : LruTileStorage(max_tiles)
{
  this->root_dir = make_unique_temp_dir("va_lru");
}

DiskLruTileStorage::~DiskLruTileStorage()
{
  std::error_code ec;
  std::filesystem::remove_all(this->root_dir, ec);
}

std::unique_ptr<TileStorage> DiskLruTileStorage::clone() const
{
  return std::make_unique<DiskLruTileStorage>(this->max_tiles);
}

Array &DiskLruTileStorage::get_tile(const TileRegion &region)
{
  std::lock_guard<std::mutex> lock(mutex);

  const TileKey &key = region.key;

  // ---- RAM hit → delegate to base

  auto it = tiles.find(key);
  if (it != tiles.end()) return LruTileStorage::get_tile_no_mutex_lock(region);

  // ---- disk hit

  auto path = tile_path(key);
  if (std::filesystem::exists(path))
  {
    // eviction if needed
    if (this->max_tiles > 0 && tiles.size() >= max_tiles)
    {
      this->evict_unpinned_no_mutex_lock(this->max_tiles - 1);
    }

    Array tile = load_tile_from_disk(region);

    lru.push_front(key);
    auto [inserted_it,
          _] = tiles.emplace(key,
                             LruTileEntry{std::move(tile), lru.begin(), 1});

    return inserted_it->second.value;
  }

  // ---- miss everywhere → empty tile

  return LruTileStorage::get_tile_no_mutex_lock(region);
}

Array DiskLruTileStorage::load_tile_from_disk(const TileRegion &region)
{
  Array tile(region.shape);

  auto path = tile_path(region.key);

  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("DiskLRU: cannot read tile");

  const size_t count = tile.size();
  in.read(reinterpret_cast<char *>(tile.vector.data()), count * sizeof(float));

  return tile;
}

size_t DiskLruTileStorage::max_live_tiles() const
{
  return this->max_tiles;
}

void DiskLruTileStorage::on_evict(const TileKey &key, Array &tile)
{
  auto path = tile_path(key);

  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw std::runtime_error("DiskLRU: cannot write tile");

  const size_t count = tile.size();
  out.write(reinterpret_cast<const char *>(tile.vector.data()),
            count * sizeof(float));
}

std::filesystem::path DiskLruTileStorage::tile_path(const TileKey &key) const
{
  return root_dir / ("tile_" + std::to_string(key.tx) + "_" +
                     std::to_string(key.ty) + ".bin");
}

void DiskLruTileStorage::trim()
{
  std::lock_guard<std::mutex> lock(mutex);

  for (auto &[key, entry] : tiles)
    this->on_evict(key, entry.value); // write-back

  this->tiles.clear();
  this->lru.clear();
}

} // namespace hmap
