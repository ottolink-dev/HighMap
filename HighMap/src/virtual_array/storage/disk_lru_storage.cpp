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
  auto new_storage = std::make_unique<DiskLruTileStorage>(this->max_tiles);
  std::error_code ec;
  if (std::filesystem::exists(this->root_dir))
  {
    std::filesystem::copy(this->root_dir,
                          new_storage->root_dir,
                          std::filesystem::copy_options::recursive |
                              std::filesystem::copy_options::overwrite_existing,
                          ec);
  }
  return new_storage;
}

Array &DiskLruTileStorage::get_tile(const TileRegion &region)
{
  std::lock_guard<std::mutex> lock(mutex);

  const TileKey &key = region.key;

  // ---- RAM hit → delegate to base

  auto it = tiles.find(key);
  if (it != tiles.end())
  {
    return LruTileStorage::get_tile_no_mutex_lock(region);
  }

  // ---- disk hit

  auto path = tile_path(key);
  if (std::filesystem::exists(path))
  {
    // eviction if needed
    if (this->max_tiles > 0 && this->lru.size() >= this->max_tiles)
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
  if (!in) throw std::runtime_error("DiskLRU: cannot open tile for reading");

  const std::streamsize expected_bytes = static_cast<std::streamsize>(
      tile.size() * sizeof(float));
  in.read(reinterpret_cast<char *>(tile.vector.data()), expected_bytes);
  if (in.gcount() != expected_bytes)
  {
    throw std::runtime_error("DiskLRU: incomplete tile read from disk");
  }

  return tile;
}

size_t DiskLruTileStorage::max_live_tiles() const
{
  return this->max_tiles;
}

void DiskLruTileStorage::on_evict(const TileKey &key, Array &tile)
{
  auto path = tile_path(key);
  auto tmp_path = path;
  tmp_path += ".tmp";

  {
    std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
    if (!out)
      throw std::runtime_error("DiskLRU: cannot write tile to temp file");

    const std::streamsize expected_bytes = static_cast<std::streamsize>(
        tile.size() * sizeof(float));
    out.write(reinterpret_cast<const char *>(tile.vector.data()),
              expected_bytes);
    out.flush();
    if (!out.good())
    {
      throw std::runtime_error("DiskLRU: failed to flush tile to disk");
    }
  }

  std::error_code ec;
  std::filesystem::rename(tmp_path, path, ec);
  if (ec)
  {
    std::filesystem::remove(tmp_path, ec);
    throw std::runtime_error("DiskLRU: cannot rename temp tile to destination");
  }
}

std::filesystem::path DiskLruTileStorage::tile_path(const TileKey &key) const
{
  return root_dir / ("tile_" + std::to_string(key.tx) + "_" +
                     std::to_string(key.ty) + ".bin");
}

void DiskLruTileStorage::trim()
{
  std::lock_guard<std::mutex> lock(mutex);

  this->evict_unpinned_no_mutex_lock(0);
}

} // namespace hmap
