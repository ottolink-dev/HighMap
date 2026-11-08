/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
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

DiskLruTileStorage::DiskLruTileStorage(StorageMode mode) : mode(mode)
{
  this->root_dir = make_unique_temp_dir("va_disk");
}

DiskLruTileStorage::~DiskLruTileStorage()
{
  std::error_code ec;
  std::filesystem::remove_all(this->root_dir, ec);
}

std::unique_ptr<TileStorage> DiskLruTileStorage::clone() const
{
  std::lock_guard<std::mutex> lock(this->mutex);
  // Flush any in-memory tiles to disk first so copied directory is up to date
  const_cast<DiskLruTileStorage *>(this)->evict_all_unpinned_no_mutex_lock();

  auto new_storage = std::make_unique<DiskLruTileStorage>(this->mode);
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

void DiskLruTileStorage::evict_entry_no_mutex_lock(
    std::unordered_map<TileKey, LruTileEntry, TileKeyHash>::iterator it)
{
  if (it == this->tiles.end() || it->second.pin_count > 0) return;

  if (it->second.dirty)
  {
    this->save_tile_to_disk(it->first, it->second.value);
  }

  this->lru.erase(it->second.lru_it);
  this->tiles.erase(it);
}

void DiskLruTileStorage::evict_all_unpinned_no_mutex_lock()
{
  for (auto it = this->tiles.begin(); it != this->tiles.end();)
  {
    if (it->second.pin_count == 0)
    {
      if (it->second.dirty)
      {
        this->save_tile_to_disk(it->first, it->second.value);
      }
      this->lru.erase(it->second.lru_it);
      it = this->tiles.erase(it);
    }
    else
    {
      ++it;
    }
  }
}

StorageMode DiskLruTileStorage::get_mode() const
{
  return this->mode;
}

Array &DiskLruTileStorage::get_tile(const TileRegion &region)
{
  std::lock_guard<std::mutex> lock(this->mutex);
  const TileKey              &key = region.key;

  // 1. RAM hit
  auto it = this->tiles.find(key);
  if (it != this->tiles.end())
  {
    this->lru.splice(this->lru.begin(), this->lru, it->second.lru_it);
    it->second.pin_count++;
    it->second.dirty = true;
    return it->second.value;
  }

  // 2. Disk hit
  auto path = this->tile_path(key);
  if (std::filesystem::exists(path))
  {
    Array tile = this->load_tile_from_disk(region);
    this->lru.push_front(key);
    auto [inserted_it, _] = this->tiles.emplace(
        key,
        LruTileEntry{std::move(tile), this->lru.begin(), 1, true});
    return inserted_it->second.value;
  }

  // 3. New tile allocation
  Array tile(region.shape);
  this->lru.push_front(key);
  auto [inserted_it, _] = this->tiles.emplace(
      key,
      LruTileEntry{std::move(tile), this->lru.begin(), 1, true});
  return inserted_it->second.value;
}

void DiskLruTileStorage::release_tile(const TileRegion &region)
{
  std::lock_guard<std::mutex> lock(this->mutex);
  auto                        it = this->tiles.find(region.key);
  if (it == this->tiles.end()) return;

  if (it->second.pin_count > 0)
  {
    it->second.pin_count--;
  }

  // Strict stream: as soon as pin count drops to 0, immediately evict to disk
  if (this->mode == StorageMode::VA_DISK_LRU_STREAM &&
      it->second.pin_count == 0)
  {
    this->evict_entry_no_mutex_lock(it);
  }
}

size_t DiskLruTileStorage::max_live_tiles() const
{
  return (this->mode == StorageMode::VA_DISK_LRU_STREAM)
             ? 0
             : std::numeric_limits<size_t>::max();
}

size_t DiskLruTileStorage::live_tile_count() const
{
  std::lock_guard<std::mutex> lock(this->mutex);
  return this->tiles.size();
}

size_t DiskLruTileStorage::live_memory_bytes() const
{
  std::lock_guard<std::mutex> lock(this->mutex);
  size_t                      bytes = 0;
  for (const auto &[key, entry] : this->tiles)
  {
    bytes += entry.value.vector.capacity() * sizeof(float) +
             sizeof(LruTileEntry) + sizeof(key);
  }
  return bytes;
}

void DiskLruTileStorage::trim()
{
  std::lock_guard<std::mutex> lock(this->mutex);
  this->evict_all_unpinned_no_mutex_lock();
}

std::string DiskLruTileStorage::info_string() const
{
  return (this->mode == StorageMode::VA_DISK_LRU_STREAM) ? "DiskLRU (Stream)"
                                                         : "DiskLRU (Full)";
}

Array DiskLruTileStorage::load_tile_from_disk(const TileRegion &region)
{
  Array tile(region.shape);
  auto  path = this->tile_path(region.key);

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

void DiskLruTileStorage::save_tile_to_disk(const TileKey &key,
                                           const Array   &tile)
{
  auto path = this->tile_path(key);
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
  return this->root_dir / ("tile_" + std::to_string(key.tx) + "_" +
                           std::to_string(key.ty) + ".bin");
}

} // namespace hmap
