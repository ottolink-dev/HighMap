/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file tile_storage.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2025
 */
#pragma once
#include <filesystem>
#include <fstream>
#include <list>
#include <mutex>
#include <optional>

#include "highmap/array.hpp"
#include "highmap/virtual_array/tile_region.hpp"

#include <unordered_map>

namespace hmap
{

// =====================================
// Storage Modes
// =====================================

enum StorageMode : int
{
  VA_RAM = 0,
  VA_DISK_LRU_STREAM = 1,
  VA_DISK_LRU_FULL = 2
};

std::string to_string(StorageMode m);

struct TileKeyHash
{
  size_t operator()(const TileKey &k) const
  {
    if constexpr (sizeof(size_t) >= 8)
      return (static_cast<size_t>(k.tx) << 32) ^ static_cast<size_t>(k.ty);
    else
      return (static_cast<size_t>(k.tx) << 16) ^ static_cast<size_t>(k.ty);
  }
};

// =====================================
// Abstract TileStorage
// =====================================

class TileStorage
{
public:
  virtual ~TileStorage() = default;

  virtual std::unique_ptr<TileStorage> clone() const = 0;

  virtual Array      &get_tile(const TileRegion &region) = 0;
  virtual std::string info_string() const = 0;
  virtual size_t      live_tile_count() const = 0;
  virtual size_t      live_memory_bytes() const = 0;
  virtual size_t      max_live_tiles() const = 0;
  virtual void        release_tile(const TileRegion &region) = 0;
  virtual void        trim();
};

// =====================================
// RAM Storage
// =====================================

class RamTileStorage : public TileStorage
{
public:
  RamTileStorage() = default;
  RamTileStorage(const RamTileStorage &other);

  std::unique_ptr<TileStorage> clone() const override;

  Array      &get_tile(const TileRegion &region) override;
  std::string info_string() const override;
  size_t      live_tile_count() const override;
  size_t      live_memory_bytes() const override;
  size_t      max_live_tiles() const override;
  void        release_tile(const TileRegion &region) override;

private:
  mutable std::mutex                              mtx;
  std::unordered_map<TileKey, Array, TileKeyHash> tiles;
};

// =====================================
// Disk LRU Storage
// =====================================

struct LruTileEntry
{
  Array                        value;
  std::list<TileKey>::iterator lru_it;
  size_t                       pin_count = 0;
  bool                         dirty = true;
};

class DiskLruTileStorage : public TileStorage
{
public:
  explicit DiskLruTileStorage(
      StorageMode mode = StorageMode::VA_DISK_LRU_STREAM);
  ~DiskLruTileStorage() override;

  std::unique_ptr<TileStorage> clone() const override;

  Array      &get_tile(const TileRegion &region) override;
  std::string info_string() const override;
  size_t      live_tile_count() const override;
  size_t      live_memory_bytes() const override;
  size_t      max_live_tiles() const override;
  void        release_tile(const TileRegion &region) override;
  void        trim() override;

  StorageMode get_mode() const;

private:
  StorageMode                                            mode;
  std::filesystem::path                                  root_dir;
  std::list<TileKey>                                     lru;
  std::unordered_map<TileKey, LruTileEntry, TileKeyHash> tiles;
  mutable std::mutex                                     mutex;

  std::filesystem::path tile_path(const TileKey &key) const;
  Array                 load_tile_from_disk(const TileRegion &region);
  void save_tile_to_disk(const TileKey &key, const Array &tile);
  void evict_entry_no_mutex_lock(
      std::unordered_map<TileKey, LruTileEntry, TileKeyHash>::iterator it);
  void evict_all_unpinned_no_mutex_lock();
};

// =====================================
// Factory function
// =====================================

std::unique_ptr<TileStorage> make_storage(glm::ivec2  shape,
                                          glm::ivec2  tile_shape,
                                          StorageMode storage_mode);

} // namespace hmap
