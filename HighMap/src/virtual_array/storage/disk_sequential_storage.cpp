/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/string_utils.hpp"
#include "highmap/logger.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/tile_storage.hpp"

namespace hmap
{

DiskSequentialTileStorage::DiskSequentialTileStorage()
{
  root_dir = make_unique_temp_dir("va_seq");
}

DiskSequentialTileStorage::DiskSequentialTileStorage(
    const DiskSequentialTileStorage &other)
{
  root_dir = make_unique_temp_dir("va_seq");
  current_key = other.current_key;
  if (other.current_tile)
  {
    current_tile = *other.current_tile;
  }
  std::error_code ec;
  if (std::filesystem::exists(other.root_dir))
  {
    std::filesystem::copy(other.root_dir,
                          root_dir,
                          std::filesystem::copy_options::recursive |
                              std::filesystem::copy_options::overwrite_existing,
                          ec);
  }
}

DiskSequentialTileStorage::~DiskSequentialTileStorage()
{
  std::error_code ec;
  std::filesystem::remove_all(root_dir, ec);
}

std::unique_ptr<TileStorage> DiskSequentialTileStorage::clone() const
{
  return std::make_unique<DiskSequentialTileStorage>(*this);
}

Array &DiskSequentialTileStorage::get_tile(const TileRegion &region)
{
  if (current_tile)
  {
    hmap::log::error("tile already loaded");
    throw std::logic_error("DiskSequentialTileStorage misuse");
  }

  current_key = region.key;
  current_tile = load_or_create(region);
  return *current_tile;
}

Array DiskSequentialTileStorage::load_or_create(const TileRegion &region)
{
  auto path = tile_path(region.key);

  Array tile(region.shape);

  if (std::filesystem::exists(path))
  {
    std::ifstream in(path, std::ios::binary);
    if (!in)
      throw std::runtime_error("DiskSequential: cannot open tile for reading");

    const std::streamsize expected_bytes = static_cast<std::streamsize>(
        tile.size() * sizeof(float));
    in.read(reinterpret_cast<char *>(tile.vector.data()), expected_bytes);
    if (in.gcount() != expected_bytes)
    {
      throw std::runtime_error(
          "DiskSequential: incomplete tile read from disk");
    }
  }

  return tile;
}

size_t DiskSequentialTileStorage::max_live_tiles() const
{
  return 1;
}

size_t DiskSequentialTileStorage::live_tile_count() const
{
  return current_tile ? 1 : 0;
}

size_t DiskSequentialTileStorage::live_memory_bytes() const
{
  return current_tile
             ? (current_tile->vector.capacity() * sizeof(float) + sizeof(Array))
             : 0;
}

void DiskSequentialTileStorage::release_tile(const TileRegion &region)
{
  if (!current_tile)
  {
    hmap::log::error("no tile loaded on release");
    return;
  }

  if (region.key != current_key) hmap::log::error("releasing wrong tile");

  save_tile(region.key, *current_tile);
  current_tile.reset();
}

void DiskSequentialTileStorage::save_tile(const TileKey &key, const Array &tile)
{
  auto path = tile_path(key);
  auto tmp_path = path;
  tmp_path += ".tmp";

  {
    std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
    if (!out)
      throw std::runtime_error(
          "DiskSequential: cannot write tile to temp file");

    const std::streamsize expected_bytes = static_cast<std::streamsize>(
        tile.size() * sizeof(float));
    out.write(reinterpret_cast<const char *>(tile.vector.data()),
              expected_bytes);
    out.flush();
    if (!out.good())
    {
      throw std::runtime_error("DiskSequential: failed to flush tile to disk");
    }
  }

  std::error_code ec;
  std::filesystem::rename(tmp_path, path, ec);
  if (ec)
  {
    std::filesystem::remove(tmp_path, ec);
    throw std::runtime_error(
        "DiskSequential: cannot rename temp tile to destination");
  }
}

std::filesystem::path DiskSequentialTileStorage::tile_path(
    const TileKey &key) const
{
  return root_dir / ("tile_" + std::to_string(key.tx) + "_" +
                     std::to_string(key.ty) + ".bin");
}

void DiskSequentialTileStorage::trim()
{
  if (current_tile)
  {
    save_tile(current_key, *current_tile);
    current_tile.reset();
  }
}

} // namespace hmap
