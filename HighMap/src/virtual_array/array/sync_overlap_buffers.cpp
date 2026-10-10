/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <memory>

#include "highmap/array.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/core.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/tile_storage.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

void VirtualArray::sync_overlap_buffers(SyncOperation op)
{
  if (this->storage->max_live_tiles() < 2)
  {
    hmap::log::error(
        "sync_overlap_buffers requires at least 2 tiles in memory, skipping");
    return;
  }

  int nx = ceil_div(this->shape.x, this->tile_shape.x);
  int ny = ceil_div(this->shape.y, this->tile_shape.y);

  // --- x-direction

  for (int ty = 0; ty < ny; ++ty)
    for (int tx = 0; tx < nx - 1; ++tx)
    {
      // load
      TileRegion region0 = this->tile_region_from_tile_coords(tx, ty);
      TileRegion region1 = this->tile_region_from_tile_coords(tx + 1, ty);
      Array     &tile0 = this->storage->get_tile(region0);
      Array     &tile1 = this->storage->get_tile(region1);

      switch (op)
      {
      case SyncOperation::Average:
      case SyncOperation::Mean:
        for (int p = 0; p < this->halo; p++)
        {
          int pbuf = tile0.shape.x - 2 * this->halo + p;
          for (int q = 0; q < tile0.shape.y; q++)
          {
            float avg = 0.5f * (tile0(pbuf, q) + tile1(p, q));
            tile0(pbuf, q) = avg;
            tile1(p, q) = avg;
          }
        }
        break;

      case SyncOperation::CopyFirst:
        for (int p = 0; p < this->halo; p++)
        {
          int pbuf = tile0.shape.x - 2 * this->halo + p;
          for (int q = 0; q < tile0.shape.y; q++)
            tile1(p, q) = tile0(pbuf, q);
        }
        break;

      case SyncOperation::CopySecond:
        for (int p = 0; p < this->halo; p++)
        {
          int pbuf = tile0.shape.x - 2 * this->halo + p;
          for (int q = 0; q < tile0.shape.y; q++)
            tile0(pbuf, q) = tile1(p, q);
        }
        break;

      case SyncOperation::Min:
        for (int p = 0; p < this->halo; p++)
        {
          int pbuf = tile0.shape.x - 2 * this->halo + p;
          for (int q = 0; q < tile0.shape.y; q++)
          {
            float m = std::min(tile0(pbuf, q), tile1(p, q));
            tile0(pbuf, q) = m;
            tile1(p, q) = m;
          }
        }
        break;

      case SyncOperation::Max:
        for (int p = 0; p < this->halo; p++)
        {
          int pbuf = tile0.shape.x - 2 * this->halo + p;
          for (int q = 0; q < tile0.shape.y; q++)
          {
            float m = std::max(tile0(pbuf, q), tile1(p, q));
            tile0(pbuf, q) = m;
            tile1(p, q) = m;
          }
        }
        break;

      case SyncOperation::SmoothBlend:
        for (int p = 0; p < this->halo; p++)
        {
          float r = float(p) / float(this->halo - 1);
          r = smoothstep5(r);
          int pbuf = tile0.shape.x - 2 * this->halo + p;

          for (int q = 0; q < tile0.shape.y; q++)
          {
            tile1(p, q) = lerp(tile0(pbuf, q), tile1(p, q), r);
            tile0(pbuf, q) = tile1(p, q);
          }
        }
        break;
      }

      // release
      this->storage->release_tile(region0);
      this->storage->release_tile(region1);
    }

  // --- y-direction

  for (int ty = 0; ty < ny - 1; ++ty)
    for (int tx = 0; tx < nx; ++tx)
    {
      // load
      TileRegion region0 = this->tile_region_from_tile_coords(tx, ty);
      TileRegion region1 = this->tile_region_from_tile_coords(tx, ty + 1);
      Array     &tile0 = this->storage->get_tile(region0);
      Array     &tile1 = this->storage->get_tile(region1);

      switch (op)
      {
      case SyncOperation::Average:
      case SyncOperation::Mean:
        for (int q = 0; q < this->halo; q++)
        {
          int qbuf = tile0.shape.y - 2 * this->halo + q;
          for (int p = 0; p < tile0.shape.x; p++)
          {
            float avg = 0.5f * (tile0(p, qbuf) + tile1(p, q));
            tile0(p, qbuf) = avg;
            tile1(p, q) = avg;
          }
        }
        break;

      case SyncOperation::CopyFirst:
        for (int q = 0; q < this->halo; q++)
        {
          int qbuf = tile0.shape.y - 2 * this->halo + q;
          for (int p = 0; p < tile0.shape.x; p++)
            tile1(p, q) = tile0(p, qbuf);
        }
        break;

      case SyncOperation::CopySecond:
        for (int q = 0; q < this->halo; q++)
        {
          int qbuf = tile0.shape.y - 2 * this->halo + q;
          for (int p = 0; p < tile0.shape.x; p++)
            tile0(p, qbuf) = tile1(p, q);
        }
        break;

      case SyncOperation::Min:
        for (int q = 0; q < this->halo; q++)
        {
          int qbuf = tile0.shape.y - 2 * this->halo + q;
          for (int p = 0; p < tile0.shape.x; p++)
          {
            float m = std::min(tile0(p, qbuf), tile1(p, q));
            tile0(p, qbuf) = m;
            tile1(p, q) = m;
          }
        }
        break;

      case SyncOperation::Max:
        for (int q = 0; q < this->halo; q++)
        {
          int qbuf = tile0.shape.y - 2 * this->halo + q;
          for (int p = 0; p < tile0.shape.x; p++)
          {
            float m = std::max(tile0(p, qbuf), tile1(p, q));
            tile0(p, qbuf) = m;
            tile1(p, q) = m;
          }
        }
        break;

      case SyncOperation::SmoothBlend:
        for (int q = 0; q < this->halo; q++)
        {
          float r = float(q) / float(this->halo - 1);
          r = smoothstep5(r);
          int qbuf = tile0.shape.y - 2 * this->halo + q;

          for (int p = 0; p < tile0.shape.x; p++)
          {
            tile1(p, q) = lerp(tile0(p, qbuf), tile1(p, q), r);
            tile0(p, qbuf) = tile1(p, q);
          }
        }
        break;
      }

      // release
      this->storage->release_tile(region0);
      this->storage->release_tile(region1);
    }
}

} // namespace hmap
