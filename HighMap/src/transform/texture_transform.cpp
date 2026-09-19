/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <vector>

#include "highmap/internal/validation.hpp"
#include "highmap/texture.hpp"
#include "highmap/transform.hpp"

namespace hmap
{

void flip_lr(Texture &texture)
{
  if (!validate_non_empty(texture)) return;

  for (auto &ch : texture.channels)
  {
    flip_lr(ch);
  }
}

void flip_ud(Texture &texture)
{
  if (!validate_non_empty(texture)) return;

  for (auto &ch : texture.channels)
  {
    flip_ud(ch);
  }
}

void rot90(Texture &texture)
{
  if (!validate_non_empty(texture)) return;

  for (auto &ch : texture.channels)
  {
    rot90(ch);
  }
  texture.shape = glm::ivec2(texture.shape.y, texture.shape.x);
}

void rot180(Texture &texture)
{
  if (!validate_non_empty(texture)) return;

  for (auto &ch : texture.channels)
  {
    rot180(ch);
  }
}

void rot270(Texture &texture)
{
  if (!validate_non_empty(texture)) return;

  for (auto &ch : texture.channels)
  {
    rot270(ch);
  }
  texture.shape = glm::ivec2(texture.shape.y, texture.shape.x);
}

void scale_uv(Texture &texture, glm::vec2 uv_scale)
{
  if (!validate_non_empty(texture)) return;

  for (auto &ch : texture.channels)
  {
    scale_uv(ch, uv_scale);
  }
}

Texture transpose(const Texture &texture)
{
  if (!validate_non_empty(texture)) return Texture();

  std::vector<Array> transposed_channels;
  transposed_channels.reserve(texture.channels.size());
  for (const auto &ch : texture.channels)
  {
    transposed_channels.push_back(transpose(ch));
  }
  return Texture(transposed_channels);
}

} // namespace hmap
