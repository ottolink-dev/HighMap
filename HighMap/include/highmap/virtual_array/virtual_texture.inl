/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

template <typename Func>
void for_each_tile(VirtualTexture &tex, Func &&func, const ComputeMode &cm)
{
  std::vector<VirtualArray *> p_vas;
  p_vas.reserve(tex.channels());

  for (auto &ch : tex.get_arrays())
    p_vas.push_back(&ch);

  for_each_tile(p_vas, std::forward<Func>(func), cm);
}

template <typename Func>
void for_each_tile(const VirtualTexture &tex,
                   Func                &&func,
                   const ComputeMode    &cm)
{
  auto &mutable_tex = const_cast<VirtualTexture &>(tex);
  for_each_tile(mutable_tex, std::forward<Func>(func), cm);
}

template <typename Func>
void for_each_pixel(VirtualTexture &tex, Func &&func, const ComputeMode &cm)
{
  for_each_tile(
      tex,
      [&](std::vector<Array *> &chs, const TileRegion &region)
      {
        int w = region.shape.x;
        int h = region.shape.y;
        int c = int(chs.size());

        std::vector<float> px(c);

        for (int y = 0; y < h; ++y)
          for (int x = 0; x < w; ++x)
          {
            for (int k = 0; k < c; ++k)
              px[k] = (*chs[k])(x, y);

            func(px, x, y, region);

            for (int k = 0; k < c; ++k)
              (*chs[k])(x, y) = px[k];
          }
      },
      cm);
}
