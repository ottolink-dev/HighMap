/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>

#include "highmap/array.hpp"
#include "highmap/export/image.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/texture.hpp"

namespace hmap
{

Texture compute_splatmap(const Array *p_r,
                         const Array *p_g,
                         const Array *p_b,
                         const Array *p_a)
{
  if (!p_r || !validate_non_empty(*p_r)) return Texture();
  if (p_g && !validate_same_shape(*p_r, *p_g)) return Texture();
  if (p_b && !validate_same_shape(*p_r, *p_b)) return Texture();
  if (p_a && !validate_same_shape(*p_r, *p_a)) return Texture();

  Texture smap = Texture(p_r->shape, 4);

  smap[0] = *p_r;

  if (p_g) smap[1] = *p_g;
  if (p_b) smap[2] = *p_b;
  if (p_a) smap[3] = *p_a;
  return smap;
}

void export_splatmap_png(const std::string &fname,
                         const Array       *p_r,
                         const Array       *p_g,
                         const Array       *p_b,
                         const Array       *p_a,
                         int                depth)
{
  if (!p_r || !validate_non_empty(*p_r)) return;

  Texture smap = compute_splatmap(p_r, p_g, p_b, p_a);
  smap.to_png(fname, depth);
}

} // namespace hmap
