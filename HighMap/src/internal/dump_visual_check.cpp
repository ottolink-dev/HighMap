/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>

#include "highmap/array.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/gradient.hpp"
#include "highmap/math/array.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap
{

void dump_visual_check(const Array &array, const std::string &fname)
{
  Array out = array;
  remap(out);

  // C1 continuity
  {
    Array gn = gradient_norm(array);
    remap(gn);
    out = hstack(out, gn);
  }

  // C2 continuity
  {
    Array dn = abs(laplacian(array));
    remap(dn);
    out = hstack(out, dn);
  }

  // tiling
  {
    // line 3 x 1
    Array a = array;
    a = hstack(a, array);
    a = hstack(a, array);

    // stack vertically
    Array b = a;
    b = vstack(b, a);
    b = vstack(b, a);

    // reduce 3 x 3 to initial shape
    Array bs = b.resample_to_shape(array.shape);
    out = hstack(out, bs);

    Array gn = gradient_norm(bs);
    remap(gn);
    out = hstack(out, gn);
  }

  out.to_png(fname, hmap::Cmap::JET);
}

} // namespace hmap
