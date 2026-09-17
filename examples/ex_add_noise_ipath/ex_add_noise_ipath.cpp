#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};

  float kw = 4.f;
  int   seed = 0;
  float amp = 128.f;

  hmap::CellPath ipath;

  int i0 = 32;
  int i1 = shape.x - i0;
  int j0 = 32;
  int j1 = 3 * j0;

  for (int i = i0; i < i1; ++i)
  {
    int j = int(float(i - i0) / (i1 - i0) * (j1 - j0) + j0);
    ipath.push_back({i, j});
  }

  hmap::Array a0(shape);
  for (const auto &p : ipath)
    a0(p) = 1.f;

  hmap::add_noise(ipath, seed, kw, amp, shape);

  hmap::Array a1(shape);
  for (const auto &p : ipath)
    a1(p) = 1.f;

  hmap::export_banner_png("ex_add_noise_ipath.png", {a0, a1}, hmap::Cmap::GRAY);
}
