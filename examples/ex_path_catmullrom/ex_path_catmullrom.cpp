#include "highmap.hpp"

int main(void)
{
  hmap::Path path;
  path.push_back({10.f, 10.f});
  path.push_back({50.f, 100.f});
  path.push_back({200.f, 50.f});
  path.push_back({240.f, 240.f});

  hmap::Path bpath = hmap::catmullrom(path, 10);

  hmap::Array z({256, 256}, 0.f);
  bpath.to_array(z);

  z.to_png("ex_path_catmullrom.png", hmap::Cmap::VIRIDIS);
}
