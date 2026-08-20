#include "highmap.hpp"

int main(void)
{
  hmap::Path path;
  path.add_point({10.f, 10.f});
  path.add_point({50.f, 100.f});
  path.add_point({200.f, 50.f});
  path.add_point({240.f, 240.f});

  hmap::Path bpath = hmap::catmullrom(path, 10);

  hmap::Array z({256, 256}, 0.f);
  bpath.to_array(z);

  z.to_png("ex_path_catmullrom.png", hmap::Cmap::VIRIDIS);
}
