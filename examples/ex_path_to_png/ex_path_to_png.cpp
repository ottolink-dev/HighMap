#include "highmap.hpp"

int main(void)
{
  hmap::Path path;
  path.add_point({10.f, 10.f});
  path.add_point({240.f, 240.f});

  path.to_png("ex_path_to_png_output.png");
  return 0;
}
