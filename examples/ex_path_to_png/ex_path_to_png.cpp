#include "highmap.hpp"

int main(void)
{
  hmap::Path path;
  path.push_back({10.f, 10.f});
  path.push_back({240.f, 240.f});

  path.to_png("ex_path_to_png_output.png");
  return 0;
}
