#include "highmap.hpp"

int main(void)
{
  hmap::Path path;
  path.add_point({10.f, 10.f});
  path.add_point({240.f, 240.f});

  hmap::Array z({256, 256}, 0.f);
  path.to_array(z);
  return 0;
}
