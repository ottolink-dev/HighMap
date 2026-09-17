#include "highmap.hpp"

int main(void)
{
  hmap::Path path;
  path.push_back({10.f, 10.f});
  path.push_back({240.f, 240.f});

  hmap::Array z({256, 256}, 0.f);
  path.to_array(z);
  return 0;
}
