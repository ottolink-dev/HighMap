#include "highmap.hpp"

int main(void)
{
  hmap::Path path;
  path.add_point({10.f, 10.f});
  path.add_point({100.f, 10.f});
  path.add_point({100.f, 100.f});
  path.add_point({10.f, 100.f});
  path.add_point({10.f, 10.f}); // loop

  hmap::Path clean_path = hmap::remove_geometric_loops(path);
  return 0;
}
