#include "highmap.hpp"

int main(void)
{
  hmap::Path path;
  path.push_back({10.f, 10.f});
  path.push_back({100.f, 10.f});
  path.push_back({100.f, 100.f});
  path.push_back({10.f, 100.f});
  path.push_back({10.f, 10.f}); // loop

  hmap::Path clean_path = hmap::remove_geometric_loops(path);
  return 0;
}
