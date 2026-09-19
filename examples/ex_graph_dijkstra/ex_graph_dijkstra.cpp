#include <iostream>

#include "highmap.hpp"

int main(void)
{
  int       seed = 1;
  glm::vec4 bbox = {-1.f, 2.f, 0.f, 5.f};

  // create a cloud of points and convert it to a graph using Delaunay
  // triangulation
  int         npoints = 15;
  hmap::Cloud cloud = hmap::Cloud(npoints, seed, bbox);
  hmap::Graph graph = cloud.to_graph_delaunay();

  graph.print();

  int              i_point_start = 0;
  int              i_point_end = 4;
  std::vector<int> route = graph.dijkstra(i_point_start, i_point_end);

  std::cout << "Route (point indices):\n";
  for (auto &p : route)
    std::cout << p << "\n";
}
