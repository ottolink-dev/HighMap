#include "highmap.hpp"

int main(void)
{
  int seed = 1;

  glm::vec4 bbox = {-1.f, 2.f, 0.f, 5.f};

  // create a cloud of points and convert it to a graph using Delaunay
  // triangulation
  int         npoints = 10;
  hmap::Cloud cloud = hmap::Cloud(npoints, seed, bbox);
  hmap::Graph graph = cloud.to_graph_delaunay();

  graph.print();
  graph.to_png("ex_graph0.png");
  graph.to_csv("ex_graph_nodes.csv", "ex_graph_edges.csv");
}
