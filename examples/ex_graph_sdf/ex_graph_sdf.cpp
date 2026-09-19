#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 1;

  glm::vec4 bbox = {0.2f, 0.8f, 0.2f, 0.8f};
  glm::vec4 bbox_array = {0.f, 1.f, 0.f, 1.f};

  // create a cloud of points and convert it to a graph using Delaunay
  // triangulation
  int         npoints = 10;
  hmap::Cloud cloud = hmap::Cloud(npoints, seed, bbox);
  hmap::Graph graph_delaunay = cloud.to_graph_delaunay();

  hmap::Array zd = graph_delaunay.to_array_sdf(shape, bbox_array);

  // convert to a minimum spanning tree
  hmap::Graph graph_mst = graph_delaunay.minimum_spanning_tree_prim();

  hmap::Array zm = graph_mst.to_array_sdf(shape, bbox_array);

  hmap::export_banner_png("ex_graph_sdf.png", {zd, zm}, hmap::Cmap::JET);
}
