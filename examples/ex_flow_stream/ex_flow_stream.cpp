#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  // shape = {1024, 1024};
  glm::vec2 kw = {4.f, 4.f};
  int       seed = 1;

  // base terrain with Perlin noise
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  // z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_SMOOTH_COSINE, 1.f);

  z += hmap::circus(shape);

  hmap::remap(z);

  // generate triangular terrain mesh from heightmap
  hmap::TerrainTriMesh mesh =
      hmap::generate_terrain_tri_mesh_from_heightmap_random(z, 2 * 4096, seed);

  // multiple starting points across the slopes of the valley head
  std::vector<glm::vec2> start_points = {
      {0.39f, 0.78f},
      {0.39f, 0.31f},
      {0.31f, 0.39f},
      {0.35f, 0.43f},
      {0.27f, 0.49f},
      {0.27f, 0.59f},
      {0.24f, 0.67f},
      {0.69f, 0.63f},
      {0.61f, 0.82f},
      {0.40f, 0.45f},
      {0.50f, 0.70f},
  };

  std::vector<hmap::Path> streams;

  for (const auto &start : start_points)
  {
    hmap::Path stream = hmap::flow_stream(mesh,
                                          start,
                                          0.95f,
                                          2.f,
                                          100.f,
                                          true,
                                          streams,
                                          0.02f,
                                          0.4f,
                                          true);

    if (!stream.empty()) streams.push_back(stream);
  }

  // convert merged streams into a connected Graph
  hmap::Graph stream_graph = hmap::hydrology_paths_to_graph(streams, 1e-3f);

  // enforce strict monotonic downflow along all stream branches
  stream_graph.enforce_downhill_elevation(0.01f);

  // visualize graph on an empty grid
  glm::vec4   bbox(0.f, 1.f, 0.f, 1.f);
  hmap::Array z_graph(shape, 0.f);
  stream_graph.to_array(z_graph, bbox, false);

  // carve stream pixel-wise for visual checking
  hmap::Array z_carved = z;
  float       depth = 0.05f;

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      if (z_graph(i, j) > 0.f)
        z_carved(i, j) = std::max(0.f, z(i, j) - depth);

  hmap::export_banner_png("ex_flow_stream.png",
                          {z, z_graph, z_carved},
                          hmap::Cmap::JET,
                          true);
}
