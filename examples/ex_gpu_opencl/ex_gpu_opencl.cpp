#include <iostream>

#include "highmap.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/internal/timer.hpp"

template <typename F1, typename F2>
void compare(F1 fct1, F2 fct2, float tolerance, const std::string &fname)
{
  glm::ivec2 shape = {256, 512};
  shape = {512, 512};
  // shape = {1024, 1024};
  // shape = {2048, 2048};
  // shape = {4096 * 2, 4096 * 2};
  glm::vec2 kw = {2.f, 2.f};
  int       seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);
  // hmap::zeroed_edges(z);

  z.to_png_grayscale("out0.png", CV_16U);

  hmap::Array z1 = z;
  hmap::Array z2 = z;

  // host
  hmap::Timer::Start("host");
  fct1(z1);
  hmap::Timer::Stop("host");

  z1.to_png_grayscale("out1.png", CV_16U);

  // GPU
  hmap::Timer::Start("GPU");
  fct2(z2);
  hmap::Timer::Stop("GPU");

  z2.to_png_grayscale("out2.png", CV_16U);

  hmap::AssertResults res;
  hmap::assert_almost_equal(z1, z2, tolerance, fname, &res);
  res.msg += "[" + fname + "]";
  res.print();
}

// ---

int main(void)
{

  hmap::gpu::init_opencl();

  // int ir = 32;

  // clwrapper::KernelManager::get_instance().set_block_size(32);

  // compare([](hmap::Array &z) { hmap::median_3x3(z); },
  //         [](hmap::Array &z) { hmap::gpu::median_3x3(z); },
  //         1e-3f,
  //         "diff_median_3x3.png");

  // compare([ir](hmap::Array &z)
  //         { hmap::hydraulic_stream_log(z, 0.1f, 5.f / 512.f, 64); },
  //         [ir](hmap::Array &z)
  //         { hmap::gpu::hydraulic_stream_log(z, 0.1f, 5.f / 512.f, 64); },
  //         1e-3f,
  //         "hydraulic_stream_log.png");

  // {
  //   glm::vec4 bbox = {1.f, 2.f, -0.5f, 0.5f};
  //   hmap::Path path = hmap::Path(200, 0, bbox.adjust(0.2f, -0.2f, 0.2f,
  //   -0.2f)); path.reorder_nns();

  //   compare([bbox, path](hmap::Array &z)
  //           { z = hmap::sdf_2d_polyline(path, z.shape, bbox); },
  //           [bbox, path](hmap::Array &z)
  //           { z = hmap::gpu::sdf_2d_polyline(path, z.shape, bbox); },
  //           1e-3f,
  //           "sdf_2d_polyline.png");
  // }

  // {
  //   std::vector<hmap::NoiseType> types = {
  //       // hmap::NoiseType::PARBERRY,
  //       hmap::NoiseType::PERLIN,
  //       hmap::NoiseType::PERLIN_BILLOW,
  //       hmap::NoiseType::PERLIN_HALF,
  //       hmap::NoiseType::SIMPLEX2,
  //       // hmap::NoiseType::SIMPLEX2S,
  //       hmap::NoiseType::VALUE,
  //       hmap::NoiseType::VALUE_CUBIC,
  //       // hmap::NoiseType::VALUE_DELAUNAY,
  //       hmap::NoiseType::VALUE_LINEAR,
  //       hmap::NoiseType::WORLEY
  //       // hmap::NoiseType::WORLEY_DOUBLE,
  //       // hmap::NoiseType::WORLEY_VALUE
  //   };

  //   glm::vec2 kw = {32.f, 32.f};

  //   for (auto type : types)
  //   {
  //     compare([type, shape, kw](hmap::Array &z)
  //             { z = hmap::noise(type, shape, kw, 1); },
  //             [type, shape, kw](hmap::Array &z)
  //             { z = hmap::gpu::noise(type, shape, kw, 1); },
  //             1e-3f,
  //             "noise" + std::to_string(type) + ".png");
  //   }

  //   // fbm
  //   kw = {4.f, 4.f};
  //   for (auto type : types)
  //   {
  //     compare(
  //         [type, shape, kw](hmap::Array &z)
  //         {
  //           if (type != hmap::NoiseType::VALUE_LINEAR)
  //             z = hmap::noise_fbm(type, shape, kw, 1);
  //           else
  //             z = 0.f;
  //         },
  //         [type, shape, kw](hmap::Array &z)
  //         { z = hmap::gpu::noise_fbm(type, shape, kw, 1); },
  //         1e-3f,
  //         "noise_fbm" + std::to_string(type) + ".png");
  //   }
  // }
}
