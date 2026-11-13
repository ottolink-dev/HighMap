#include <chrono>
#include <iostream>
#include <string>

#include "highmap.hpp"

int main(int argc, char **argv)
{
  // usage: ex_aeolian_erosion_multiscale [baseline|multiscale]
  const std::string mode = (argc > 1) ? argv[1] : "multiscale";
  const bool        is_baseline = (mode == "baseline");

  hmap::gpu::init_opencl();

  const std::uint32_t seed = 0;
  const int           steps = 500;

  for (int n : {256, 512, 1024})
  {
    glm::ivec2 shape = {n, n};

    // transverse dunes setup, same as ex_aeolian_erosion
    hmap::Array z(shape, 0.05f);

    auto t0 = std::chrono::steady_clock::now();

    if (is_baseline)
    {
      hmap::aeolian_erosion(z,
                            steps,
                            seed,
                            {3.f, 0.f},
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            -1.f,
                            0.001f,
                            shape.x * shape.y,
                            33.f,
                            5.f,
                            15.f,
                            12.f,
                            0.f,
                            0.35f);
    }
    else
    {
      hmap::aeolian_erosion_multiscale(z, steps, seed, {3.f, 0.f});
    }

    auto t1 = std::chrono::steady_clock::now();

    double seconds = std::chrono::duration<double>(t1 - t0).count();

    std::cout << mode + " n=" << n << " time=" << seconds << " s" << std::endl;

    z.dump(mode + "_" + std::to_string(n) + ".png");
    z.to_file(mode + "_" + std::to_string(n) + ".bin");
  }

  return 0;
}
