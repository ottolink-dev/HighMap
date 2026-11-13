#include <chrono>
#include <iostream>
#include <string>

#include "highmap.hpp"

int main(int argc, char **argv)
{
  // usage: ex_aeolian_erosion_multiscale [baseline|multiscale|shadowmap1|shadowmap4]
  const std::string mode = (argc > 1) ? argv[1] : "multiscale";
  const bool        is_baseline = (mode == "baseline") ||
                            (mode == "baseline_s1") ||
                            (mode.rfind("shadowmap", 0) == 0);
  int               shadow_map_interval = 0;

  if (mode == "shadowmap1") shadow_map_interval = 1;
  if (mode == "shadowmap4") shadow_map_interval = 4;

  hmap::gpu::init_opencl();

  // "baseline_s1" is the baseline with another seed (noise floor reference)
  const std::uint32_t seed = (mode == "baseline_s1") ? 1u : 0u;
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
                            0.35f,
                            0,
                            hmap::TileRegion::dummy(),
                            hmap::SingleArrayComm::instance(),
                            shadow_map_interval);
    }
    else
    {
      hmap::aeolian_erosion_multiscale(z, steps, seed, {3.f, 0.f});
    }

    auto t1 = std::chrono::steady_clock::now();

    double seconds = std::chrono::duration<double>(t1 - t0).count();

    std::cout << mode + " n=" << n << " time=" << seconds << " s" << std::endl;

    // RMS difference against the saved baseline, if any
    {
      hmap::Array z_ref(z.shape);
      z_ref.from_file("baseline_" + std::to_string(n) + ".bin");

      if (z_ref.shape == z.shape && mode != "baseline")
      {
        hmap::Array diff = z - z_ref;
        diff = diff * diff;

        std::cout << "  rmse vs baseline = " << std::sqrt(diff.mean())
                  << " (baseline ptp = " << z_ref.ptp() << ")" << std::endl;
      }
    }

    z.dump(mode + "_" + std::to_string(n) + ".png");
    z.to_file(mode + "_" + std::to_string(n) + ".bin");
  }

  return 0;
}
