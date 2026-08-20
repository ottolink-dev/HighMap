/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <iostream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "highmap/openmp.hpp"

namespace hmap
{

bool init_openmp(int num_threads)
{
#ifdef _OPENMP
  if (num_threads <= 0) num_threads = omp_get_num_procs();
  omp_set_num_threads(num_threads);
  log_openmp_info();
  return true;
#else
  std::cout << "[HighMap/OpenMP] Not enabled (no _OPENMP macro)\n";
  return false;
#endif
}

void log_openmp_info()
{
#ifdef _OPENMP
  std::cout << "[HighMap/OpenMP] Enabled\n";

  // Version
  std::cout << "[HighMap/OpenMP] Version: " << _OPENMP << "\n";

  // Number of processors
  std::cout << "[HighMap/OpenMP] Num processors: " << omp_get_num_procs()
            << "\n";

  // Max threads
  std::cout << "[HighMap/OpenMP] Max threads: " << omp_get_max_threads()
            << "\n";

  // Dynamic threads
  std::cout << "[HighMap/OpenMP] Dynamic threads: "
            << (omp_get_dynamic() ? "ON" : "OFF") << "\n";

  // Nested parallelism
  std::cout << "[HighMap/OpenMP] Nested parallelism: "
            << (omp_get_nested() ? "ON" : "OFF") << "\n";

// Parallel region test
#pragma omp parallel
  {
#pragma omp single
    {
      std::cout << "[HighMap/OpenMP] Actual threads in parallel region: "
                << omp_get_num_threads() << "\n";
    }
  }

#else
  std::cout << "[HighMap/OpenMP] Not enabled (no _OPENMP macro)\n";
#endif
}

} // namespace hmap
