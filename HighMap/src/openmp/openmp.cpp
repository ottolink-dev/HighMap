/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "highmap/logger.hpp"
#include "highmap/openmp.hpp"

namespace hmap
{

bool init_openmp(int num_threads)
{
#ifdef _OPENMP
  if (num_threads > 0)
    omp_set_num_threads(num_threads);
  else
    omp_set_num_threads(omp_get_num_procs());

  log_openmp_info();
  return true;
#else
  hmap::log::info("OpenMP not enabled (no _OPENMP macro)");
  return false;
#endif
}

void log_openmp_info()
{
#ifdef _OPENMP
  hmap::log::info("OpenMP enabled");

  // Version
  hmap::log::info("Version: {}", _OPENMP);

  // Number of processors
  hmap::log::info("Num processors: {}", omp_get_num_procs());

  // Max threads
  hmap::log::info("Max threads: {}", omp_get_max_threads());

  // Dynamic threads
  hmap::log::info("Dynamic threads: {}", (omp_get_dynamic() ? "ON" : "OFF"));

  // Nested parallelism
  hmap::log::info("Nested parallelism: {}", (omp_get_nested() ? "ON" : "OFF"));

// Parallel region test
#pragma omp parallel
  {
#pragma omp single
    {
      hmap::log::info("Actual threads in parallel region: {}",
                      omp_get_num_threads());
    }
  }

#else
  hmap::log::info("OpenMP not enabled (no _OPENMP macro)");
#endif
}

} // namespace hmap
