/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file openmp.hpp
 */
#pragma once

namespace hmap
{

/**
 * @brief Initialize OpenMP and set the number of threads to use.
 *
 * @param num_threads Number of threads; pass 0 (default) to use every
 *                    available processor (`omp_get_num_procs()`).
 * @return true if OpenMP is enabled, false otherwise.
 */
bool init_openmp(int num_threads = 0);

void log_openmp_info();

} // namespace hmap