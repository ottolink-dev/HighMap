/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file openmp.hpp
 */
#pragma once

namespace hmap
{

bool init_openmp(int num_threads = -1);

void log_openmp_info();

} // namespace hmap