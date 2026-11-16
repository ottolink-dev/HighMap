/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file assert.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
#include <iostream>

#include "highmap/array.hpp"

namespace hmap
{

struct AssertResults
{
  std::string msg = "";
  float       diff = 0.f;
  float       tolerance = 0.f;
  float       count = 0.f;
  bool        ret = true;

  void print() const
  {
    std::string prefix = ret ? "ok  - " : "NOK - ";
    std::cout << prefix << "diff: " << std::to_string(diff)
              << ", tolerance:" << std::to_string(tolerance)
              << ", count:" << std::to_string(count) << "\n";
  }
};

bool assert_almost_equal(const Array       &a,
                         const Array       &b,
                         float              tolerance = 1e-6f,
                         const std::string &fname = "",
                         AssertResults     *p_results = nullptr);

// true if all non-zero pixels of a are also non-zero in b
bool is_subset(const Array &a, const Array &b);

void print_array(const Array &arr,
                 bool         show_coords = false,
                 int          width = 6,
                 int          precision = 2);

void dump_visual_check(const Array       &array,
                       const std::string &fname = "check.png");

} // namespace hmap