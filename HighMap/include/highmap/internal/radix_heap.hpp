/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file radix_heap.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Monotone priority queue (radix heap) and IEEE-754 float-to-uint32
 * order-preserving key helpers.
 * @copyright Copyright (c) 2023
 */
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

#include <bit>

namespace hmap
{

// order-preserving mapping between IEEE-754 float and uint32_t for radix
// sorting
inline uint32_t float_key(float f)
{
  uint32_t u = 0;
  std::memcpy(&u, &f, 4);
  return (u & 0x80000000u) ? ~u : (u | 0x80000000u);
}

inline float key_float(uint32_t k)
{
  uint32_t u = (k & 0x80000000u) ? (k & 0x7fffffffu) : ~k;
  float    f = 0.f;
  std::memcpy(&f, &u, 4);
  return f;
}

inline float next_up(float f)
{
  return key_float(float_key(f) + 1u);
}

// monotone priority queue (radix heap) for fast exact flood / Dijkstra
template <typename T = int> struct RadixHeap
{
  std::vector<std::pair<uint32_t, T>> bucket[33];
  uint32_t                            last = 0;
  size_t                              count = 0;

  bool empty() const
  {
    return count == 0;
  }

  size_t size() const
  {
    return count;
  }

  void clear()
  {
    for (auto &b : bucket)
      b.clear();

    last = 0;
    count = 0;
  }

  static int slot(uint32_t x)
  {
    return (int)std::bit_width(x);
  }

  void push(uint32_t key, const T &val)
  {
    bucket[slot(key ^ last)].emplace_back(key, val);
    ++count;
  }

  void push(uint32_t key, T &&val)
  {
    bucket[slot(key ^ last)].emplace_back(key, std::move(val));
    ++count;
  }

  std::pair<uint32_t, T> pop()
  {
    if (bucket[0].empty())
    {
      int i = 1;
      while (bucket[i].empty())
        ++i;

      uint32_t mn = 0xffffffffu;
      for (const auto &e : bucket[i])
        mn = std::min(mn, e.first);

      last = mn;
      for (auto &e : bucket[i])
        bucket[slot(e.first ^ last)].push_back(std::move(e));

      bucket[i].clear();
    }

    auto e = std::move(bucket[0].back());
    bucket[0].pop_back();
    --count;
    return e;
  }
};

} // namespace hmap
