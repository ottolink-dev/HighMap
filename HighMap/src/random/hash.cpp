/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "highmap/random.hpp"

namespace hmap
{

float fast_hash32_to_unit_float(unsigned int seed, size_t k)
{
  uint32_t x = static_cast<uint32_t>(k) ^ seed;
  x ^= x >> 16;
  x *= 0x45d9f3bu;
  x ^= x >> 16;
  return static_cast<float>((x >> 8) * (1.f / float(1 << 24)));
}

float hash01(uint32_t a, uint32_t b)
{
  return (float)(hash32(a ^ hash32(b + 0x9e3779b9u)) >> 8) * (1.f / 16777216.f);
}

uint32_t hash32(uint32_t x)
{
  x ^= x >> 16;
  x *= 0x7feb352dU;
  x ^= x >> 15;
  x *= 0x846ca68bU;
  x ^= x >> 16;
  return x;
}

uint64_t splitmix64(uint64_t x)
{
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

float splitmix64_to_unit_float(unsigned int seed, size_t k)
{
  // combine seed + index into 64-bit to avoid truncation
  uint64_t x = static_cast<uint64_t>(seed) ^
               (static_cast<uint64_t>(k) + 0x9e3779b97f4a7c15ull);

  // 64-bit mix (SplitMix64-inspired)
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ull;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebull;
  x ^= x >> 31;

  // convert to float in [0,1)
  return static_cast<float>((x >> 40) & 0xFFFFFF) / static_cast<float>(1 << 24);
}

float uniform01(uint64_t h)
{
  return (float)(h >> 40) / 16777216.f; // 24-bit mantissa in [0, 1)
}

float vnoise(float x, float y, uint32_t seed)
{
  float fx = std::floor(x), fy = std::floor(y);
  int   xi = (int)fx, yi = (int)fy;
  float tx = x - fx, ty = y - fy;
  tx = tx * tx * (3.f - 2.f * tx);
  ty = ty * ty * (3.f - 2.f * ty);
  auto h = [&](int i, int j)
  {
    return 2.f * hash01((uint32_t)i * 0x9E3779B1u + (uint32_t)j * 0x85EBCA77u,
                        seed) -
           1.f;
  };
  float a = h(xi, yi), b = h(xi + 1, yi), c = h(xi, yi + 1),
        d = h(xi + 1, yi + 1);
  return (a + tx * (b - a)) + ty * ((c + tx * (d - c)) - (a + tx * (b - a)));
}

} // namespace hmap
