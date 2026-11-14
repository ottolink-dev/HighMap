/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <random>
#include <utility>
#include <vector>

#include "highmap/functions.hpp"
#include "highmap/math/core.hpp"

namespace hmap
{

// methods

ParberryFunction::ParberryFunction(glm::vec2 kw, std::uint32_t seed, float mu)
    : NoiseFunction(kw, seed), mu(mu)
{
  this->p.resize(this->perlin_b);
  this->g2.resize(this->perlin_b, std::vector<float>(2));
  this->m.resize(this->perlin_b);

  this->initialize();
}

void ParberryFunction::initialize()
{
  std::mt19937                          gen(this->seed);
  std::uniform_real_distribution<float> dis(0.f, 1.f);

  // --- initialize the permutation, gradient, and magnitude tables

  // random gradient vectors
  for (int i = 0; i < this->perlin_b; i++)
  {
    g2[i][0] = 2.f * dis(gen) - 1.f;
    g2[i][1] = 2.f * dis(gen) - 1.f;

    // normalize
    float s = std::hypot(g2[i][0], g2[i][1]);
    g2[i][0] /= s;
    g2[i][1] /= s;
  }

  // random permutation
  for (int i = 0; i < this->perlin_b; i++)
    p[i] = i;

  for (int i = this->perlin_b - 1; i > 0; i--)
    std::swap(p[i], p[(int)(RAND_MAX * dis(gen)) % (i + 1)]);

  // gradient magnitude array
  float s = 1.f;

  for (int i = 0; i < this->perlin_b; i++)
  {
    m[i] = s;
    s /= this->mu;
  }

  // --- noise function
}

float ParberryFunction::operator()(float x, float y, float) const
{
  int   bx0, bx1, by0, by1;
  float rx0, rx1, ry0, ry1;

  {
    float t = this->kw.x * x + this->perlin_n;
    bx0 = ((int)t) & this->perlin_bm;
    bx1 = (bx0 + 1) & this->perlin_bm;
    rx0 = t - (int)t;
    rx1 = rx0 - 1.f;
  }

  {
    float t = this->kw.y * y + this->perlin_n;
    by0 = ((int)t) & this->perlin_bm;
    by1 = (by0 + 1) & this->perlin_bm;
    ry0 = t - (int)t;
    ry1 = ry0 - 1.f;
  }

  int b00 = this->p[(this->p[bx0] + by0) & this->perlin_bm];
  int b10 = this->p[(this->p[bx1] + by0) & this->perlin_bm];
  int b01 = this->p[(this->p[bx0] + by1) & this->perlin_bm];
  int b11 = this->p[(this->p[bx1] + by1) & this->perlin_bm];

  float sx = smoothstep5(rx0);

  float u, v;
  u = this->m[b00] * (rx0 * this->g2[b00][0] + ry0 * this->g2[b00][1]);
  v = this->m[b10] * (rx1 * this->g2[b10][0] + ry0 * this->g2[b10][1]);
  float a = lerp(u, v, sx);

  u = this->m[b01] * (rx0 * this->g2[b01][0] + ry1 * this->g2[b01][1]);
  v = this->m[b11] * (rx1 * this->g2[b11][0] + ry1 * this->g2[b11][1]);
  float b = lerp(u, v, sx);

  float sy = smoothstep5(ry0);

  return lerp(a, b, sy);
}

} // namespace hmap
