/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace hmap
{

//----------------------------------------
// GenericFractalFunction & Fbm*Function
//----------------------------------------

template <typename BaseFn>
GenericFractalFunction<BaseFn>::GenericFractalFunction(BaseFn base,
                                                       int    octaves,
                                                       float  weight,
                                                       float  persistence,
                                                       float  lacunarity)
    : NoiseFunction(base.get_kw(), base.get_seed()),
      base(std::move(base)),
      octaves(octaves),
      weight(weight),
      persistence(persistence),
      lacunarity(lacunarity)
{
  this->update_amp0();
}

template <typename BaseFn> void GenericFractalFunction<BaseFn>::update_amp0()
{
  float sum = 0.f;
  float amp = 1.f;
  for (int k = 0; k < this->octaves; k++)
  {
    sum += amp;
    amp *= this->persistence;
  }
  this->amp0 = 1.f / sum;
}

template <typename BaseFn>
FbmFunction<BaseFn>::FbmFunction(BaseFn base,
                                 int    octaves,
                                 float  weight,
                                 float  persistence,
                                 float  lacunarity)
    : GenericFractalFunction<BaseFn>(std::move(base),
                                     octaves,
                                     weight,
                                     persistence,
                                     lacunarity)
{
}

template <typename BaseFn>
float FbmFunction<BaseFn>::operator()(float x,
                                      float y,
                                      float /* ctrl_param */) const
{
  float sum = 0.f;
  float amp = this->amp0;
  float ki = 1.f;
  float kj = 1.f;

  for (int k = 0; k < this->octaves; k++)
  {
    float value = this->base(ki * x, kj * y, 0.f);
    sum += value * amp;
    amp *= (1.f - this->weight) +
           this->weight * std::clamp(value + 1.f, 0.f, 1.f);
    amp *= this->persistence;
    ki *= this->lacunarity;
    kj *= this->lacunarity;
  }

  return sum;
}

template <typename BaseFn>
FbmIqFunction<BaseFn>::FbmIqFunction(BaseFn base,
                                     int    octaves,
                                     float  weight,
                                     float  persistence,
                                     float  lacunarity,
                                     float  gradient_scale)
    : GenericFractalFunction<BaseFn>(std::move(base),
                                     octaves,
                                     weight,
                                     persistence,
                                     lacunarity),
      gradient_scale(gradient_scale)
{
}

template <typename BaseFn>
float FbmIqFunction<BaseFn>::operator()(float x,
                                        float y,
                                        float /* ctrl_param */) const
{
  float sum = 0.f;
  float amp = this->amp0;
  float ki = 1.f;
  float kj = 1.f;
  float dx_sum = 0.f;
  float dy_sum = 0.f;

  for (int k = 0; k < this->octaves; k++)
  {
    float xw = ki * x;
    float yw = kj * y;

    float value = this->base(xw, yw, 0.f);
    float dx = (this->base(xw + HMAP_GRADIENT_OFFSET, yw, 0.f) -
                this->base(xw - HMAP_GRADIENT_OFFSET, yw, 0.f)) /
               (2.f * HMAP_GRADIENT_OFFSET);
    float dy = (this->base(xw, yw + HMAP_GRADIENT_OFFSET, 0.f) -
                this->base(xw, yw - HMAP_GRADIENT_OFFSET, 0.f)) /
               (2.f * HMAP_GRADIENT_OFFSET);

    dx_sum += dx;
    dy_sum += dy;

    sum += value * amp /
           (1.f + this->gradient_scale * (dx_sum * dx_sum + dy_sum * dy_sum));
    amp *= (1.f - this->weight) +
           this->weight * std::clamp(value + 1.f, 0.f, 1.f);
    amp *= this->persistence;
    ki *= this->lacunarity;
    kj *= this->lacunarity;
  }

  return sum;
}

template <typename BaseFn>
FbmJordanFunction<BaseFn>::FbmJordanFunction(BaseFn base,
                                             int    octaves,
                                             float  weight,
                                             float  persistence,
                                             float  lacunarity,
                                             float  warp0,
                                             float  damp0,
                                             float  warp_scale,
                                             float  damp_scale)
    : GenericFractalFunction<BaseFn>(std::move(base),
                                     octaves,
                                     weight,
                                     persistence,
                                     lacunarity),
      warp0(warp0),
      damp0(damp0),
      warp_scale(warp_scale),
      damp_scale(damp_scale)
{
}

template <typename BaseFn>
float FbmJordanFunction<BaseFn>::operator()(float x,
                                            float y,
                                            float /* ctrl_param */) const
{
  float value = this->base(x, y, 0.f);
  float dx = (this->base(x + HMAP_GRADIENT_OFFSET, y, 0.f) -
              this->base(x - HMAP_GRADIENT_OFFSET, y, 0.f)) /
             (2.f * HMAP_GRADIENT_OFFSET);
  float dy = (this->base(x, y + HMAP_GRADIENT_OFFSET, 0.f) -
              this->base(x, y - HMAP_GRADIENT_OFFSET, 0.f)) /
             (2.f * HMAP_GRADIENT_OFFSET);

  float sum = value * value;
  float dsum_warp_x = this->warp0 * dx;
  float dsum_warp_y = this->warp0 * dy;
  float dsum_damp_x = this->damp0 * dx;
  float dsum_damp_y = this->damp0 * dy;

  float amp = this->persistence;
  float ki = this->lacunarity;
  float kj = this->lacunarity;
  float damped_amp = amp * this->persistence;

  for (int k = 1; k < this->octaves; k++)
  {
    float xw = ki * x + dsum_warp_x;
    float yw = kj * y + dsum_warp_y;

    float v = this->base(xw, yw, 0.f);
    float dkx = (this->base(xw + HMAP_GRADIENT_OFFSET, yw, 0.f) -
                 this->base(xw - HMAP_GRADIENT_OFFSET, yw, 0.f)) /
                (2.f * HMAP_GRADIENT_OFFSET);
    float dky = (this->base(xw, yw + HMAP_GRADIENT_OFFSET, 0.f) -
                 this->base(xw, yw - HMAP_GRADIENT_OFFSET, 0.f)) /
                (2.f * HMAP_GRADIENT_OFFSET);

    sum += damped_amp * v * v;
    dsum_warp_x += this->warp_scale * dkx;
    dsum_warp_y += this->warp_scale * dky;
    dsum_damp_x += this->damp_scale * dkx;
    dsum_damp_y += this->damp_scale * dky;

    amp *= this->persistence;
    damped_amp = amp *
                 (1.f - this->damp_scale / (1.f + dsum_damp_x * dsum_damp_x +
                                            dsum_damp_y * dsum_damp_y));
    ki *= this->lacunarity;
    kj *= this->lacunarity;
  }

  return sum;
}

template <typename BaseFn>
FbmPingpongFunction<BaseFn>::FbmPingpongFunction(BaseFn base,
                                                 int    octaves,
                                                 float  weight,
                                                 float  persistence,
                                                 float  lacunarity)
    : GenericFractalFunction<BaseFn>(std::move(base),
                                     octaves,
                                     weight,
                                     persistence,
                                     lacunarity)
{
}

template <typename BaseFn>
float FbmPingpongFunction<BaseFn>::operator()(float x,
                                              float y,
                                              float /* ctrl_param */) const
{
  float sum = 0.f;
  float amp = this->amp0;
  float ki = 1.f;
  float kj = 1.f;

  for (int k = 0; k < this->octaves; k++)
  {
    float value = (this->base(ki * x, kj * y, 0.f) + 1.f) * 0.5f;
    value = 1.f - std::abs(std::fmod(value * 2.f, 2.f) - 1.f);
    value = 2.f * value - 1.f;
    sum += value * amp;
    amp *= (1.f - this->weight) +
           this->weight * std::clamp(value + 1.f, 0.f, 1.f);
    amp *= this->persistence;
    ki *= this->lacunarity;
    kj *= this->lacunarity;
  }

  return sum;
}

template <typename BaseFn>
FbmRidgedFunction<BaseFn>::FbmRidgedFunction(BaseFn base,
                                             int    octaves,
                                             float  weight,
                                             float  persistence,
                                             float  lacunarity,
                                             float  k_smoothing)
    : GenericFractalFunction<BaseFn>(std::move(base),
                                     octaves,
                                     weight,
                                     persistence,
                                     lacunarity),
      k_smoothing(k_smoothing)
{
}

template <typename BaseFn>
float FbmRidgedFunction<BaseFn>::operator()(float x,
                                            float y,
                                            float /* ctrl_param */) const
{
  float sum = 0.f;
  float amp = this->amp0;
  float ki = 1.f;
  float kj = 1.f;

  for (int k = 0; k < this->octaves; k++)
  {
    float value = 1.f - std::abs(this->base(ki * x, kj * y, 0.f));
    value = value * value;
    sum += value * amp;
    amp *= (1.f - this->weight) + this->weight * std::clamp(value, 0.f, 1.f);
    amp *= this->persistence;
    ki *= this->lacunarity;
    kj *= this->lacunarity;
  }

  return sum;
}

template <typename BaseFn>
FbmSwissFunction<BaseFn>::FbmSwissFunction(BaseFn base,
                                           int    octaves,
                                           float  weight,
                                           float  persistence,
                                           float  lacunarity,
                                           float  warp_scale)
    : GenericFractalFunction<BaseFn>(std::move(base),
                                     octaves,
                                     weight,
                                     persistence,
                                     lacunarity),
      warp_scale(warp_scale)
{
}

template <typename BaseFn>
float FbmSwissFunction<BaseFn>::operator()(float x,
                                           float y,
                                           float /* ctrl_param */) const
{
  float sum = 0.f;
  float amp = this->amp0;
  float ki = 1.f;
  float kj = 1.f;
  float dx_sum = 0.f;
  float dy_sum = 0.f;

  for (int k = 0; k < this->octaves; k++)
  {
    float xw = ki * x + this->warp_scale * dx_sum;
    float yw = kj * y + this->warp_scale * dy_sum;

    float value = this->base(xw, yw, 0.f);
    float dx = (this->base(xw + HMAP_GRADIENT_OFFSET, yw, 0.f) -
                this->base(xw - HMAP_GRADIENT_OFFSET, yw, 0.f)) /
               (2.f * HMAP_GRADIENT_OFFSET);
    float dy = (this->base(xw, yw + HMAP_GRADIENT_OFFSET, 0.f) -
                this->base(xw, yw - HMAP_GRADIENT_OFFSET, 0.f)) /
               (2.f * HMAP_GRADIENT_OFFSET);

    dx_sum += amp * dx;
    dy_sum += amp * dy;

    sum += value * amp;
    amp *= (1.f - this->weight) +
           this->weight * std::clamp(value + 1.f, 0.f, 1.f);
    amp *= this->persistence;
    ki *= this->lacunarity;
    kj *= this->lacunarity;
  }

  return sum;
}

//----------------------------------------
// FieldFunction
//----------------------------------------

template <typename BaseFn>
FieldFunction<BaseFn>::FieldFunction(BaseFn base)
    : xr({0.5f}), yr({0.5f}), zr({1.f}), base(std::move(base))
{
}

template <typename BaseFn>
FieldFunction<BaseFn>::FieldFunction(BaseFn             base,
                                     std::vector<float> xr,
                                     std::vector<float> yr,
                                     std::vector<float> zr)
    : xr(std::move(xr)),
      yr(std::move(yr)),
      zr(std::move(zr)),
      base(std::move(base))
{
}

template <typename BaseFn>
float FieldFunction<BaseFn>::operator()(float x,
                                        float y,
                                        float ctrl_param) const
{
  float value = 0.f;
  for (size_t k = 0; k < this->xr.size(); k++)
  {
    float xs = (x - this->xr[k]) * this->zr[k];
    float ys = (y - this->yr[k]) * this->zr[k];
    value += this->base(xs, ys, ctrl_param) / (1.f + this->zr[k]);
  }
  return value;
}

//----------------------------------------
// Helpers
//----------------------------------------

template <typename Fn>
decltype(auto) with_noise_function(NoiseType     noise_type,
                                   glm::vec2     kw,
                                   std::uint32_t seed,
                                   Fn          &&fn)
{
  switch (noise_type)
  {
  case NoiseType::PERLIN: return fn(PerlinFunction(kw, seed));
  case NoiseType::PERLIN_BILLOW: return fn(PerlinBillowFunction(kw, seed));
  case NoiseType::PERLIN_HALF: return fn(PerlinHalfFunction(kw, seed, 0.f));
  case NoiseType::SIMPLEX2: return fn(Simplex2Function(kw, seed));
  case NoiseType::SIMPLEX2S: return fn(Simplex2SFunction(kw, seed));
  case NoiseType::VALUE: return fn(ValueNoiseFunction(kw, seed));
  case NoiseType::VALUE_CUBIC: return fn(ValueCubicNoiseFunction(kw, seed));
  case NoiseType::VALUE_DELAUNAY:
    return fn(ValueDelaunayNoiseFunction(kw, seed));
  case NoiseType::VALUE_LINEAR: return fn(ValueLinearNoiseFunction(kw, seed));
  case NoiseType::WORLEY: return fn(WorleyFunction(kw, seed));
  case NoiseType::WORLEY_DOUBLE:
    return fn(WorleyDoubleFunction(kw, seed, 0.5f, 0.f));
  default:
    throw std::invalid_argument("with_noise_function: unsupported NoiseType");
  }
}

} // namespace hmap
