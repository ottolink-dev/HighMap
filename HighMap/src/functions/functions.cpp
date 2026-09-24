/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/math/core.hpp"

namespace hmap
{

std::function<float(float, float)> make_xy_function_from_array(
    const Array     &array,
    const glm::vec4 &bbox)
{
  return [&array, &bbox](float x, float y) -> float
  {
    x = (x - bbox.x) / (bbox.y - bbox.x);
    y = (y - bbox.z) / (bbox.w - bbox.z);

    x = std::clamp(x, 0.f, 1.f);
    y = std::clamp(y, 0.f, 1.f);

    float xn = x * (array.shape.x - 1);
    float yn = y * (array.shape.y - 1);

    int   i = static_cast<int>(xn);
    int   j = static_cast<int>(yn);
    float u = xn - i;
    float v = yn - j;

    return array.get_value_bilinear_at(i, j, u, v);
  };
}

//----------------------------------------------------------------------
// base class
//----------------------------------------------------------------------

HMAP_FCT_XY_TYPE Function::get_delegate() const
{
  return this->delegate;
}

float Function::get_value(float x, float y, float ctrl_param) const
{
  return this->delegate(x, y, ctrl_param);
}

void Function::set_delegate(HMAP_FCT_XY_TYPE new_delegate)
{
  this->delegate = std::move(new_delegate);
}

//----------------------------------------------------------------------
// derived from Function class
//----------------------------------------------------------------------

ArrayFunction::ArrayFunction(hmap::Array array, glm::vec2 kw, bool periodic)
    : Function(), kw(kw), array(array)
{
  if (periodic)
    this->set_delegate(
        [this](float x, float y, float)
        {
          float xp = 0.5f * this->kw.x * x;
          float yp = 0.5f * this->kw.y * y;

          xp = 2.f * std::abs(xp - int(xp));
          yp = 2.f * std::abs(yp - int(yp));

          xp = xp < 1.f ? xp : 2.f - xp;
          yp = yp < 1.f ? yp : 2.f - yp;

          smoothstep5(xp);
          smoothstep5(yp);

          float xg = xp * (this->array.shape.x - 1);
          float yg = yp * (this->array.shape.y - 1);
          int   i = (int)xg;
          int   j = (int)yg;
          return this->array.get_value_bilinear_at(i, j, xg - i, yg - j);
        });
  else
    this->set_delegate(
        [this](float x, float y, float)
        {
          float xp = this->kw.x * x;
          float yp = this->kw.y * y;

          xp = xp < 0.f ? 0.f : (xp >= 1.f ? 1.f : xp - int(xp));
          yp = yp < 0.f ? 0.f : (yp >= 1.f ? 1.f : yp - int(yp));

          float xg = xp * (this->array.shape.x - 1);
          float yg = yp * (this->array.shape.y - 1);
          int   i = (int)xg;
          int   j = (int)yg;

          return this->array.get_value_bilinear_at(i, j, xg - i, yg - j);
        });
}

BiquadFunction::BiquadFunction(float gain, glm::vec2 center)
    : Function(), gain(gain), center(center)
{
  this->set_delegate(
      [this](float x, float y, float ctrl_param)
      {
        const float dx = x - this->center.x;
        const float dy = y - this->center.y;

        if (std::abs(dx) > 0.5f || std::abs(dy) > 0.5f) return 0.0f;

        float vx = 1.f - 4.f * dx * dx;
        float vy = 1.f - 4.f * dy * dy;
        float v = vx * vy;
        v = std::clamp(v, 0.f, 1.f);
        return std::pow(v, 1.f / (this->gain * ctrl_param));
      });
}

BumpFunction::BumpFunction(float gain, glm::vec2 center)
    : Function(), gain(gain), center(center)
{
  this->set_delegate(
      [this](float x, float y, float ctrl_param)
      {
        const float dx = x - this->center.x;
        const float dy = y - this->center.y;
        const float r2 = dx * dx + dy * dy;

        if (r2 > 0.25f) return 0.0f;

        const float denom = 1.0f - 4.0f * r2;
        const float exponent = -1.0f / denom;
        const float base = std::exp(exponent);
        const float power = 1.0f / (this->gain * ctrl_param);

        return std::pow(base / std::exp(-1.f), power);
      });
}

CraterFunction::CraterFunction(float     radius,
                               float     depth,
                               float     lip_decay,
                               float     lip_height_ratio,
                               glm::vec2 center)
    : Function(),
      radius(radius),
      depth(depth),
      lip_decay(lip_decay),
      lip_height_ratio(lip_height_ratio),
      center(center)
{
  this->set_delegate(
      [this](float x, float y, float ctrl_param)
      {
        float dx = x - this->center.x;
        float dy = y - this->center.y;
        float r = std::hypot(dx, dy);

        float value = std::min(
            r * r / (this->radius * this->radius),
            1.f + this->lip_height_ratio * ctrl_param *
                      std::exp(-(r - this->radius) / this->lip_decay));
        value -= 1.f;
        value *= this->depth;
        return value;
      });
}

DiskFunction::DiskFunction(float radius, float slope, glm::vec2 center)
    : Function(), radius(radius), slope(slope), center(center)
{
  this->set_delegate(
      [this](float x, float y, float ctrl_param)
      {
        float dx = x - this->center.x;
        float dy = y - this->center.y;
        float r = std::hypot(dx, dy);

        if (r < this->radius)
          return ctrl_param;
        else
        {
          float t = std::max(0.f, 1.f - this->slope * (r - this->radius));
          return ctrl_param * smoothstep3(t);
        }
      });
}

GaussianPulseFunction::GaussianPulseFunction(float sigma, glm::vec2 center)
    : Function(), center(center)
{
  this->set_sigma(sigma);
  this->set_delegate(
      [this](float x, float y, float ctrl_param)
      {
        float dx = x - this->center.x;
        float dy = y - this->center.y;
        float r2 = dx * dx + dy * dy;
        return std::exp(-0.5f * r2 * this->inv_sigma2 * ctrl_param);
      });
}

QuadSurfaceFunction::QuadSurfaceFunction(float c00,
                                         float c10,
                                         float c01,
                                         float c11)
    : Function(), c00(c00), c10(c10), c01(c01), c11(c11)
{
  this->set_delegate(
      [this](float x, float y, float ctrl_param)
      {
        float val = (1.f - x) * (1.f - y) * this->c00 +
                    x * (1.f - y) * this->c10 + (1.f - x) * y * this->c01 +
                    x * y * this->c11;
        return val * ctrl_param;
      });
}

RectangleFunction::RectangleFunction(float     rx,
                                     float     ry,
                                     float     angle,
                                     float     slope,
                                     glm::vec2 center)
    : Function(), rx(rx), ry(ry), slope(slope), center(center)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float ctrl_param)
      {
        x = x - this->center.x;
        y = y - this->center.y;

        float xc = this->ca * x + this->sa * y;
        float yc = -this->sa * x + this->ca * y;

        xc = std::abs(xc);
        yc = std::abs(yc);

        float ax;
        float ay;

        if (xc < this->rx)
          ax = ctrl_param;
        else
        {
          float t = std::max(0.f, 1.f - this->slope * (xc - this->rx));
          ax = ctrl_param * smoothstep3(t);
        }

        if (yc < this->ry)
          ay = ctrl_param;
        else
        {
          float t = std::max(0.f, 1.f - this->slope * (yc - this->ry));
          ay = ctrl_param * smoothstep3(t);
        }

        return ax * ay;
      });
}

RiftFunction::RiftFunction(float     angle,
                           float     slope,
                           float     width,
                           bool      sharp_bottom,
                           glm::vec2 center)
    : Function(), slope(slope), width(width), center(center)
{
  this->set_angle(angle);

  if (sharp_bottom)
    this->set_delegate(
        [this](float x, float y, float ctrl_param)
        {
          float local_width = 0.5f * this->width * ctrl_param;

          float r = this->ca * (x - this->center.x) +
                    this->sa * (y - this->center.y);
          r = std::abs(r);

          if (r > local_width + 1.f / this->slope)
            return 1.f;
          else if (r < local_width)
            return 0.f;
          else
          {
            r = (r - local_width) * this->slope;
            return smoothstep3_upper(r);
          }
        });
  else
    this->set_delegate(
        [this](float x, float y, float ctrl_param)
        {
          float local_width = 0.5f * this->width * ctrl_param;

          float r = this->ca * (x - this->center.x) +
                    this->sa * (y - this->center.y);
          r = std::abs(r);

          if (r > local_width + 1.f / this->slope)
            return 1.f;
          else if (r < local_width)
            return 0.f;
          else
          {
            r = (r - local_width) * this->slope;
            return smoothstep3(r);
          }
        });
}

SlopeFunction::SlopeFunction(float angle, float slope, glm::vec2 center)
    : Function(), slope(slope), center(center)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float ctrl_param)
      {
        float r = this->ca * (x - this->center.x) +
                  this->sa * (y - this->center.y);
        return this->slope * ctrl_param * r;
      });
}

StepFunction::StepFunction(float angle, float slope, glm::vec2 center)
    : Function(), slope(slope), center(center)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float ctrl_param)
      {
        float local_slope = this->slope * ctrl_param;

        float r = this->ca * (x - this->center.x) +
                  this->sa * (y - this->center.y);
        float dt = 0.5f / local_slope;
        if (r > dt)
          return 1.f;
        else if (r > -dt)
        {
          r = local_slope * (r + dt);
          return smoothstep3(r);
        }
        else
          return 0.f;
      });
}

WaveDuneFunction::WaveDuneFunction(glm::vec2 kw,
                                   float     angle,
                                   float     xtop,
                                   float     xbottom,
                                   float     phase_shift,
                                   glm::vec2 center)
    : Function(),
      kw(kw),
      xtop(xtop),
      xbottom(xbottom),
      phase_shift(phase_shift),
      center(center)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float)
      {
        float r = ca * this->kw.x * (x - this->center.x) +
                  sa * this->kw.y * (y - this->center.y);
        float xp = std::fmod(r + this->phase_shift +
                                 10.f * (this->kw.x + this->kw.y),
                             1.f);
        float yp = 0.f;

        if (xp < this->xtop)
        {
          float r = xp / this->xtop;
          yp = r * r * (3.f - 2.f * r);
        }
        else if (xp < this->xbottom)
        {
          float r = (xp - this->xbottom) / (this->xtop - this->xbottom);
          yp = r * r * (2.f - r);
        }
        return yp;
      });
}

WaveSineFunction::WaveSineFunction(glm::vec2 kw,
                                   float     angle,
                                   float     phase_shift,
                                   glm::vec2 center)
    : Function(), kw(kw), phase_shift(phase_shift), center(center)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float)
      {
        float r = ca * this->kw.x * (x - this->center.x) +
                  sa * this->kw.y * (y - this->center.y);
        return std::cos(2.f * M_PI * r + this->phase_shift);
      });
}

WaveSquareFunction::WaveSquareFunction(glm::vec2 kw,
                                       float     angle,
                                       float     phase_shift,
                                       glm::vec2 center)
    : Function(), kw(kw), phase_shift(phase_shift), center(center)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float)
      {
        float r = ca * this->kw.x * (x - this->center.x) +
                  sa * this->kw.y * (y - this->center.y) + this->phase_shift;
        return r = 2.f * std::floor(r) - std::floor(2.f * r) + 1.f;
      });
}

WaveTriangularFunction::WaveTriangularFunction(glm::vec2 kw,
                                               float     angle,
                                               float     slant_ratio,
                                               float     phase_shift,
                                               glm::vec2 center)
    : Function(),
      kw(kw),
      slant_ratio(slant_ratio),
      phase_shift(phase_shift),
      center(center)
{
  this->set_angle(angle);

  this->set_delegate(
      [this](float x, float y, float)
      {
        float r = ca * this->kw.x * (x - this->center.x) +
                  sa * this->kw.y * (y - this->center.y) + this->phase_shift;

        r = r - std::floor(r);
        if (r < this->slant_ratio)
          r /= this->slant_ratio;
        else
          r = 1.f - (r - this->slant_ratio) / (1.f - this->slant_ratio);
        return r * r * (3.f - 2.f * r);
      });
}

} // namespace hmap
