/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <utility>
#include <vector>

#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"

extern "C"
{
#include "config.h"
//
#include "nn.h"
//
#include "nncommon.h"
//
#include "delaunay.h"
}

namespace hmap
{

NaturalNeighborInterpolator::~NaturalNeighborInterpolator()
{
  if (this->handle) nnai_destroy(this->handle);
  if (this->d) delaunay_destroy(this->d);
}

NaturalNeighborInterpolator::NaturalNeighborInterpolator(
    NaturalNeighborInterpolator &&other) noexcept
    : handle(other.handle),
      d(other.d),
      xout(std::move(other.xout)),
      yout(std::move(other.yout)),
      nout(other.nout),
      nin(other.nin)
{
  other.handle = nullptr;
  other.d = nullptr;
  other.nout = 0;
  other.nin = 0;
}

NaturalNeighborInterpolator &NaturalNeighborInterpolator::operator=(
    NaturalNeighborInterpolator &&other) noexcept
{
  if (this != &other)
  {
    if (this->handle) nnai_destroy(this->handle);
    if (this->d) delaunay_destroy(this->d);

    this->handle = other.handle;
    this->d = other.d;
    this->xout = std::move(other.xout);
    this->yout = std::move(other.yout);
    this->nout = other.nout;
    this->nin = other.nin;

    other.handle = nullptr;
    other.d = nullptr;
    other.nout = 0;
    other.nin = 0;
  }
  return *this;
}

void NaturalNeighborInterpolator::build(const std::vector<float> &xin,
                                        const std::vector<float> &yin)
{
  if (this->handle)
  {
    nnai_destroy(this->handle);
    this->handle = nullptr;
  }
  if (this->d)
  {
    delaunay_destroy(this->d);
    this->d = nullptr;
  }
  this->nin = 0;

  if (!validate_min_size(xin, 3, "Input points xin") ||
      xin.size() != yin.size() || this->nout == 0)
    return;

  int npoints = static_cast<int>(xin.size());

  // create delaunay triangulation
  std::vector<point> pin(npoints);
  for (int i = 0; i < npoints; ++i)
  {
    pin[i].x = double(xin[i]);
    pin[i].y = double(yin[i]);
    pin[i].z = double(0.0); // can be set later via values array
  }

  this->d = delaunay_build(npoints, pin.data(), 0, nullptr, 0, nullptr);
  if (!this->d) return;

  this->nin = xin.size();

  // create nnai interpolator
  this->handle = nnai_build(this->d,
                            int(this->nout),
                            this->xout.data(),
                            this->yout.data());
  if (this->handle)
  {
    nnai_setwmin(this->handle, double(-1e30));
  }
}

void NaturalNeighborInterpolator::interpolate(
    const std::vector<float> &values_in,
    std::vector<float>       &values_out) const
{
  if (!this->handle || !this->d) return;
  if (values_in.size() != this->nin) return;

  values_out.resize(this->nout);

  // convert float input to double
  std::vector<double> values_in_d(values_in.begin(), values_in.end());
  std::vector<double> values_out_d(this->nout);

  // call the nnai C function
  nnai_interpolate(this->handle, values_in_d.data(), values_out_d.data());

  // convert double output back to float
  for (size_t k = 0; k < this->nout; ++k)
    values_out[k] = static_cast<float>(values_out_d[k]);
}

void NaturalNeighborInterpolator::setup_output_points(
    const std::vector<float> &x,
    const std::vector<float> &y)
{
  if (x.size() != y.size()) return;

  this->xout.clear();
  this->yout.clear();

  this->xout.reserve(x.size());
  this->yout.reserve(y.size());

  for (const auto &v : x)
    this->xout.push_back(double(v));

  for (const auto &v : y)
    this->yout.push_back(double(v));

  this->nout = this->xout.size();
}

void NaturalNeighborInterpolator::setup_output_points(
    const std::vector<double> &x,
    const std::vector<double> &y)
{
  if (x.size() != y.size()) return;

  this->xout = x;
  this->yout = y;
  this->nout = this->xout.size();
}

void NaturalNeighborInterpolator::setup_output_points(std::vector<double> &&x,
                                                      std::vector<double> &&y)
{
  if (x.size() != y.size()) return;

  this->xout = std::move(x);
  this->yout = std::move(y);
  this->nout = this->xout.size();
}

} // namespace hmap
