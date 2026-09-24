/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file functions.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Defines modular function objects for procedural generation, including
 * noise algorithms (Perlin, Simplex, Worley), geometric primitives (disks,
 * waves, slopes), and fractal compositions (FBM, IQ, Swiss). Provides parameter
 * control for frequency, seeding, and layer blending to construct complex
 * procedural patterns.
 *
 * @copyright Copyright (c) 2023
 */
#pragma once
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "FastNoiseLite.h"

#include "highmap/array.hpp"
#include "highmap/math/core.hpp"

/**
 * @typedef HMAP_FCT_XY_TYPE
 * @brief Type alias for a function taking three floats and returning a float.
 */
#define HMAP_FCT_XY_TYPE std::function<float(float, float, float)>

#define HMAP_GRADIENT_OFFSET 0.001f

namespace hmap
{

/**
 * @brief Create a continuous 2D function from a sampled array.
 *
 * This adaptor takes a discrete 2D array (indexed by i, j) defined over a
 * bounding box and returns a continuous function \f$f(x, y)\f$ that
 * interpolates values at arbitrary 2D positions inside the bounding box.
 *
 * @param  array Input array containing sampled values on a regular grid.
 * @param  bbox  Bounding box of the array in world coordinates, given as
 *             `(xmin, xmax, ymin, ymax)`.
 * @return       A function `f(point)` that returns the interpolated value at a
 *               2D position.
 */
std::function<float(float, float)> make_xy_function_from_array(
    const Array     &array,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @enum NoiseType
 * @brief Enumeration of various noise types used for procedural generation.
 *
 * This enumeration defines different types of noise algorithms that can be used
 * for procedural generation tasks such as terrain generation, texture
 * synthesis, and other applications where pseudo-random patterns are required.
 */
enum NoiseType : int
{
  PARBERRY,       ///< Parberry (Perlin variant)
  PERLIN,         ///< Perlin
  PERLIN_BILLOW,  ///< Perlin billow
  PERLIN_HALF,    ///< Perlin half
  SIMPLEX2,       ///< OpenSimplex2
  SIMPLEX2S,      ///< OpenSimplex2S
  VALUE,          ///< Value
  VALUE_CUBIC,    ///< Value (cubic)
  VALUE_DELAUNAY, ///< Value (delaunay)
  VALUE_LINEAR,   ///< Value (linear)
  WORLEY,         ///< Worley
  WORLEY_DOUBLE,  ///< Worley double
  WORLEY_VALUE,   ///< Worley (cell value return)
};

//----------------------------------------
// Base Function class and derived
//----------------------------------------

/**
 * @class Function
 * @brief A class that wraps a callable entity taking three floats and returning
 * a float.
 */
class Function
{
public:
  /**
   * @brief Default constructor. Initializes the delegate function to a default
   * that returns 0.
   */
  Function() : delegate([](float, float, float) { return 0.f; })
  {
  }

  /**
   * @brief Virtual destructor to ensure proper cleanup in derived classes.
   */
  virtual ~Function() = default;

  /**
   * @brief Constructor to initialize with a specific delegate function.
   * @param delegate The delegate function to initialize with.
   */
  explicit Function(HMAP_FCT_XY_TYPE delegate) : delegate(std::move(delegate))
  {
  }

  /**
   * @brief Get the current delegate function.
   * @return The current delegate function.
   */
  HMAP_FCT_XY_TYPE get_delegate() const;

  /**
   * @brief Call the delegate function with given arguments.
   * @param  x          The first float parameter.
   * @param  y          The second float parameter.
   * @param  ctrl_param The third float parameter.
   * @return            The result of the delegate function call.
   */
  float get_value(float x, float y, float ctrl_param) const;

  /**
   * @brief Set a new delegate function.
   * @param new_delegate The new delegate function to set.
   */
  void set_delegate(HMAP_FCT_XY_TYPE new_delegate);

private:
  HMAP_FCT_XY_TYPE delegate; ///< The stored delegate function object.
};

/**
 * @class ArrayFunction
 * @brief Array (x, y) function class.
 *
 * This class functions like an image sampler with normalized coordinates found
 * on a GPU. It is based on a 2D array, allowing interpolation of any value at a
 * given (x, y) coordinate. Specific boundary conditions, such as periodicity,
 * can also be used.
 */
class ArrayFunction : public Function
{
public:
  /**
   * @brief Construct a new ArrayFunction object.
   *
   * @param array    Data array.
   * @param kw       Noise wavenumbers {kx, ky} for each direction, with respect
   *                 to a unit domain.
   * @param periodic Whether the domain is periodic or not.
   */
  ArrayFunction(Array array, glm::vec2 kw, bool periodic = true);

  /**
   * @brief Set the array object.
   *
   * @param new_array New data array.
   */
  void set_array(Array new_array)
  {
    this->array = new_array;
  }

protected:
  glm::vec2 kw;  ///< Frequency scaling vector.
  bool periodic; ///< Flag indicating whether the domain is periodic or not.

private:
  Array array; ///< Data array representing the 2D field.
};

/**
 * @class BiquadFunction
 * @brief Biquad (x, y) function class.
 *
 * This class models a biquad function with a specified gain and reference
 * center. It extends the base Function class.
 */
class BiquadFunction : public Function
{
public:
  /**
   * @brief Construct a new Biquad Function object.
   *
   * @param gain   Gain that controls the steepness of the bump.
   * @param center Primitive reference center.
   */
  BiquadFunction(float gain, glm::vec2 center);

protected:
  float     gain;   ///< Gain value that controls the steepness of the bump.
  glm::vec2 center; ///< Primitive reference center.
};

/**
 * @class BumpFunction
 * @brief Bump (x, y) function class.
 *
 * This class models a bump function with a specified gain and reference center.
 * It extends the base Function class.
 */
class BumpFunction : public Function
{
public:
  /**
   * @brief Construct a new Bump Function object.
   *
   * @param gain   Gain that controls the steepness of the bump.
   * @param center Primitive reference center.
   */
  BumpFunction(float gain, glm::vec2 center);

protected:
  float     gain;   ///< Gain value that controls the steepness of the bump.
  glm::vec2 center; ///< Primitive reference center.
};

/**
 * @class CraterFunction
 * @brief Crater (x, y) function class.
 *
 * This class models a crater function with specified radius, depth, lip decay,
 * lip height ratio, and reference center. It extends the base Function class.
 */
class CraterFunction : public Function
{
public:
  /**
   * @brief Construct a new Crater Function object.
   *
   * @param radius           Radius of the crater (with respect to a unit
   * domain).
   * @param depth            Depth of the crater (with respect to a unit
   * domain).
   * @param lip_decay        Decay rate of the crater's lip (with respect to a
   *                         unit domain).
   * @param lip_height_ratio Height ratio of the crater's lip.
   * @param center           Primitive reference center.
   */
  CraterFunction(float     radius,
                 float     depth,
                 float     lip_decay,
                 float     lip_height_ratio,
                 glm::vec2 center);

protected:
  float     radius;           ///< Radius of the crater.
  float     depth;            ///< Depth of the crater.
  float     lip_decay;        ///< Decay rate of the crater's lip.
  float     lip_height_ratio; ///< Height ratio of the crater's lip.
  glm::vec2 center;           ///< Primitive reference center.
};

/**
 * @class DiskFunction
 * @brief DiskFunction (x, y) function class.
 */
class DiskFunction : public Function
{
public:
  /**
   * @brief Construct a new Disk Function object.
   *
   * @param radius Radius of the disk (with respect to a unit domain).
   * @param slope  Side slope (with respect to a unit domain).
   * @param center Primitive reference center.
   */
  DiskFunction(float radius, float slope, glm::vec2 center);

protected:
  float     radius; ///< Radius of the disk.
  float     slope;  ///< Slope of the disk.
  glm::vec2 center; ///< Primitive reference center.
};

/**
 * @class GaussianPulseFunction
 * @brief GaussianPulse (x, y) function class.
 *
 * This class models a gaussian pulse function. It extends the base Function
 * class.
 */
class GaussianPulseFunction : public Function
{
public:
  /**
   * @brief Construct a new GaussianPulse Function object.
   *
   * @param sigma  Pulse half-width (with respect to a unit domain).
   * @param center Primitive reference center.
   */
  GaussianPulseFunction(float sigma, glm::vec2 center);

  /**
   * @brief Set the half-width.
   *
   * @param new_sigma New half-width.
   */
  void set_sigma(float new_sigma)
  {
    this->sigma = new_sigma;
    this->inv_sigma2 = 1.f / (this->sigma * this->sigma);
  }

protected:
  float     sigma;  ///< Pulse half-width.
  glm::vec2 center; ///< Primitive reference center.

private:
  float inv_sigma2; ///< Cached squared inverse of the half-width.
};

/**
 * @class QuadSurfaceFunction
 * @brief Quadratic surface function defined by four corner values.
 *
 * This class interpolates a surface across the bounding box domain
 * using bilinear/quadratic interpolation from four corner elevations:
 * c00 (bottom-left), c10 (bottom-right), c01 (top-left), and c11 (top-right).
 */
class QuadSurfaceFunction : public Function
{
public:
  /**
   * @brief Construct a new QuadSurface Function object.
   *
   * @param c00 Value at bottom-left corner (x=xmin, y=ymin).
   * @param c10 Value at bottom-right corner (x=xmax, y=ymin).
   * @param c01 Value at top-left corner (x=xmin, y=ymax).
   * @param c11 Value at top-right corner (x=xmax, y=ymax).
   */
  QuadSurfaceFunction(float c00, float c10, float c01, float c11);

protected:
  float     c00;  ///< Bottom-left value.
  float     c10;  ///< Bottom-right value.
  float     c01;  ///< Top-left value.
  float     c11;  ///< Top-right value.
  glm::vec4 bbox; ///< Domain bounding box.
};

/**
 * @class RectangleFunction
 * @brief RectangleFunction (x, y) function class.
 */
class RectangleFunction : public Function
{
public:
  /**
   * @brief Construct a new Rectangle Function object.
   *
   * @param rx     Radius of the rectangle (with respect to a unit domain).
   * @param ry     Radius of the rectangle (with respect to a unit domain).
   * @param angle  Overall rotation angle (in degrees).
   * @param slope  Side slope (with respect to a unit domain).
   * @param center Primitive reference center.
   */
  RectangleFunction(float     rx,
                    float     ry,
                    float     angle,
                    float     slope,
                    glm::vec2 center);

  /**
   * @brief Set the angle.
   *
   * @param new_angle New angle in degrees.
   */
  void set_angle(float new_angle)
  {
    this->angle = new_angle;
    this->ca = std::cos(angle / 180.f * M_PI);
    this->sa = std::sin(angle / 180.f * M_PI);
  }

protected:
  float     rx, ry; ///< Radius of the rectangle.
  float     angle;  ///< Angle of the rectangle.
  float     slope;  ///< Slope of the rectangle.
  glm::vec2 center; ///< Primitive reference center.

private:
  float ca; ///< Cached cosine of the angle.
  float sa; ///< Cached sine of the angle.
};

/**
 * @class RiftFunction
 * @brief Rift (x, y) function class.
 *
 * This class models a step function with a specified angle, slope, and
 * reference center. It extends the base Function class.
 */
class RiftFunction : public Function
{
public:
  /**
   * @brief Construct a new Rift Function object.
   *
   * @param angle        Overall rotation angle (in degrees).
   * @param slope        Rift slope.
   * @param width        Rift width.
   * @param sharp_bottom Decide whether the rift bottom is sharp or not.
   * @param center       Primitive reference center.
   */
  RiftFunction(float     angle,
               float     slope,
               float     width,
               bool      sharp_bottom,
               glm::vec2 center);

  /**
   * @brief Set the angle.
   *
   * @param new_angle New angle in degrees.
   */
  void set_angle(float new_angle)
  {
    this->angle = new_angle;
    this->ca = std::cos(angle / 180.f * M_PI);
    this->sa = std::sin(angle / 180.f * M_PI);
  }

protected:
  float     angle;        ///< Overall rotation angle (in degrees).
  float     slope;        ///< Rift slope.
  float     width;        ///< Rift width.
  bool      sharp_bottom; ///< Rift bottom sharpness.
  glm::vec2 center;       ///< Primitive reference center.

private:
  float ca; ///< Cached cosine of the angle.
  float sa; ///< Cached sine of the angle.
};

/**
 * @class SlopeFunction
 * @brief Slope (x, y) function class.
 *
 * This class models a slope function with a specified angle, slope, and
 * reference center. It extends the base Function class.
 */
class SlopeFunction : public Function
{
public:
  /**
   * @brief Construct a new Slope Function object.
   *
   * @param angle  Overall rotation angle (in degrees).
   * @param slope  Step slope.
   * @param center Primitive reference center.
   */
  SlopeFunction(float angle, float slope, glm::vec2 center);

  /**
   * @brief Set the angle.
   *
   * @param new_angle New angle in degrees.
   */
  void set_angle(float new_angle)
  {
    this->angle = new_angle;
    this->ca = std::cos(angle / 180.f * M_PI);
    this->sa = std::sin(angle / 180.f * M_PI);
  }

protected:
  float     angle;  ///< Overall rotation angle (in degrees).
  float     slope;  ///< Step slope.
  glm::vec2 center; ///< Primitive reference center.

private:
  float ca; ///< Cached cosine of the angle.
  float sa; ///< Cached sine of the angle.
};

/**
 * @class StepFunction
 * @brief Step (x, y) function class.
 *
 * This class models a step function with a specified angle, slope, and
 * reference center. It extends the base Function class.
 */
class StepFunction : public Function
{
public:
  /**
   * @brief Construct a new Step Function object.
   *
   * @param angle  Overall rotation angle (in degrees).
   * @param slope  Step slope.
   * @param center Primitive reference center.
   */
  StepFunction(float angle, float slope, glm::vec2 center);

  /**
   * @brief Set the angle.
   *
   * @param new_angle New angle in degrees.
   */
  void set_angle(float new_angle)
  {
    this->angle = new_angle;
    this->ca = std::cos(angle / 180.f * M_PI);
    this->sa = std::sin(angle / 180.f * M_PI);
  }

protected:
  float     angle;  ///< Overall rotation angle (in degrees).
  float     slope;  ///< Step slope.
  glm::vec2 center; ///< Primitive reference center.

private:
  float ca; ///< Cached cosine of the angle.
  float sa; ///< Cached sine of the angle.
};

/**
 * @class WaveDuneFunction
 * @brief Wave dune (x, y) function class.
 *
 * This class models a wave dune function with specific noise wavenumbers,
 * rotation angle, slant ratio, and phase shift. It extends the base Function
 * class.
 */
class WaveDuneFunction : public Function
{
public:
  /**
   * @brief Construct a new Wave Dune Function object.
   *
   * @param kw          Noise wavenumbers {kx, ky} for each direction, with
   *                    respect to a unit domain.
   * @param angle       Overall rotation angle (in degrees).
   * @param xtop        Relative location of the top of the dune profile (in [0,
   *                    1]).
   * @param xbottom     Relative location of the foot of the dune profile (in
   * [0, 1]).
   * @param phase_shift Phase shift (in radians).
   * @param center      Primitive reference center.
   */
  WaveDuneFunction(glm::vec2 kw,
                   float     angle,
                   float     xtop,
                   float     xbottom,
                   float     phase_shift,
                   glm::vec2 center = {0.5f, 0.5f});

  /**
   * @brief Set the angle.
   *
   * @param new_angle New angle.
   */
  void set_angle(float new_angle)
  {
    this->angle = new_angle;
    this->ca = std::cos(angle / 180.f * M_PI);
    this->sa = std::sin(angle / 180.f * M_PI);
  }

protected:
  glm::vec2 kw;    ///< Frequency scaling vector.
  float     angle; ///< Overall rotation angle (in degrees).
  float xtop; ///< Relative location of the top of the dune profile (in [0, 1]).
  float xbottom; ///< Relative location of the foot of the dune profile (in [0,
  ///< 1]).
  float     phase_shift; ///< Phase shift (in radians).
  glm::vec2 center;      ///< Primitive reference center.

private:
  float ca; ///< Cached cosine of the angle.
  float sa; ///< Cached sine of the angle.
};

/**
 * @class WaveSineFunction
 * @brief Wave sine (x, y) function class.
 *
 * This class models a wave sine function with specific noise wavenumbers,
 * rotation angle, and phase shift. It extends the base Function class.
 */
class WaveSineFunction : public Function
{
public:
  /**
   * @brief Construct a new Wave Sine Function object.
   *
   * @param kw          Noise wavenumbers {kx, ky} for each direction, with
   *                    respect to a unit domain.
   * @param angle       Overall rotation angle (in degrees).
   * @param phase_shift Phase shift (in radians).
   * @param center      Primitive reference center.
   */
  WaveSineFunction(glm::vec2 kw,
                   float     angle,
                   float     phase_shift,
                   glm::vec2 center = {0.5f, 0.5f});

  /**
   * @brief Set the angle.
   *
   * @param new_angle New angle.
   */
  void set_angle(float new_angle)
  {
    this->angle = new_angle;
    this->ca = std::cos(angle / 180.f * M_PI);
    this->sa = std::sin(angle / 180.f * M_PI);
  }

protected:
  glm::vec2 kw;          ///< Frequency scaling vector.
  float     angle;       ///< Overall rotation angle (in degrees).
  float     phase_shift; ///< Phase shift (in radians).
  glm::vec2 center;      ///< Primitive reference center.

private:
  float ca; ///< Cached cosine of the angle.
  float sa; ///< Cached sine of the angle.
};

/**
 * @class WaveSquareFunction
 * @brief Wave square (x, y) function class.
 *
 * This class models a wave square function with specific noise wavenumbers,
 * rotation angle, and phase shift. It extends the base Function class.
 */
class WaveSquareFunction : public Function
{
public:
  /**
   * @brief Construct a new Wave Square Function object.
   *
   * @param kw          Noise wavenumbers {kx, ky} for each direction, with
   *                    respect to a unit domain.
   * @param angle       Overall rotation angle (in degrees).
   * @param phase_shift Phase shift (in radians).
   * @param center      Primitive reference center.
   */
  WaveSquareFunction(glm::vec2 kw,
                     float     angle,
                     float     phase_shift,
                     glm::vec2 center = {0.5f, 0.5f});

  /**
   * @brief Set the angle.
   *
   * @param new_angle New angle.
   */
  void set_angle(float new_angle)
  {
    this->angle = new_angle;
    this->ca = std::cos(angle / 180.f * M_PI);
    this->sa = std::sin(angle / 180.f * M_PI);
  }

protected:
  glm::vec2 kw;          ///< Frequency scaling vector.
  float     angle;       ///< Overall rotation angle (in degrees).
  float     phase_shift; ///< Phase shift (in radians).
  glm::vec2 center;      ///< Primitive reference center.

private:
  float ca; ///< Cached cosine of the angle.
  float sa; ///< Cached sine of the angle.
};

/**
 * @class WaveTriangularFunction
 * @brief Wave triangular (x, y) function class.
 *
 * This class models a wave triangular function with specific noise wavenumbers,
 * rotation angle, slant ratio, and phase shift. It extends the base Function
 * class.
 */
class WaveTriangularFunction : public Function
{
public:
  /**
   * @brief Construct a new Wave Triangular Function object.
   *
   * @param kw          Noise wavenumbers {kx, ky} for each direction, with
   *                    respect to a unit domain.
   * @param angle       Overall rotation angle (in degrees).
   * @param slant_ratio Relative location of the triangle apex, in [0, 1].
   * @param phase_shift Phase shift (in radians).
   * @param center      Primitive reference center.
   */
  WaveTriangularFunction(glm::vec2 kw,
                         float     angle,
                         float     slant_ratio,
                         float     phase_shift,
                         glm::vec2 center = {0.5f, 0.5f});

  /**
   * @brief Set the angle.
   *
   * @param new_angle New angle.
   */
  void set_angle(float new_angle)
  {
    this->angle = new_angle;
    this->ca = std::cos(angle / 180.f * M_PI);
    this->sa = std::sin(angle / 180.f * M_PI);
  }

protected:
  glm::vec2 kw;          ///< Frequency scaling vector.
  float     angle;       ///< Overall rotation angle (in degrees).
  float     slant_ratio; ///< Relative location of the triangle apex, in [0, 1].
  float     phase_shift; ///< Phase shift (in radians).
  glm::vec2 center;      ///< Primitive reference center.

private:
  float ca; ///< Cached cosine of the angle.
  float sa; ///< Cached sine of the angle.
};

//----------------------------------------
// NoiseFunction class and derived
//----------------------------------------

/**
 * @class NoiseFunction
 * @brief A class for generating noise functions.
 *
 * The `NoiseFunction` class provides an interface for generating noise-based
 * functions. It inherits from the `Function` class and adds parameters for
 * frequency scaling (`kw`) and a random seed (`seed`) to customize the noise
 * generation.
 */
class NoiseFunction : public Function
{
public:
  /**
   * @brief Default constructor. Initializes with default frequency scaling and
   * seed.
   */
  NoiseFunction() : Function(), kw(glm::vec2(0.f, 0.f)), seed(0)
  {
  }

  /**
   * @brief Constructor to initialize with specific frequency scaling.
   * @param kw Frequency scaling vector.
   */
  NoiseFunction(glm::vec2 kw) : Function(), kw(kw), seed(0)
  {
  }

  /**
   * @brief Constructor to initialize with specific frequency scaling and seed.
   * @param kw   Frequency scaling vector.
   * @param seed Random seed for noise generation.
   */
  NoiseFunction(glm::vec2 kw, std::uint32_t seed)
      : Function(), kw(kw), seed(seed)
  {
  }

  /**
   * @brief Get the frequency scaling vector.
   * @return The current frequency scaling vector.
   */
  glm::vec2 get_kw() const
  {
    return this->kw;
  }

  /**
   * @brief Get the random seed.
   * @return The current random seed.
   */
  std::uint32_t get_seed() const
  {
    return this->seed;
  }

  /**
   * @brief Set a new random seed for noise generation.
   * @param new_seed The new random seed.
   */
  virtual void set_seed(std::uint32_t new_seed)
  {
    this->seed = new_seed;
  }

  /**
   * @brief Set a new frequency scaling vector.
   * @param new_kw The new frequency scaling vector.
   */
  virtual void set_kw(glm::vec2 new_kw)
  {
    this->kw = new_kw;
  }

protected:
  glm::vec2     kw;   ///< Frequency scaling vector.
  std::uint32_t seed; ///< Random seed for noise generation.
};

//----------------------------------------
// Actual functions
//----------------------------------------

/**
 * @brief Parberry (x, y) function class.
 */
class ParberryFunction : public NoiseFunction
{
public:
  /**
   * @brief Gradient magnitude exponent.
   */
  float mu;

  /**
   * @brief Construct a new Perlin Function object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   * @param mu   Gradient magnitude exponent.
   */
  ParberryFunction(glm::vec2 kw, std::uint32_t seed, float mu);

  /**
   * @brief Initialize generator.
   */
  void initialize();

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t /*new_seed*/)
  {
    // FIX ME seed cannot be changed with current implemtantion.
  }

private:
  /**
   * @brief Perlin's B, a power of 2 usually equal to 256.
   */
  int perlin_b = 0X100;

  /**
   * @brief A bit mask, one less than B.
   */
  int perlin_bm = 0xff;

  /**
   * @brief Perlin's N.
   */
  int perlin_n = 0x100;

  /**
   * @brief Perlin's permutation table.
   */
  std::vector<int> p;

  /**
   * @brief Perlin's gradient table.
   */
  std::vector<std::vector<float>> g2;

  /**
   * @brief Ian Parberry's gradient magnitude table.
   */
  std::vector<float> m;
};

/**
 * @brief Perlin (x, y) function class.
 */
class PerlinFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new Perlin Function object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   */
  PerlinFunction(glm::vec2 kw, std::uint32_t seed);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise.SetSeed(new_seed);
  }

private:
  /**
   * @brief FastNoiseLite noise generator object.
   */
  FastNoiseLite noise;
};

/**
 * @brief Perlin 'billow' (x, y) function class.
 */
class PerlinBillowFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new Perlin Function object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   */
  PerlinBillowFunction(glm::vec2 kw, std::uint32_t seed);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise.SetSeed(new_seed);
  }

private:
  /**
   * @brief FastNoiseLite noise generator object.
   */
  FastNoiseLite noise;
};

/**
 * @brief Perlin 'half' (x, y) function class.
 */
class PerlinHalfFunction : public NoiseFunction
{
public:
  float k;

  /**
   * @brief Construct a new Perlin Function object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   * @param k    Smoothing factor.
   */
  PerlinHalfFunction(glm::vec2 kw, std::uint32_t seed, float k);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise.SetSeed(new_seed);
  }

private:
  /**
   * @brief FastNoiseLite noise generator object.
   */
  FastNoiseLite noise;
};

/**
 * @brief Perlin 'mix' (x, y) function class.
 */
class PerlinMixFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new Perlin Function object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   */
  PerlinMixFunction(glm::vec2 kw, std::uint32_t seed);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise.SetSeed(new_seed);
  }

private:
  /**
   * @brief FastNoiseLite noise generator object.
   */
  FastNoiseLite noise;
};

/**
 * @brief OpenSimplex2 (x, y) function class.
 */
class Simplex2Function : public NoiseFunction
{
public:
  /**
   * @brief Construct a new Perlin Function object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   */
  Simplex2Function(glm::vec2 kw, std::uint32_t seed);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise.SetSeed(new_seed);
  }

private:
  /**
   * @brief FastNoiseLite noise generator object.
   */
  FastNoiseLite noise;
};

/**
 * @brief OpenSimplex2S (x, y) function class.
 */
class Simplex2SFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new Perlin Function object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   */
  Simplex2SFunction(glm::vec2 kw, std::uint32_t seed);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise.SetSeed(new_seed);
  }

private:
  /**
   * @brief FastNoiseLite noise generator object.
   */
  FastNoiseLite noise;
};

/**
 * @brief Value noise (x, y) function class.
 */
class ValueNoiseFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new ValueNoiseFunction object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   */
  ValueNoiseFunction(glm::vec2 kw, std::uint32_t seed);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise.SetSeed(new_seed);
  }

private:
  /**
   * @brief FastNoiseLite noise generator object.
   */
  FastNoiseLite noise;
};

/**
 * @brief Value Cubic noise (x, y) function class.
 */
class ValueCubicNoiseFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new ValueCubicNoiseFunction object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   */
  ValueCubicNoiseFunction(glm::vec2 kw, std::uint32_t seed);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise.SetSeed(new_seed);
  }

private:
  /**
   * @brief FastNoiseLite noise generator object.
   */
  FastNoiseLite noise;
};

/**
 * @brief ValueDelaunayNoise (x, y) function class.
 */
class ValueDelaunayNoiseFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new ValueNoiseFunction object.
   *
   * @param kw   Noise wavenumber, with respect to a unit domain.
   * @param seed Random seed number.
   */
  ValueDelaunayNoiseFunction(glm::vec2 kw, std::uint32_t seed);

  /**
   * @brief Set the wavenumber attribute.
   *
   * @param new_kw New kw.
   */
  void set_kw(glm::vec2 new_kw)
  {
    NoiseFunction::set_kw(new_kw);
    this->update_interpolation_function();
  }

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->update_interpolation_function();
  }

  /**
   * @brief Update base interpolation.
   */
  void update_interpolation_function();
};

/**
 * @brief ValueLinearNoiseFunction (x, y) function class.
 */
class ValueLinearNoiseFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new ValueLinearNoiseFunction object.
   *
   * @param kw   Noise wavenumber, with respect to a unit domain.
   * @param seed Random seed number.
   */
  ValueLinearNoiseFunction(glm::vec2 kw, std::uint32_t seed);

  /**
   * @brief Set the wavenumber attribute.
   *
   * @param new_kw New kw.
   */
  void set_kw(glm::vec2 new_kw)
  {
    NoiseFunction::set_kw(new_kw);
    this->update_interpolation_function();
  }

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->update_interpolation_function();
  }

  /**
   * @brief Update base interpolation.
   */
  void update_interpolation_function();
};

/**
 * @brief Worley (x, y) function class.
 */
class WorleyFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new Worley Function object.
   *
   * @param kw   Noise wavenumbers {kx, ky} for each directions, with respect to
   * a unit domain.
   * @param seed Random seed number.
   */
  WorleyFunction(glm::vec2     kw,
                 std::uint32_t seed,
                 bool          return_cell_value = false);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise.SetSeed(new_seed);
  }

private:
  /**
   * @brief FastNoiseLite noise generator object.
   */
  FastNoiseLite noise;
};

/**
 * @brief Worley (x, y) function class.
 */
class WorleyDoubleFunction : public NoiseFunction
{
public:
  /**
   * @brief Amplitude ratio between each Worley noise.
   */
  float ratio;

  /**
   * @brief Transition smoothing parameter.
   */
  float k;

  /**
   * @brief Construct a new Worley Double Function object.
   *
   * @param kw    Noise wavenumbers {kx, ky} for each directions, with respect
   * to a unit domain.
   * @param seed  Random seed number.
   * @param ratio Amplitude ratio between each Worley noise.
   * @param k     Transition smoothing parameter.
   */
  WorleyDoubleFunction(glm::vec2 kw, std::uint32_t seed, float ratio, float k);

  /**
   * @brief Set the seed attribute.
   *
   * @param new_seed New seed.
   */
  void set_seed(std::uint32_t new_seed)
  {
    NoiseFunction::set_seed(new_seed);
    this->noise1.SetSeed(new_seed);
    this->noise2.SetSeed(new_seed + 1);
  }

private:
  /**
   * @brief FastNoiseLite noise generator objects.
   */
  FastNoiseLite noise1, noise2;
};

//----------------------------------------
// Composite functions
//----------------------------------------

/**
 * @class GenericFractalFunction
 * @brief A class for generating fractal noise functions based on an underlying
 * noise function.
 *
 * The `GenericFractalFunction` class generates fractal noise using an
 * underlying base noise function. It allows customization of the fractal
 * properties such as octaves, weight, persistence, and lacunarity.
 */
class GenericFractalFunction : public NoiseFunction
{
public:
  /**
   * @brief Construct a new GenericFractalFunction object.
   *
   * @param p_base      Unique pointer to the base noise function.
   * @param octaves     Number of octaves in the fractal noise.
   * @param weight      Weight of the base noise function.
   * @param persistence Persistence of the fractal noise.
   * @param lacunarity  Lacunarity of the fractal noise.
   */
  explicit GenericFractalFunction(std::unique_ptr<NoiseFunction> p_base,
                                  int                            octaves,
                                  float                          weight,
                                  float                          persistence,
                                  float                          lacunarity);

  /**
   * @brief Set the frequency scaling vector.
   *
   * @param new_kw The new frequency scaling vector.
   */
  void set_kw(glm::vec2 new_kw) override
  {
    NoiseFunction::set_kw(new_kw);
    this->p_base->set_kw(new_kw);
  }

  /**
   * @brief Set the lacunarity of the fractal noise.
   *
   * @param new_lacunarity The new lacunarity value.
   */
  void set_lacunarity(float new_lacunarity)
  {
    this->lacunarity = new_lacunarity;
  }

  /**
   * @brief Set the number of octaves in the fractal noise.
   *
   * @param new_octaves The new number of octaves.
   */
  void set_octaves(int new_octaves)
  {
    this->octaves = new_octaves;
    this->update_amp0();
  }

  /**
   * @brief Set the persistence of the fractal noise.
   *
   * @param new_persistence The new persistence value.
   */
  void set_persistence(float new_persistence)
  {
    this->persistence = new_persistence;
    this->update_amp0();
  }

  /**
   * @brief Set a new random seed for the noise generation.
   *
   * @param new_seed The new random seed.
   */
  void set_seed(std::uint32_t new_seed) override
  {
    NoiseFunction::set_seed(new_seed);
    this->p_base->set_seed(new_seed);
  }

  /**
   * @brief Scale the initial amplitude of the fractal noise.
   *
   * @param scale The scale factor for the initial amplitude.
   */
  void scale_amp0(float scale)
  {
    this->amp0 *= scale;
  }

  /**
   * @brief Get the lacunarity of the fractal noise.
   *
   * @return The current lacunarity value.
   */
  float get_lacunarity() const
  {
    return this->lacunarity;
  }

  /**
   * @brief Get the number of octaves in the fractal noise.
   *
   * @return The current number of octaves.
   */
  int get_octaves() const
  {
    return this->octaves;
  }

  /**
   * @brief Get the persistence of the fractal noise.
   *
   * @return The current persistence value.
   */
  float get_persistence() const
  {
    return this->persistence;
  }

  /**
   * @brief Get the weight of the fractal noise.
   *
   * @return The current weight value.
   */
  float get_weight() const
  {
    return this->weight;
  }

protected:
  /**
   * @brief Update the initial amplitude (amp0) based on the current octaves and
   * persistence.
   */
  void update_amp0();

protected:
  std::unique_ptr<NoiseFunction>
        p_base;      ///< Unique pointer to the base noise function.
  int   octaves;     ///< Number of octaves in the fractal noise.
  float weight;      ///< Weight of the base noise function.
  float persistence; ///< Persistence of the fractal noise.
  float lacunarity;  ///< Lacunarity of the fractal noise.
  float amp0;        ///< Initial amplitude of the fractal noise.
};

/**
 * @class FbmFunction
 * @brief Fractional Brownian Motion (FBM) function class.
 */
class FbmFunction : public GenericFractalFunction
{
public:
  /**
   * @brief Construct a new Fbm Function object.
   *
   * @param p_base      Unique pointer to the base noise function.
   * @param octaves     Number of octaves.
   * @param weight      Octave weighting.
   * @param persistence Octave persistence.
   * @param lacunarity  Defines the wavenumber ratio between each octave.
   */
  FbmFunction(std::unique_ptr<NoiseFunction> p_base,
              int                            octaves,
              float                          weight,
              float                          persistence,
              float                          lacunarity);
};
/**
 * @class FbmIqFunction
 * @brief IQ layering function class.
 */
class FbmIqFunction : public GenericFractalFunction
{
public:
  /**
   * @brief Construct a new Fbm Iq Function object.
   *
   * @param p_base         Unique pointer to the base noise function.
   * @param octaves        Number of octaves.
   * @param weight         Octave weighting.
   * @param persistence    Octave persistence.
   * @param lacunarity     Defines the wavenumber ratio between each octave.
   * @param gradient_scale Gradient scale influence.
   */
  FbmIqFunction(std::unique_ptr<NoiseFunction> p_base,
                int                            octaves,
                float                          weight,
                float                          persistence,
                float                          lacunarity,
                float                          gradient_scale);

  /**
   * @brief Set the gradient scale.
   *
   * @param new_gradient_scale New gradient scale.
   */
  void set_gradient_scale(float new_gradient_scale)
  {
    this->gradient_scale = new_gradient_scale;
  }

protected:
  float gradient_scale; ///< Gradient scale influence.
};

/**
 * @class FbmJordanFunction
 * @brief Jordan layering function class.
 */
class FbmJordanFunction : public GenericFractalFunction
{
public:
  /**
   * @brief Construct a new Fbm Jordan Function object.
   *
   * @param p_base      Unique pointer to the base noise function.
   * @param octaves     Number of octaves.
   * @param weight      Octave weighting.
   * @param persistence Octave persistence.
   * @param lacunarity  Defines the wavenumber ratio between each octave.
   * @param warp0       Initial warp.
   * @param damp0       Initial damp.
   * @param warp_scale  Warp scale.
   * @param damp_scale  Damp scale.
   */
  FbmJordanFunction(std::unique_ptr<NoiseFunction> p_base,
                    int                            octaves,
                    float                          weight,
                    float                          persistence,
                    float                          lacunarity,
                    float                          warp0,
                    float                          damp0,
                    float                          warp_scale,
                    float                          damp_scale);

  /**
   * @brief Set the initial warp.
   *
   * @param new_warp0 New initial warp.
   */
  void set_warp0(float new_warp0)
  {
    this->warp0 = new_warp0;
  }

  /**
   * @brief Set the initial damp.
   *
   * @param new_damp0 New initial damp.
   */
  void set_damp0(float new_damp0)
  {
    this->damp0 = new_damp0;
  }

  /**
   * @brief Set the warp scale.
   *
   * @param new_warp_scale New warp scale.
   */
  void set_warp_scale(float new_warp_scale)
  {
    this->warp_scale = new_warp_scale;
  }

  /**
   * @brief Set the damp scale.
   *
   * @param new_damp_scale New damp scale.
   */
  void set_damp_scale(float new_damp_scale)
  {
    this->damp_scale = new_damp_scale;
  }

protected:
  float warp0;      ///< Initial warp.
  float damp0;      ///< Initial damp.
  float warp_scale; ///< Warp scale.
  float damp_scale; ///< Damp scale.
};

/**
 * @class FbmPingpongFunction
 * @brief Pingpong layering function class.
 */
class FbmPingpongFunction : public GenericFractalFunction
{
public:
  /**
   * @brief Construct a new Fbm Pingpong Function object.
   *
   * @param p_base      Unique pointer to the base noise function.
   * @param octaves     Number of octaves.
   * @param weight      Octave weighting.
   * @param persistence Octave persistence.
   * @param lacunarity  Defines the wavenumber ratio between each octave.
   */
  FbmPingpongFunction(std::unique_ptr<NoiseFunction> p_base,
                      int                            octaves,
                      float                          weight,
                      float                          persistence,
                      float                          lacunarity);

  /**
   * @brief Set the smoothing parameter.
   *
   * @param new_k_smoothing New smoothing parameter.
   */
  void set_k_smoothing(float new_k_smoothing)
  {
    this->k_smoothing = new_k_smoothing;
  }

protected:
  float k_smoothing; ///< Smoothing parameter.
};

/**
 * @class FbmRidgedFunction
 * @brief Ridged layering function class.
 */
class FbmRidgedFunction : public GenericFractalFunction
{
public:
  /**
   * @brief Construct a new Fbm Ridged Function object.
   *
   * @param p_base      Unique pointer to the base noise function.
   * @param octaves     Number of octaves.
   * @param weight      Octave weighting.
   * @param persistence Octave persistence.
   * @param lacunarity  Defines the wavenumber ratio between each octave.
   * @param k_smoothing Smoothing parameter.
   */
  FbmRidgedFunction(std::unique_ptr<NoiseFunction> p_base,
                    int                            octaves,
                    float                          weight,
                    float                          persistence,
                    float                          lacunarity,
                    float                          k_smoothing);

  /**
   * @brief Set the smoothing parameter.
   *
   * @param new_k_smoothing New smoothing parameter.
   */
  void set_k_smoothing(float new_k_smoothing)
  {
    this->k_smoothing = new_k_smoothing;
  }

protected:
  float k_smoothing; ///< Smoothing parameter.
};

/**
 * @class FbmSwissFunction
 * @brief Swiss layering function class.
 */
class FbmSwissFunction : public GenericFractalFunction
{
public:
  /**
   * @brief Construct a new Fbm Swiss Function object.
   *
   * @param p_base      Unique pointer to the base noise function.
   * @param octaves     Number of octaves.
   * @param weight      Octave weighting.
   * @param persistence Octave persistence.
   * @param lacunarity  Defines the wavenumber ratio between each octave.
   * @param warp_scale  Warping scale.
   */
  FbmSwissFunction(std::unique_ptr<NoiseFunction> p_base,
                   int                            octaves,
                   float                          weight,
                   float                          persistence,
                   float                          lacunarity,
                   float                          warp_scale);

  /**
   * @brief Set the warp scale.
   *
   * @param new_warp_scale New warp scale.
   */
  void set_warp_scale(float new_warp_scale)
  {
    this->warp_scale = new_warp_scale;
    this->warp_scale_normalized = new_warp_scale / this->kw.x;
  }

protected:
  float warp_scale;            ///< Warping scale.
  float warp_scale_normalized; ///< Normalized warping scale.
};

//----------------------------------------
// Field functions
//----------------------------------------

/**
 * @class FieldFunction
 * @brief Field function class.
 */
class FieldFunction : public Function
{
public:
  /**
   * @brief Construct a new FieldFunction object.
   *
   * @param p_base Unique pointer to the base function.
   */
  FieldFunction(std::unique_ptr<Function> p_base);

  /**
   * @brief Construct a new FieldFunction object with given vectors.
   *
   * @param p_base Unique pointer to the base function.
   * @param xr     Vector of x coordinates representing the centers of the
   *               primitive.
   * @param yr     Vector of y coordinates representing the centers of the
   *               primitive.
   * @param zr     Vector of z coordinates used to scale the primitive in x and
   * y directions, and also to scale the primitive amplitude if requested.
   */
  FieldFunction(std::unique_ptr<Function> p_base,
                std::vector<float>        xr,
                std::vector<float>        yr,
                std::vector<float>        zr);

  /**
   * @brief Set the x coordinates representing the centers of the primitive.
   *
   * @param new_xr New vector of x coordinates.
   */
  void set_xr(std::vector<float> new_xr)
  {
    this->xr = new_xr;
  }

  /**
   * @brief Set the y coordinates representing the centers of the primitive.
   *
   * @param new_yr New vector of y coordinates.
   */
  void set_yr(std::vector<float> new_yr)
  {
    this->yr = new_yr;
  }

  /**
   * @brief Set the z coordinates used to scale the primitive in x and y
   * directions, and also to scale the primitive amplitude if requested.
   *
   * @param new_zr New vector of z coordinates.
   */
  void set_zr(std::vector<float> new_zr)
  {
    this->zr = new_zr;
  }

protected:
  std::vector<float> xr; ///< Vector of x coordinates representing the centers
                         ///< of the primitive.
  std::vector<float> yr; ///< Vector of y coordinates representing the centers
                         ///< of the primitive.
  std::vector<float>
      zr; ///< Vector of z coordinates used to scale the primitive in x and y
          ///< directions, and to scale the primitive amplitude if requested.

private:
  std::unique_ptr<Function> p_base; ///< Unique pointer to the base function.

  /**
   * @brief Setup the delegate function.
   */
  void setup_delegate();
};

//----------------------------------------
// Helpers
//----------------------------------------

/**
 * @brief Create a noise function based on the specified noise type.
 *
 * This function creates an instance of a noise function corresponding to the
 * specified `noise_type`. It initializes the noise function with the given
 * frequency scaling vector (`kw`) and random seed (`seed`).
 *
 * @param  noise_type The type of noise function to create (e.g., PERLIN,
 *                    SIMPLEX2, etc.).
 * @param  kw         The frequency scaling vector to be used for the noise
 *                    function.
 * @param  seed       The random seed for noise generation.
 * @return            A `std::unique_ptr` to the created noise function.
 */
std::unique_ptr<hmap::NoiseFunction> create_noise_function_from_type(
    NoiseType     noise_type,
    glm::vec2     kw,
    std::uint32_t seed);

} // namespace hmap
