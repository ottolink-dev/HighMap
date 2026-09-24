/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file functions.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2023
 */
#pragma once

#include "highmap/array.hpp"
#include "highmap/erosion.hpp"

namespace hmap
{

/**
 * @enum PrimitiveType
 * @brief Defines the primitive shape used for synthesis.
 */
enum PrimitiveType : int
{
  PRIM_BIQUAD_PULSE,
  PRIM_BUMP,
  PRIM_CONE,
  PRIM_CONE_SMOOTH,
  PRIM_CUBIC_PULSE,
  PRIM_SMOOTH_COSINE,
};

/**
 * @brief Generates a primitive shape as a 2D array.
 *
 * @param  primitive_type Type of primitive to generate.
 * @param  shape          Output array resolution.
 * @param  p_noise_x      Optional X-direction noise modulation.
 * @param  p_noise_y      Optional Y-direction noise modulation.
 * @param  center         Center position of the primitive (normalized).
 * @param  bbox           Bounding box of the primitive (normalized).
 * @return                The generated primitive as an Array.
 */
Array get_primitive_base(const PrimitiveType &primitive_type,
                         const glm::ivec2    &shape,
                         const Array         *p_noise_x = nullptr,
                         const Array         *p_noise_y = nullptr,
                         glm::vec2            center = {0.5f, 0.5f},
                         glm::vec4            bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generate a bounded band (capsule) mask.
 *
 * Creates a 0 to 1 mask shaped as a capsule: the value is 1 on a core segment
 * defined by `center`, `angle` and `length`, and falls off to 0 at a distance
 * `width` from the segment, following the requested radial profile. Typical
 * uses: mountain-range envelopes, island-chain arc envelopes and, with a length
 * larger than the domain, latitude bands / polar caps (place
 * `center` on or beyond a domain edge).
 *
 * @param  shape          Output array shape.
 * @param  angle          Band orientation angle in degrees.
 * @param  length         Core segment length in normalized domain units. A zero
 *                        length degenerates to a radial bump.
 * @param  width          Falloff half-width in normalized domain units.
 * @param  profile        Radial profile used for the falloff.
 * @param  profile_param  Additional parameter for the profile.
 * @param  p_noise_r      Optional noise array used to locally perturb the
 *                        falloff distance (edge displacement).
 * @param  p_noise_offset Optional noise array used to locally offset the band
 *                        axis position (normalized by `width`).
 * @param  center         Band center in normalized coordinates.
 * @param  bbox           Domain bounding box.
 * @return                Generated band mask in [0, 1].
 *
 * **Example**
 * @include ex_band.cpp
 *
 * **Result**
 * @image html ex_band.png
 */
Array band(glm::ivec2    shape,
           float         angle,
           float         length = 0.5f,
           float         width = 0.1f,
           RadialProfile profile = RadialProfile::RP_SMOOTHSTEP,
           float         profile_param = 0.f,
           const Array  *p_noise_r = nullptr,
           const Array  *p_noise_offset = nullptr,
           glm::vec2     center = {0.5f, 0.5f},
           glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a 'biquadratic pulse'.
 *
 * @param  shape                Array shape.
 * @param  gain                 Gain (the higher, the steeper).
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the gain parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Perlin billow noise.
 *
 * **Example**
 * @include ex_biquad_pulse.cpp
 *
 * **Result**
 * @image html ex_biquad_pulse.png
 */
Array biquad_pulse(glm::ivec2   shape,
                   float        gain = 1.f,
                   const Array *p_ctrl_param = nullptr,
                   const Array *p_noise_x = nullptr,
                   const Array *p_noise_y = nullptr,
                   glm::vec2    center = {0.5f, 0.5f},
                   glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

Array biquad_pulse_x(glm::ivec2 shape, glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f});
Array biquad_pulse_y(glm::ivec2 shape, glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a bump.
 *
 * @param  shape                Array shape.
 * @param  gain                 Gain (the higher, the steeper the bump).
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the gain parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  center               Primitive reference center.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Perlin billow noise.
 *
 * **Example**
 * @include ex_bump.cpp
 *
 * **Result**
 * @image html ex_bump.png
 */
Array bump(glm::ivec2   shape,
           float        gain = 1.f,
           const Array *p_ctrl_param = nullptr, // gain multiplier
           const Array *p_noise_x = nullptr,
           const Array *p_noise_y = nullptr,
           glm::vec2    center = {0.5f, 0.5f},
           glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a 2D Lorentzian bump pattern.
 *
 * This function fills an array with values computed from a normalized
 * Lorentzian-shaped bump function. The bump is centered at a given position,
 * has a specified radius, and can vary in width depending on a control
 * parameter. Optional noise arrays can be applied to perturb the x and y
 * coordinates before evaluation.
 *
 * The Lorentzian function is normalized so that it returns 1.0 at the center
 * and smoothly decays to 0.0 at the radius boundary. Outside the radius, the
 * value is exactly 0.
 *
 * @param  shape        Dimensions of the output array (width, height).
 * @param  width_factor Scaling factor controlling the bump width relative to
 *                      the control parameter.
 * @param  radius       Radius of the bump in the coordinate space of @p bbox.
 * @param  p_ctrl_param Optional pointer to an array of per-pixel control
 *                      parameters. If provided, it modulates the width of the
 *                      Lorentzian bump.
 * @param  p_noise_x    Optional pointer to an array of x-offset noise values.
 *                      Values are added to the x coordinate before evaluation.
 * @param  p_noise_y    Optional pointer to an array of y-offset noise values.
 *                      Values are added to the y coordinate before evaluation.
 * @param  center       Center of the bump in the coordinate space of @p bbox.
 * @param  bbox         Bounding box (xmin, ymin, xmax, ymax) defining the
 *                      coordinate space for the array.
 *
 * @return              A new Array object of size @p shape containing the bump
 *                      pattern.
 *
 * @note The function normalizes the Lorentzian curve so that the maximum value
 * is 1 at the center and smoothly decays to 0 at the boundary defined by
 *       @p radius.
 *
 * @see                 fill_array_using_xy_function
 *
 * **Example**
 * @include ex_bump.cpp
 *
 * **Result**
 * @image html ex_bump.png
 */
Array bump_lorentzian(
    glm::ivec2   shape,
    float        shape_factor = 0.5f,
    float        radius = 0.5f,
    const Array *p_ctrl_param = nullptr, // shape_factor multiplier
    const Array *p_noise_x = nullptr,
    const Array *p_noise_y = nullptr,
    glm::vec2    center = {0.5f, 0.5f},
    glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a checkerboard heightmap.
 *
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumber with respect to a unit domain.
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array New array.
 *
 * **Example**
 * @include ex_checkerboard.cpp
 *
 * **Result**
 * @image html ex_checkerboard.png
 */
Array checkerboard(glm::ivec2   shape,
                   glm::vec2    kw,
                   const Array *p_noise_x = nullptr,
                   const Array *p_noise_y = nullptr,
                   glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a synthetic conical mountain heightmap.
 *
 * This function creates a simple cone-shaped elevation field centered at a
 * given position. The cone profile follows a linear slope until it reaches the
 * talus angle, resembling volcanic cones or alluvial fans. Optional
 * displacement noise fields can be applied to perturb the cone shape.
 *
 * @param  shape     Output array shape (resolution in x and y).
 * @param  slope     Slope angle of the cone (controls steepness of the sides).
 * @param  center    Center of the cone in normalized coordinates (default =
 *                   {0.5f, 0.5f}).
 * @param  p_noise_x Optional pointer to external displacement noise field
 *                   (X-axis).
 * @param  p_noise_y Optional pointer to external displacement noise field
 *                   (Y-axis).
 * @param  bbox      Bounding box of the generation domain in normalized
 *                   coordinates (default = {0.f, 1.f, 0.f, 1.f}).
 *
 * @return           Array containing the generated conical heightmap.
 *
 * **Example**
 * @include ex_cone.cpp
 *
 * **Result**
 * @image html ex_cone.png
 */
Array cone(glm::ivec2   shape,
           float        slope,
           float        apex_elevation = 1.f,
           bool         smooth_profile = false,
           glm::vec2    center = {0.5f, 0.5f},
           const Array *p_noise_x = nullptr,
           const Array *p_noise_y = nullptr,
           glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a complex conical heightfield with valleys, directional
 * bias, and radial waviness.
 *
 * This function procedurally constructs an island-like cone shape inside a 2D
 * array. The height at each position is determined from the normalized radial
 * distance to the center, with optional modulation by valleys, directional
 * bias, and radial waviness. An erosion profile is applied to shape the
 * valleys, and an optional control map can modulate the valley amplitude
 * locally. The resulting height values are clamped to [0, 1].
 *
 * @param  shape               Output array dimensions (width, height).
 * @param  alpha               Exponent controlling the steepness of the cone
 *                             slope.
 * @param  radius              Effective radius of the island (in coordinate
 *                             space units).
 * @param  valley_amp          Global amplitude of valley depressions.
 * @param  valley_nb           Number of valleys around the island.
 * @param  valley_decay_ratio  Controls how valley amplitude decays toward the
 *                             center.
 * @param  valley_angle0       Angular offset of the first valley (in degrees).
 * @param  erosion_profile     Erosion profile type used for shaping valley
 *                             cross-sections.
 * @param  erosion_delta       Smoothing parameter for the erosion profile
 *                             function.
 * @param  radial_waviness_amp Amplitude of radial sinusoidal perturbations
 *                             (coastal waviness).
 * @param  radial_waviness_kw  Frequency multiplier for radial waviness.
 * @param  bias_angle          Direction (in degrees) of the slope bias (e.g.
 *                             wind or sun exposure).
 * @param  bias_amp            Amplitude of the directional bias effect.
 * @param  bias_exponent       Controls how bias strength varies with radius.
 * @param  center              Center position of the cone in coordinate space.
 * @param  p_ctrl_param        Optional control array; its values modulate the
 *                             local valley amplitude. If nullptr, the valley
 *                             amplitude is uniform.
 * @param  p_noise_x           Optional X-displacement noise field for
 *                             coordinate perturbation.
 * @param  p_noise_y           Optional Y-displacement noise field for
 *                             coordinate perturbation.
 * @param  bbox                Bounding box defining the world coordinates
 *                             covered by the array.
 *
 * @return                     A 2D Array containing the generated heightfield
 *                             values in the range [0, 1].
 *
 * **Example**
 * @include ex_cone.cpp
 *
 * **Result**
 * @image html ex_cone.png
 */

Array cone_complex(
    glm::ivec2            shape,
    float                 alpha,
    float                 radius = 0.5f,
    bool                  smooth_profile = true,
    float                 valley_amp = 0.2f,
    int                   valley_nb = 5,
    float                 valley_decay_ratio = 0.5f,
    float                 valley_angle0 = 15.f,
    const ErosionProfile &erosion_profile = ErosionProfile::EP_TRIANGLE_GRENIER,
    float                 erosion_delta = 0.01f,
    float                 radial_waviness_amp = 0.05f,
    float                 radial_waviness_kw = 2.f,
    float                 bias_angle = 30.f,
    float                 bias_amp = 0.75f,
    float                 bias_exponent = 1.f,
    glm::vec2             center = {0.5f, 0.5f},
    const Array          *p_ctrl_param = nullptr,
    const Array          *p_noise_x = nullptr,
    const Array          *p_noise_y = nullptr,
    glm::vec4             bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a smooth conical heightmap using a sigmoid-based profile.
 *
 * This function creates a 2D array representing a cone shape centered at a
 * given position. The height decreases radially from the apex according to a
 * sigmoid-shaped falloff, producing a smooth, rounded cone rather than a sharp
 * linear one. Optionally, displacement noise can be applied along the X and Y
 * axes to perturb the cone surface.
 *
 * @param  shape     Dimensions of the output array (width, height).
 * @param  alpha     Peak elevation value of the cone (controls maximum height).
 * @param  radius    Radius of the cone base, expressed in normalized [0, 1]
 *                   units relative to the bounding box.
 * @param  center    Normalized coordinates of the cone’s center (default =
 *                   {0.5, 0.5}).
 * @param  p_noise_x Optional pointer to an array representing horizontal
 *                   displacement noise (nullptr for none).
 * @param  p_noise_y Optional pointer to an array representing vertical
 *                   displacement noise (nullptr for none).
 * @param  bbox      Bounding box of the generated region in (xmin, xmax, ymin,
 *                   ymax) order.
 *
 * @return           A 2D Array object containing the generated conical
 *                   heightmap.
 *
 * **Example**
 * @include ex_cone.cpp
 *
 * **Result**
 * @image html ex_cone.png
 * */
Array cone_sigmoid(glm::ivec2   shape,
                   float        alpha,
                   float        radius = 0.5f,
                   glm::vec2    center = {0.5f, 0.5f},
                   const Array *p_noise_x = nullptr,
                   const Array *p_noise_y = nullptr,
                   glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a constant value array.
 *
 * @param  shape Array shape.
 * @param  value Filling value.
 * @return       Array New array.
 */
Array constant(glm::ivec2 shape, float value = 0.f);

/**
 * @brief Generates a cubic pulse array.
 */
Array cubic_pulse(glm::ivec2   shape,
                  const Array *p_noise_x,
                  const Array *p_noise_y,
                  glm::vec2    center = {0.5f, 0.5f},
                  glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a disk-shaped heightmap with optional modifications.
 *
 * This function creates a 2D array representing a disk shape with a specified
 * radius, slope, and other optional parameters such as control parameters,
 * noise, and stretching for additional customization.
 *
 * @param  shape        Dimensions of the output array (width, height).
 * @param  radius       Radius of the disk, in normalized coordinates (0.0 to
 *                      1.0).
 * @param  slope        Slope of the disk edge transition. A larger value makes
 *                      the edge transition sharper. Defaults to 1.0.
 * @param  p_ctrl_param Optional pointer to an `Array` controlling custom
 *                      parameters for the disk generation.
 * @param  p_noise_x    Optional pointer to an `Array` for adding noise in the
 *                      x-direction.
 * @param  p_noise_y    Optional pointer to an `Array` for adding noise in the
 *                      y-direction.
 * @param  p_stretching Optional pointer to an `Array` for stretching the disk
 *                      horizontally or vertically.
 * @param  center       Center of the disk in normalized coordinates (0.0 to
 *                      1.0). Defaults to {0.5, 0.5}.
 * @param  bbox         Bounding box for the disk in normalized coordinates
 *                     {x_min, x_max, y_min, y_max}. Defaults to {0.0, 1.0, 0.0,
 * 1.0}.
 *
 * @return              A 2D array representing the generated disk shape.
 *
 * * **Example**
 * @include ex_disk.cpp
 *
 * **Result**
 * @image html ex_disk.png
 */
Array disk(glm::ivec2   shape,
           float        radius,
           float        slope = 1.f,
           const Array *p_ctrl_param = nullptr,
           const Array *p_noise_x = nullptr,
           const Array *p_noise_y = nullptr,
           glm::vec2    center = {0.5f, 0.5f},
           glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a gaussian_decay pulse kernel.
 *
 * @param  shape                Array shape.
 * @param  sigma                Gaussian sigma (in pixels).
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the half-width parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  center               Primitive reference center.
 * @param  center               Primitive reference center.
 * @param  bbox                 Domain bounding box.
 * @return                      Array
 *
 * **Example**
 * @include ex_gaussian_pulse.cpp
 *
 * **Result**
 * @image html ex_gaussian_pulse.png
 */
Array gaussian_pulse(glm::ivec2   shape,
                     float        sigma,
                     const Array *p_ctrl_param = nullptr,
                     const Array *p_noise_x = nullptr,
                     const Array *p_noise_y = nullptr,
                     glm::vec2    center = {0.5f, 0.5f},
                     glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generate a multi-step height profile along a rotated axis.
 *
 * Constructs an Array where values follow a sequence of geometric steps blended
 * with an optional linear transition. Steps are rotated by the given angle, can
 * be shaped with an exponent, and may be modulated by optional control and
 * noise fields.
 *
 * @param  shape              Output array resolution.
 * @param  angle              Rotation angle in degrees.
 * @param  r                  Geometric ratio between successive step widths.
 * @param  nsteps             Number of steps.
 * @param  elevation_exponent Exponent shaping step heights.
 * @param  shape_gain         Gain used to reshape intra-step transitions.
 * @param  scale              Scale of the step axis.
 * @param  outer_slope        Slope outside the [0,1] axis interval.
 * @param  p_ctrl_param       Optional control parameter array.
 * @param  p_noise_x          Optional x-axis noise array.
 * @param  p_noise_y          Optional y-axis noise array.
 * @param  center             Center of the step axis.
 * @param  bbox               Bounding box of the domain.
 * @return                    Generated Array.
 *
 * **Example**
 * @include ex_multisteps.cpp
 *
 * **Result**
 * @image html ex_multisteps.png
 */
Array multisteps(glm::ivec2       shape,
                 float            angle,
                 float            r = 1.2f,
                 int              nsteps = 8,
                 float            elevation_exponent = 0.7f,
                 float            shape_gain = 4.f,
                 float            scale = 0.5f,
                 float            outer_slope = 0.1f,
                 const Array     *p_ctrl_param = nullptr,
                 const Array     *p_noise_x = nullptr,
                 const Array     *p_noise_y = nullptr,
                 const glm::vec2 &center = {0.5f, 0.5f},
                 const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a paraboloid.
 * @param  shape                Array shape.
 * @param  angle                Rotation angle.
 * @param  a                    Curvature parameter, first principal axis.
 * @param  b                    Curvature parameter, second principal axis.
 * @param  v0                   Value at the paraboloid center.
 * @param  reverse_x            Reverse coefficient of first principal axis.
 * @param  reverse_y            Reverse coefficient of second principal axis.
 * @param  p_base_elevation     Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  center               Primitive reference center.
 * @param  bbox                 Domain bounding box.
 * @return                      Output array.
 *
 * **Example**
 * @include ex_paraboloid.cpp
 *
 * **Result**
 * @image html ex_paraboloid.png
 */
Array paraboloid(glm::ivec2   shape,
                 float        angle,
                 float        a,
                 float        b,
                 float        v0 = 0.f,
                 bool         reverse_x = false,
                 bool         reverse_y = false,
                 const Array *p_noise_x = nullptr,
                 const Array *p_noise_y = nullptr,
                 glm::vec2    center = {0.5f, 0.5f},
                 glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generate a polar-based shape mask with optional angular and radial
 * perturbations.
 *
 * Creates a 2D scalar field defined in polar space using radial bounds, sector
 * masking, aspect ratio deformation, and optional noise offsets on radius and
 * angle.
 *
 * @param  shape           Output array shape.
 * @param  rmin            Inner radius.
 * @param  rmax            Outer radius.
 * @param  aspect_ratio    Elliptical aspect ratio applied on Y axis.
 * @param  smoothing_width Smoothing width for radial and angular transitions.
 * @param  square_base     If true, remaps the shape toward a square profile.
 * @param  angle           Global rotation angle in degrees.
 * @param  sector_angle    Angular opening removed from the shape in degrees.
 * @param  vmin            Minimum angular modulation value.
 * @param  kt_value        Angular modulation frequency.
 * @param  kr_border       Angular modulation frequency for the outer border.
 * @param  kr_border_ratio Blend factor for outer border modulation.
 * @param  p_noise_r       Optional radial noise offset field.
 * @param  p_noise_theta   Optional angular noise offset field.
 * @param  center          Shape center.
 * @param  bbox            Bounding box used for grid generation.
 *
 * @return                 Generated scalar field.
 *
 * **Example**
 * @include ex_polar_shape.cpp
 *
 * **Result**
 * @image html ex_polar_shape.png
 */
Array polar_shape(glm::ivec2   shape,
                  float        rmin = 0.1f,
                  float        rmax = 0.3f,
                  float        aspect_ratio = 1.f,
                  float        smoothing_width = 0.1f,
                  bool         square_base = false,
                  float        angle = 15.f,
                  float        sector_angle = 90.f,
                  float        vmin = 0.5f,
                  float        kt_value = 0.f,
                  float        kr_border = 0.f,
                  float        kr_border_ratio = 0.1f,
                  const Array *p_noise_r = nullptr,
                  const Array *p_noise_theta = nullptr,
                  glm::vec2    center = {0.5f, 0.5f},
                  glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a surface defined by the four corner values of the domain.
 *
 * Interpolates values across the domain using a quadratic / bilinear surface
 * function derived from the four corner elevations.
 *
 * @param  shape        Dimensions of the output array (width, height).
 * @param  c00          Elevation at bottom-left corner (xmin, ymin).
 * @param  c10          Elevation at bottom-right corner (xmax, ymin).
 * @param  c01          Elevation at top-left corner (xmin, ymax).
 * @param  c11          Elevation at top-right corner (xmax, ymax).
 * @param  p_ctrl_param Optional pointer to a control parameter array
 * (multiplier).
 * @param  p_noise_x    Optional pointer to x-direction displacement noise
 * array.
 * @param  p_noise_y    Optional pointer to y-direction displacement noise
 * array.
 * @param  bbox         Bounding box defining the domain coordinates (xmin,
 * xmax, ymin, ymax).
 * @return              Array containing the generated quadratic surface.
 *
 * **Example**
 * @include ex_quad_surface.cpp
 *
 * **Result**
 * @image html ex_quad_surface.png
 */
Array quad_surface(glm::ivec2   shape,
                   float        c00 = 0.f,
                   float        c10 = 0.f,
                   float        c01 = 0.f,
                   float        c11 = 0.f,
                   const Array *p_ctrl_param = nullptr,
                   const Array *p_noise_x = nullptr,
                   const Array *p_noise_y = nullptr,
                   glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a rectangle-shaped heightmap with optional modifications.
 *
 * This function creates a 2D array representing a rectangle with specified
 * dimensions, rotation, and optional parameters for customization such as
 * control parameters, noise, and stretching.
 *
 * @param  shape        Dimensions of the output array (width, height).
 * @param  rx           Half-width of the rectangle, in normalized coordinates
 *                      (0.0 to 1.0).
 * @param  ry           Half-height of the rectangle, in normalized coordinates
 *                      (0.0 to 1.0).
 * @param  angle        Rotation angle of the rectangle in radians. Positive
 *                      values rotate counterclockwise.
 * @param  slope        Slope of the rectangle edge transition. A larger value
 *                      makes the edge transition sharper. Defaults to 1.0.
 * @param  p_ctrl_param Optional pointer to an `Array` controlling custom
 *                      parameters for the rectangle generation.
 * @param  p_noise_x    Optional pointer to an `Array` for adding noise in the
 *                      x-direction.
 * @param  p_noise_y    Optional pointer to an `Array` for adding noise in the
 *                      y-direction.
 * @param  p_stretching Optional pointer to an `Array` for stretching the
 *                      rectangle horizontally or vertically.
 * @param  center       Center of the rectangle in normalized coordinates (0.0
 *                      to 1.0). Defaults to {0.5, 0.5}.
 * @param  bbox         Bounding box for the rectangle in normalized coordinates
 *                     {x_min, x_max, y_min, y_max}. Defaults to {0.0, 1.0, 0.0,
 * 1.0}.
 *
 * @return              A 2D array representing the generated rectangle shape.
 *    *
 * **Example**
 * @include ex_rectangle.cpp
 *
 * **Result**
 * @image html ex_rectangle.png
 */
Array rectangle(glm::ivec2   shape,
                float        rx,
                float        ry,
                float        angle,
                float        slope = 1.f,
                const Array *p_ctrl_param = nullptr,
                const Array *p_noise_x = nullptr,
                const Array *p_noise_y = nullptr,
                glm::vec2    center = {0.5f, 0.5f},
                glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});
/**
 * @brief Return an array corresponding to a slope with a given overall.
 *
 * @param  shape                Array shape.
 * @param  angle                Overall rotation angle (in degree).
 * @param  slope                Slope (assuming a unit domain).
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the slope parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local coordinate multiplier.
 * @param  center               Primitive reference center.
 * @param  bbox                 Domain bounding box.
 * @return                      Array New array.
 *
 * **Example**
 * @include ex_slope.cpp
 *
 * **Result**
 * @image html ex_slope.png
 */
Array slope(glm::ivec2   shape,
            float        angle,
            float        slope,
            const Array *p_ctrl_param = nullptr,
            const Array *p_noise_x = nullptr,
            const Array *p_noise_y = nullptr,
            glm::vec2    center = {0.5f, 0.5f},
            glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});
/**
 * @brief Generates a smooth cosine array.
 */
Array smooth_cosine(glm::ivec2   shape,
                    const Array *p_noise_x,
                    const Array *p_noise_y,
                    glm::vec2    center = {0.5f, 0.5f},
                    glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a step function (Heaviside with an optional talus slope at the
 * transition).
 *
 * @param  shape                Array shape.
 * @param  angle                Overall rotation angle (in degree).
 * @param  slope                Step slope (assuming a unit domain).
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the slope parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local coordinate multiplier.
 * @param  center               Primitive reference center.
 * @param  bbox                 Domain bounding box.
 * @return                      Array New array.
 *
 * **Example**
 * @include ex_step.cpp
 *
 * **Result**
 * @image html ex_step.png
 */
Array step(glm::ivec2   shape,
           float        angle,
           float        slope,
           const Array *p_ctrl_param = nullptr,
           const Array *p_noise_x = nullptr,
           const Array *p_noise_y = nullptr,
           glm::vec2    center = {0.5f, 0.5f},
           glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generate displacements `dx` and `dy` to apply a swirl effect to
 * another primitve.
 *
 * @param dx[out]   'x' displacement (unit domain scale).
 * @param dy[out]   'y' displacement (unit domain scale).
 * @param amplitude Displacement amplitude.
 * @param exponent  Distance exponent.
 * @param p_noise   eference to the input noise array.
 * @param bbox      Domain bounding box.
 *
 * **Example**
 * @include ex_swirl.cpp
 *
 * **Result**
 * @image html ex_swirl.png
 */
void swirl(Array       &dx,
           Array       &dy,
           float        amplitude = 1.f,
           float        exponent = 1.f,
           const Array *p_noise = nullptr,
           glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a dune shape wave.
 *
 * @param  shape                Array shape.
 * @param  kw                   Wavenumber with respect to a unit domain.
 * @param  angle                Overall rotation angle (in degree).
 * @param  xtop                 Relative location of the top of the dune profile
 *                              (in [0, 1]).
 * @param  xbottom              Relative location of the foot of the dune
 *                              profile (in [0, 1]).
 * @param  phase_shift          Phase shift (in radians).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  center               Primitive reference center.
 * @param  bbox                 Domain bounding box.
 * @return                      Array New array.
 */
Array wave_dune(glm::ivec2   shape,
                float        kw,
                float        angle,
                float        xtop,
                float        xbottom,
                float        phase_shift = 0.f,
                const Array *p_noise_x = nullptr,
                const Array *p_noise_y = nullptr,
                glm::vec2    center = {0.5f, 0.5f},
                glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a sine wave.
 *
 * @param  shape                Array shape.
 * @param  kw                   Wavenumber with respect to a unit domain.
 * @param  angle                Overall rotation angle (in degree).
 * @param  phase_shift          Phase shift (in radians).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  center               Primitive reference center.
 * @param  bbox                 Domain bounding box.
 * @return                      Array New array.
 *
 * **Example**
 * @include ex_wave.cpp
 *
 * **Result**
 * @image html ex_wave0.png
 * @image html ex_wave1.png
 */

Array wave_sine(glm::ivec2   shape,
                float        kw,
                float        angle,
                float        phase_shift = 0.f,
                const Array *p_noise_x = nullptr,
                const Array *p_noise_y = nullptr,
                glm::vec2    center = {0.5f, 0.5f},
                glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a square wave.
 *
 * @param  shape                Array shape.
 * @param  kw                   Wavenumber with respect to a unit domain.
 * @param  angle                Overall rotation angle (in degree).
 * @param  phase_shift          Phase shift (in radians).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  center               Primitive reference center.
 * @param  bbox                 Domain bounding box.
 * @return                      Array New array.
 *
 * **Example**
 * @include ex_wave.cpp
 *
 * **Result**
 * @image html ex_wave0.png
 * @image html ex_wave1.png
 */
Array wave_square(glm::ivec2   shape,
                  float        kw,
                  float        angle,
                  float        phase_shift = 0.f,
                  const Array *p_noise_x = nullptr,
                  const Array *p_noise_y = nullptr,
                  glm::vec2    center = {0.5f, 0.5f},
                  glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return a triangular wave.
 *
 * @param  shape                Array shape.
 * @param  kw                   Wavenumber with respect to a unit domain.
 * @param  angle                Overall rotation angle (in degree).
 * @param  slant_ratio          Relative location of the triangle apex, in [0,
 *                              1].
 * @param  phase_shift          Phase shift (in radians).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  center               Primitive reference center.
 * @param  bbox                 Domain bounding box.
 * @return                      Array New array.
 *
 * **Example**
 * @include ex_wave.cpp
 *
 * **Result**
 * @image html ex_wave0.png
 * @image html ex_wave1.png
 */
Array wave_triangular(glm::ivec2   shape,
                      float        kw,
                      float        angle,
                      float        slant_ratio,
                      float        phase_shift = 0.f,
                      const Array *p_noise_x = nullptr,
                      const Array *p_noise_y = nullptr,
                      glm::vec2    center = {0.5f, 0.5f},
                      glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap

namespace hmap::gpu
{

/**
 * @brief Generates a scalar field representing the signed distance to randomly
 * generated hemispheres.
 *
 * @param  shape     Resolution of the output field (width, height).
 * @param  kw        Scaling factor for field coordinates (world units per
 *                   pixel).
 * @param  seed      Random seed generation.
 * @param  rmin      Minimum sphere radius (relative to bbox size).
 * @param  rmax      Maximum sphere radius (relative to bbox size).
 * @param  density   Fraction of pixels covered by polygons (approximate).
 * @param  jitter    Random displacement factor applied to sphere vertices.
 * @param  shift     Random position shift applied to sphere center.
 * @param  p_noise_x Optional displacement noise field in the X direction
 *                   (nullptr to disable).
 * @param  p_noise_y Optional displacement noise field in the Y direction
 *                   (nullptr to disable).
 * @param  bbox      Bounding box in normalized coordinates {xmin, xmax, ymin,
 *                   ymax}.
 * @return           Array         2D array containing the signed distance
 *                   field.
 *
 * **Example**
 * @include ex_phemisphere_field.cpp
 *
 * **Result**
 * @image html ex_phemisphere_field.png
 */
Array hemisphere_field(glm::ivec2    shape,
                       glm::vec2     kw,
                       std::uint32_t seed,
                       float         rmin = 0.05f,
                       float         rmax = 0.8f,
                       float         amplitude_random_ratio = 1.f,
                       float         density = 0.1f,
                       glm::vec2     jitter = {1.f, 1.f},
                       float         shift = 0.f,
                       const Array  *p_noise_x = nullptr,
                       const Array  *p_noise_y = nullptr,
                       const Array  *p_noise_distance = nullptr,
                       const Array  *p_density_multiplier = nullptr,
                       const Array  *p_size_multiplier = nullptr,
                       glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/*! @brief See hmap::hemisphere_field */
Array hemisphere_field_fbm(glm::ivec2    shape,
                           glm::vec2     kw,
                           std::uint32_t seed,
                           float         rmin = 0.05f,
                           float         rmax = 0.8f,
                           float         amplitude_random_ratio = 1.f,
                           float         density = 0.1f,
                           glm::vec2     jitter = {0.5f, 0.5f},
                           float         shift = 0.1f,
                           int           octaves = 8,
                           float         persistence = 0.5f,
                           float         lacunarity = 2.f,
                           const Array  *p_noise_x = nullptr,
                           const Array  *p_noise_y = nullptr,
                           const Array  *p_noise_distance = nullptr,
                           const Array  *p_density_multiplier = nullptr,
                           const Array  *p_size_multiplier = nullptr,
                           glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief GPU-accelerated multi-step height generation with procedural noise.
 *
 * Builds the same multi-step profile as the CPU version but injects a
 * Voronoi-FBM noise field to modulate the step positions. Noise can be inflated
 * or deflated and is projected onto the step axis before being passed to the
 * CPU multisteps() implementation.
 *
 * @param  shape              Output array resolution.
 * @param  angle              Rotation angle in degrees.
 * @param  seed               Noise seed.
 * @param  kw                 Base wave numbers for noise.
 * @param  noise_amp          Amplitude of injected noise.
 * @param  noise_rugosity     Roughness of the FBM noise.
 * @param  noise_inflate      Whether to invert noise direction.
 * @param  r                  Geometric ratio between steps.
 * @param  nsteps             Number of steps.
 * @param  elevation_exponent Exponent shaping step heights.
 * @param  shape_gain         Gain used to reshape intra-step transitions.
 * @param  scale              Scale of the step axis.
 * @param  outer_slope        Slope outside the [0,1] axis interval.
 * @param  p_ctrl_param       Optional control parameter array.
 * @param  center             Center of the step axis.
 * @param  bbox               Bounding box of the domain.
 * @return                    Generated Array.
 *
 * **Example**
 * @include ex_multisteps.cpp
 *
 * **Result**
 * @image html ex_multisteps.png
 */
Array multisteps(glm::ivec2       shape,
                 float            angle,
                 std::uint32_t    seed,
                 glm::vec2        kw = {2.f, 2.f},
                 float            noise_amp = 0.1f,
                 float            noise_rugosity = 0.f,
                 bool             noise_inflate = true,
                 float            r = 1.2f,
                 int              nsteps = 8,
                 float            elevation_exponent = 0.7f,
                 float            shape_gain = 4.f,
                 float            scale = 0.5f,
                 float            outer_slope = 0.1f,
                 const Array     *p_ctrl_param = nullptr,
                 const glm::vec2 &center = {0.5f, 0.5f},
                 const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a scalar field representing the signed distance to randomly
 * generated polygons.
 *
 * This function procedurally generates a field where each pixel stores the
 * signed distance to the nearest polygon, using a smooth signed distance
 * function (SDF) with optional clamping. The polygons are generated randomly
 * within the specified bounding box, with configurable vertex counts, sizes,
 * jitter, and density. Optionally, displacement noise fields can be applied to
 * polygon positions.
 *
 * @param  shape          Resolution of the output field (width, height).
 * @param  kw             Scaling factor for field coordinates (world units per
 *                        pixel).
 * @param  seed           Random seed for polygon generation.
 * @param  rmin           Minimum polygon radius (relative to bbox size).
 * @param  rmax           Maximum polygon radius (relative to bbox size).
 * @param  clamping_dist  Distance threshold for clamping the SDF value (used to
 *                        soften edges).
 * @param  clamping_k     Smoothness parameter for clamping.
 * @param  n_vertices_min Minimum number of vertices per polygon.
 * @param  n_vertices_max Maximum number of vertices per polygon.
 * @param  density        Fraction of pixels covered by polygons (approximate).
 * @param  jitter         Random displacement factor applied to polygon
 *                        vertices.
 * @param  shift          Random position shift applied to polygon center.
 * @param  p_noise_x      Optional displacement noise field in the X direction
 *                        (nullptr to disable).
 * @param  p_noise_y      Optional displacement noise field in the Y direction
 *                        (nullptr to disable).
 * @param  bbox           Bounding box in normalized coordinates {xmin, xmax,
 *                        ymin, ymax}.
 * @return                Array         2D array containing the signed distance
 *                        field.
 *
 * @note Polygons are randomly generated per call and are not guaranteed to be
 * convex.
 * @note The sign of the SDF is negative inside polygons and positive outside.
 *
 * **Example**
 * @include ex_polygon_field_fbm.cpp
 *
 * **Result**
 * @image html ex_polygon_field_fbm.png
 */
Array polygon_field(glm::ivec2    shape,
                    glm::vec2     kw,
                    std::uint32_t seed,
                    float         rmin = 0.05f,
                    float         rmax = 0.8f,
                    float         clamping_dist = 0.1f,
                    float         clamping_k = 0.1f,
                    int           n_vertices_min = 3,
                    int           n_vertices_max = 16,
                    float         density = 0.5f,
                    glm::vec2     jitter = {0.5f, 0.5f},
                    float         shift = 0.1f,
                    const Array  *p_noise_x = nullptr,
                    const Array  *p_noise_y = nullptr,
                    const Array  *p_noise_distance = nullptr,
                    const Array  *p_density_multiplier = nullptr,
                    const Array  *p_size_multiplier = nullptr,
                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a scalar field representing the signed distance to randomly
 * generated polygons combined with fractal Brownian motion (fBm) noise
 * modulation.
 *
 * Similar to polygon_field(), but the generated SDF is modulated by an fBm
 * noise function to create more natural, irregular shapes. The fBm parameters
 * allow control over the noise persistence, lacunarity, and number of octaves.
 *
 * @param  shape          Resolution of the output field (width, height).
 * @param  kw             Scaling factor for field coordinates (world units per
 *                        pixel).
 * @param  seed           Random seed for polygon generation and fBm noise.
 * @param  rmin           Minimum polygon radius (relative to bbox size).
 * @param  rmax           Maximum polygon radius (relative to bbox size).
 * @param  clamping_dist  Distance threshold for clamping the SDF value (used to
 *                        soften edges).
 * @param  clamping_k     Smoothness parameter for clamping.
 * @param  n_vertices_min Minimum number of vertices per polygon.
 * @param  n_vertices_max Maximum number of vertices per polygon.
 * @param  density        Fraction of pixels covered by polygons (approximate).
 * @param  jitter         Random displacement factor applied to polygon
 *                        vertices.
 * @param  shift          Random position shift applied to polygon center.
 * @param  octaves        Number of fBm octaves.
 * @param  persistence    Amplitude decay per octave in fBm.
 * @param  lacunarity     Frequency multiplier per octave in fBm.
 * @param  p_noise_x      Optional displacement noise field in the X direction
 *                        (nullptr to disable).
 * @param  p_noise_y      Optional displacement noise field in the Y direction
 *                        (nullptr to disable).
 * @param  bbox           Bounding box in normalized coordinates {xmin, xmax,
 *                        ymin, ymax}.
 * @return                Array         2D array containing the signed distance
 *                        field modulated by fBm.
 *
 * @note The sign of the SDF is negative inside polygons and positive outside.
 *
 * **Example**
 * @include ex_polygon_field_fbm.cpp
 *
 * **Result**
 * @image html ex_polygon_field_fbm.png
 */
Array polygon_field_fbm(glm::ivec2    shape,
                        glm::vec2     kw,
                        std::uint32_t seed,
                        float         rmin = 0.05f,
                        float         rmax = 0.8f,
                        float         clamping_dist = 0.1f,
                        float         clamping_k = 0.1f,
                        int           n_vertices_min = 3,
                        int           n_vertices_max = 16,
                        float         density = 0.1f,
                        glm::vec2     jitter = {0.5f, 0.5f},
                        float         shift = 0.1f,
                        int           octaves = 8,
                        float         persistence = 0.5f,
                        float         lacunarity = 2.f,
                        const Array  *p_noise_x = nullptr,
                        const Array  *p_noise_y = nullptr,
                        const Array  *p_noise_distance = nullptr,
                        const Array  *p_density_multiplier = nullptr,
                        const Array  *p_size_multiplier = nullptr,
                        glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap::gpu