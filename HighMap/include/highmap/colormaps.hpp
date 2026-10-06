/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file colormaps.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Predefined color mapping schemes for data visualization, including
 * perceptual uniform (Viridis, Inferno), scientific (Seismic, Terrain), and
 * gradient-based (Jet, Hot) colormaps. Provides enum mappings and RGB data for
 * seamless integration with visualization pipelines.
 *
 * @copyright Copyright (c) 2023
 */
#pragma once
#include <map>
#include <string>
#include <vector>

#include <glm/vec3.hpp>

namespace hmap
{

/**
 * @enum Cmap
 * @brief Enumeration of colormap types for image processing.
 *
 * This enumeration defines various colormap options that can be used for
 * visualizing data in image processing. Each colormap provides a different
 * color mapping scheme, which can be applied to grayscale images or data to
 * enhance visual interpretation.
 *
 * @var Cmap::BONE A colormap that resembles human bone, providing a
 * grayscale-to-blue gradient.
 *
 * @var Cmap::GRAY A simple grayscale colormap, mapping data values directly to
 * shades of gray.
 *
 * @var Cmap::HOT A colormap that transitions from black through red and yellow
 * to white, often used to represent heatmaps.
 *
 * @var Cmap::INFERNO A perceptually uniform colormap that transitions from dark
 * blue to yellow.
 *
 * @var Cmap::JET A widely used colormap that transitions from blue to red via
 * cyan and yellow, often criticized for creating misleading visual
 * interpretations.
 *
 * @var Cmap::MAGMA A perceptually uniform colormap that transitions from dark
 * purple to yellow.
 *
 * @var Cmap::NIPY_SPECTRAL A colormap that spans a wide range of colors, useful
 * for scientific visualization.
 *
 * @var Cmap::SEISMIC A diverging colormap that transitions from blue to red,
 * typically used for visualizing differences from a central value.
 *
 * @var Cmap::TERRAIN A colormap resembling natural terrain, with shades of
 * green, brown, and white, often used in topographical data.
 *
 * @var Cmap::VIRIDIS A perceptually uniform colormap that transitions from dark
 * blue to yellow-green.
 *
 * @var Cmap::WHITE_UNIFORM A uniform white colormap.
 *
 * **Example**
 * @include ex_colormaps.cpp
 *
 * **Result**
 * @image html ex_colormaps1.png
 * @image html ex_colormaps2.png
 * @image html ex_colormaps3.png
 * @image html ex_colormaps4.png
 * @image html ex_colormaps5.png
 * @image html ex_colormaps6.png
 * @image html ex_colormaps7.png
 * @image html ex_colormaps8.png
 * @image html ex_colormaps9.png
 * @image html ex_colormaps10.png
 */
enum Cmap : int
{
  BONE,
  GRAY,
  HOT,
  INFERNO,
  JET,
  MAGMA,
  NIPY_SPECTRAL,
  SEISMIC,
  TERRAIN,
  TURBO,
  VIRIDIS,
  WHITE_UNIFORM
};

/**
 * @brief Mapping between the `Cmap` enum and its corresponding lowercase string
 * representation.
 *
 * This `std::map` provides a lookup table that associates each colormap defined
 * in the `Cmap` enum with its plain text name in lowercase. This mapping allows
 * for easy conversion between the enum values and their lowercase string
 * equivalents, which can be useful for debugging, logging, or displaying the
 * colormap names in a user interface.
 *
 * Example usage:
 * @code Cmap selected_cmap = Cmap::JET;
 * std::cout << "Selected colormap: " << cmap_as_string[selected_cmap] <<
 * std::endl;
 * @endcode
 */
static std::map<Cmap, std::string> cmap_as_string = {
    {BONE, "bone"},
    {GRAY, "gray"},
    {HOT, "hot"},
    {INFERNO, "inferno"},
    {JET, "jet"},
    {MAGMA, "magma"},
    {NIPY_SPECTRAL, "nipy_spectral"},
    {SEISMIC, "seismic"},
    {TERRAIN, "terrain"},
    {TURBO, "turbo"},
    {VIRIDIS, "viridis"},
    {WHITE_UNIFORM, "white_uniform"}};

#define CMAP_BONE                                                              \
  {                                                                            \
      {0.000f, 0.000f, 0.000f},                                                \
      {0.058f, 0.058f, 0.081f},                                                \
      {0.117f, 0.117f, 0.162f},                                                \
      {0.175f, 0.175f, 0.243f},                                                \
      {0.233f, 0.233f, 0.325f},                                                \
      {0.292f, 0.292f, 0.406f},                                                \
      {0.350f, 0.361f, 0.475f},                                                \
      {0.408f, 0.442f, 0.533f},                                                \
      {0.467f, 0.522f, 0.592f},                                                \
      {0.525f, 0.602f, 0.650f},                                                \
      {0.583f, 0.682f, 0.708f},                                                \
      {0.642f, 0.762f, 0.767f},                                                \
      {0.727f, 0.825f, 0.825f},                                                \
      {0.818f, 0.883f, 0.883f},                                                \
      {0.909f, 0.942f, 0.942f},                                                \
      {1.000f, 1.000f, 1.000f},                                                \
  }

#define CMAP_GRAY                                                              \
  {                                                                            \
      {0.000f, 0.000f, 0.000f},                                                \
      {0.067f, 0.067f, 0.067f},                                                \
      {0.133f, 0.133f, 0.133f},                                                \
      {0.200f, 0.200f, 0.200f},                                                \
      {0.267f, 0.267f, 0.267f},                                                \
      {0.333f, 0.333f, 0.333f},                                                \
      {0.400f, 0.400f, 0.400f},                                                \
      {0.467f, 0.467f, 0.467f},                                                \
      {0.533f, 0.533f, 0.533f},                                                \
      {0.600f, 0.600f, 0.600f},                                                \
      {0.667f, 0.667f, 0.667f},                                                \
      {0.733f, 0.733f, 0.733f},                                                \
      {0.800f, 0.800f, 0.800f},                                                \
      {0.867f, 0.867f, 0.867f},                                                \
      {0.933f, 0.933f, 0.933f},                                                \
      {1.000f, 1.000f, 1.000f},                                                \
  }

#define CMAP_HOT                                                               \
  {                                                                            \
      {0.042f, 0.000f, 0.000f},                                                \
      {0.217f, 0.000f, 0.000f},                                                \
      {0.392f, 0.000f, 0.000f},                                                \
      {0.567f, 0.000f, 0.000f},                                                \
      {0.742f, 0.000f, 0.000f},                                                \
      {0.917f, 0.000f, 0.000f},                                                \
      {1.000f, 0.092f, 0.000f},                                                \
      {1.000f, 0.267f, 0.000f},                                                \
      {1.000f, 0.442f, 0.000f},                                                \
      {1.000f, 0.617f, 0.000f},                                                \
      {1.000f, 0.792f, 0.000f},                                                \
      {1.000f, 0.967f, 0.000f},                                                \
      {1.000f, 1.000f, 0.212f},                                                \
      {1.000f, 1.000f, 0.475f},                                                \
      {1.000f, 1.000f, 0.737f},                                                \
      {1.000f, 1.000f, 1.000f},                                                \
  }

#define CMAP_INFERNO                                                           \
  {                                                                            \
      {0.001f, 0.000f, 0.014f},                                                \
      {0.047f, 0.030f, 0.150f},                                                \
      {0.142f, 0.046f, 0.309f},                                                \
      {0.258f, 0.039f, 0.406f},                                                \
      {0.367f, 0.072f, 0.432f},                                                \
      {0.472f, 0.111f, 0.428f},                                                \
      {0.578f, 0.148f, 0.404f},                                                \
      {0.683f, 0.190f, 0.361f},                                                \
      {0.781f, 0.243f, 0.300f},                                                \
      {0.865f, 0.317f, 0.226f},                                                \
      {0.930f, 0.411f, 0.145f},                                                \
      {0.971f, 0.523f, 0.058f},                                                \
      {0.988f, 0.645f, 0.040f},                                                \
      {0.979f, 0.775f, 0.176f},                                                \
      {0.950f, 0.903f, 0.380f},                                                \
      {0.988f, 0.998f, 0.645f},                                                \
  }

#define CMAP_JET                                                               \
  {                                                                            \
      {0.000f, 0.000f, 0.500f},                                                \
      {0.000f, 0.000f, 0.803f},                                                \
      {0.000f, 0.033f, 1.000f},                                                \
      {0.000f, 0.300f, 1.000f},                                                \
      {0.000f, 0.567f, 1.000f},                                                \
      {0.000f, 0.833f, 1.000f},                                                \
      {0.161f, 1.000f, 0.806f},                                                \
      {0.376f, 1.000f, 0.591f},                                                \
      {0.591f, 1.000f, 0.376f},                                                \
      {0.806f, 1.000f, 0.161f},                                                \
      {1.000f, 0.901f, 0.000f},                                                \
      {1.000f, 0.654f, 0.000f},                                                \
      {1.000f, 0.407f, 0.000f},                                                \
      {1.000f, 0.160f, 0.000f},                                                \
      {0.803f, 0.000f, 0.000f},                                                \
      {0.500f, 0.000f, 0.000f},                                                \
  }

#define CMAP_MAGMA                                                             \
  {                                                                            \
      {0.001f, 0.000f, 0.014f},                                                \
      {0.044f, 0.034f, 0.142f},                                                \
      {0.124f, 0.067f, 0.296f},                                                \
      {0.232f, 0.060f, 0.438f},                                                \
      {0.341f, 0.081f, 0.493f},                                                \
      {0.445f, 0.123f, 0.507f},                                                \
      {0.550f, 0.161f, 0.506f},                                                \
      {0.658f, 0.196f, 0.490f},                                                \
      {0.767f, 0.234f, 0.458f},                                                \
      {0.869f, 0.288f, 0.409f},                                                \
      {0.944f, 0.378f, 0.365f},                                                \
      {0.981f, 0.498f, 0.370f},                                                \
      {0.995f, 0.624f, 0.427f},                                                \
      {0.997f, 0.748f, 0.517f},                                                \
      {0.993f, 0.870f, 0.626f},                                                \
      {0.987f, 0.991f, 0.750f},                                                \
  }

#define CMAP_NIPY_SPECTRAL                                                     \
  {                                                                            \
      {0.000f, 0.000f, 0.000f},                                                \
      {0.489f, 0.000f, 0.556f},                                                \
      {0.178f, 0.000f, 0.644f},                                                \
      {0.000f, 0.000f, 0.867f},                                                \
      {0.000f, 0.511f, 0.867f},                                                \
      {0.000f, 0.644f, 0.733f},                                                \
      {0.000f, 0.667f, 0.533f},                                                \
      {0.000f, 0.644f, 0.000f},                                                \
      {0.000f, 0.822f, 0.000f},                                                \
      {0.000f, 1.000f, 0.000f},                                                \
      {0.800f, 0.978f, 0.000f},                                                \
      {0.978f, 0.844f, 0.000f},                                                \
      {1.000f, 0.600f, 0.000f},                                                \
      {0.956f, 0.000f, 0.000f},                                                \
      {0.822f, 0.000f, 0.000f},                                                \
      {0.800f, 0.800f, 0.800f},                                                \
  }

#define CMAP_SEISMIC                                                           \
  {                                                                            \
      {0.000f, 0.000f, 0.300f},                                                \
      {0.000f, 0.000f, 0.487f},                                                \
      {0.000f, 0.000f, 0.673f},                                                \
      {0.000f, 0.000f, 0.860f},                                                \
      {0.067f, 0.067f, 1.000f},                                                \
      {0.333f, 0.333f, 1.000f},                                                \
      {0.600f, 0.600f, 1.000f},                                                \
      {0.867f, 0.867f, 1.000f},                                                \
      {1.000f, 0.867f, 0.867f},                                                \
      {1.000f, 0.600f, 0.600f},                                                \
      {1.000f, 0.333f, 0.333f},                                                \
      {1.000f, 0.067f, 0.067f},                                                \
      {0.900f, 0.000f, 0.000f},                                                \
      {0.767f, 0.000f, 0.000f},                                                \
      {0.633f, 0.000f, 0.000f},                                                \
      {0.500f, 0.000f, 0.000f},                                                \
  }

#define CMAP_TERRAIN                                                           \
  {                                                                            \
      {0.200f, 0.200f, 0.600f},                                                \
      {0.111f, 0.378f, 0.778f},                                                \
      {0.022f, 0.556f, 0.956f},                                                \
      {0.000f, 0.700f, 0.700f},                                                \
      {0.067f, 0.813f, 0.413f},                                                \
      {0.333f, 0.867f, 0.467f},                                                \
      {0.600f, 0.920f, 0.520f},                                                \
      {0.867f, 0.973f, 0.573f},                                                \
      {0.933f, 0.915f, 0.564f},                                                \
      {0.800f, 0.744f, 0.492f},                                                \
      {0.667f, 0.573f, 0.420f},                                                \
      {0.533f, 0.403f, 0.348f},                                                \
      {0.600f, 0.488f, 0.464f},                                                \
      {0.733f, 0.659f, 0.643f},                                                \
      {0.867f, 0.829f, 0.821f},                                                \
      {1.000f, 1.000f, 1.000f},                                                \
  }

#define CMAP_TURBO                                                             \
  {                                                                            \
      {0.190f, 0.072f, 0.232f},                                                \
      {0.267f, 0.328f, 0.764f},                                                \
      {0.266f, 0.560f, 0.996f},                                                \
      {0.127f, 0.780f, 0.876f},                                                \
      {0.163f, 0.936f, 0.631f},                                                \
      {0.489f, 0.998f, 0.338f},                                                \
      {0.756f, 0.951f, 0.203f},                                                \
      {0.945f, 0.796f, 0.228f},                                                \
      {0.994f, 0.564f, 0.159f},                                                \
      {0.916f, 0.305f, 0.051f},                                                \
      {0.746f, 0.131f, 0.009f},                                                \
      {0.480f, 0.016f, 0.011f},                                                \
  }

#define CMAP_VIRIDIS                                                           \
  {                                                                            \
      {0.267f, 0.005f, 0.329f},                                                \
      {0.283f, 0.100f, 0.422f},                                                \
      {0.277f, 0.185f, 0.490f},                                                \
      {0.254f, 0.265f, 0.530f},                                                \
      {0.222f, 0.339f, 0.549f},                                                \
      {0.191f, 0.407f, 0.556f},                                                \
      {0.164f, 0.471f, 0.558f},                                                \
      {0.139f, 0.534f, 0.555f},                                                \
      {0.121f, 0.596f, 0.544f},                                                \
      {0.135f, 0.659f, 0.518f},                                                \
      {0.208f, 0.719f, 0.473f},                                                \
      {0.328f, 0.774f, 0.407f},                                                \
      {0.478f, 0.821f, 0.318f},                                                \
      {0.647f, 0.858f, 0.210f},                                                \
      {0.825f, 0.885f, 0.106f},                                                \
      {0.993f, 0.906f, 0.144f},                                                \
  }

#define CMAP_WHITE_UNIFORM                                                     \
  {                                                                            \
      {1.000f, 1.000f, 1.000f},                                                \
      {1.000f, 1.000f, 1.000f},                                                \
  }

std::vector<glm::vec3> get_colormap_data(int cmap);

} // namespace hmap
