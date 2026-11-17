/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file roads.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for road network generation functions.
 *
 * @copyright Copyright (c) 2023
 */
#pragma once

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Generate a large scale road network between "nodes" (or cities) using
 * the alpha model.
 *
 * Based on the paper of Molinero et al. @cite Molinero2020.
 *
 * @param  xc            `x` locations of the cities.
 * @param  yc            `y` locations of the cities.
 * @param  size          City sizes (arbitrary unit).
 * @param  bbox          Bounding box.
 * @param  z             Heightmap.
 * @param  seed          Random seed number.
 * @param  alpha         "Alpha coefficient" of the model (cost difference
 *                       between creating a new and using an old one).
 * @param  n_dummy_nodes Number of dummy nodes added for the tesselation.
 * @param  dz_weight     Weight of the elevation difference in the cost
 *                       function.
 * @param  p_weight      Map defining the local value added to the cost function
 *                       (can be used to add large penalizations to create "no
 *                       go" zones).
 * @return               Graph Road network, edge weights correspond to the road
 *                       "size".
 *
 * **Example**
 * @include ex_generate_network_alpha_model.cpp
 *
 * **Result**
 * @image html ex_generate_network_alpha_model0.png
 * @image html ex_generate_network_alpha_model1.png
 * @image html ex_generate_network_alpha_model2.png
 */
Graph generate_network_alpha_model(const std::vector<float> &xc,
                                   const std::vector<float> &yc,
                                   const std::vector<float> &size,
                                   glm::vec4                 bbox,
                                   const Array              &z,
                                   std::uint32_t             seed,
                                   float                     alpha = 0.7f,
                                   int          n_dummy_nodes = 2500,
                                   float        dz_weight = 1.f,
                                   const Array *p_weight = nullptr);

} // namespace hmap
