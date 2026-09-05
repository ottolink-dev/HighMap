/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <limits>
#include <queue>
#include <random>
#include <utility>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/transform.hpp"

namespace hmap
{

// find the cut path from top to bottom with the minimum cost using dynamic
// programming
void find_vertical_cut_path(const Array &error, std::vector<int> &path_i)
{
  if (!validate_non_empty(error)) return;

  glm::ivec2 shape = error.shape;

  path_i.resize(shape.y);

  if (shape.x == 1)
  {
    std::fill(path_i.begin(), path_i.end(), 0);
    return;
  }

  // dynamic programming table for cumulative minimum error
  std::vector<float> dp_prev(shape.x);
  std::vector<float> dp_curr(shape.x);
  // store chosen offsets: -1, 0, or +1
  std::vector<int8_t> backtrack(shape.x * shape.y, 0);

  for (int i = 0; i < shape.x; i++)
    dp_prev[i] = error(i, 0);

  for (int j = 1; j < shape.y; j++)
  {
    for (int i = 0; i < shape.x; i++)
    {
      float  min_cost = dp_prev[i];
      int8_t best_offset = 0;

      if (i > 0 && dp_prev[i - 1] < min_cost)
      {
        min_cost = dp_prev[i - 1];
        best_offset = -1;
      }
      if (i + 1 < shape.x && dp_prev[i + 1] < min_cost)
      {
        min_cost = dp_prev[i + 1];
        best_offset = 1;
      }

      dp_curr[i] = error(i, j) + min_cost;
      backtrack[j * shape.x + i] = best_offset;
    }
    dp_prev = dp_curr;
  }

  // find index with minimum cumulative error in the bottom row
  int   best_i = 0;
  float min_total_cost = dp_prev[0];
  for (int i = 1; i < shape.x; i++)
  {
    if (dp_prev[i] < min_total_cost)
    {
      min_total_cost = dp_prev[i];
      best_i = i;
    }
  }

  // backtrack from bottom to top
  path_i[shape.y - 1] = best_i;
  int curr_i = best_i;
  for (int j = shape.y - 1; j > 0; j--)
  {
    curr_i += backtrack[j * shape.x + curr_i];
    path_i[j - 1] = curr_i;
  }
}

Array generate_mask(glm::ivec2 shape, std::vector<int> cut_path_i, int ir)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_non_empty(cut_path_i, "cut_path_i")) return Array(shape);
  if ((int)cut_path_i.size() < shape.y)
  {
    hmap::log::warn(std::source_location::current(),
                    "cut_path_i size ({}) is smaller than shape height ({})",
                    cut_path_i.size(),
                    shape.y);
    return Array(shape);
  }

  Array mask = Array(shape);

  // make sure there will be a minimum transition length in the mask
  // in order to avoid visible discontinuities
  int istart = std::min(ir, (int)(0.5f * shape.x));

  for (int j = 0; j < shape.y; j++)
    for (int i = istart; i < shape.x; i++)
      mask(i, j) = (i > cut_path_i[j]) ? 1.f : 0.f;

  for (int j = 0; j < shape.y; j++)
    mask(shape.x - 1, j) = 1.f;

  // smoothing
  smooth_cpulse(mask, ir);

  // eventually make sure the mask amplitude, used for lerping, is
  // in [0, 1] on every horizontal line of the array
  for (int j = 0; j < shape.y; j++)
  {
    float vmin = std::numeric_limits<float>::max();
    float vmax = -std::numeric_limits<float>::max();
    for (int i = 0; i < shape.x; i++)
    {
      vmin = std::min(vmin, mask(i, j));
      vmax = std::max(vmax, mask(i, j));
    }

    if (vmin != vmax)
      for (int i = 0; i < shape.x; i++)
      {
        mask(i, j) -= vmin;
        mask(i, j) /= (vmax - vmin);
      }
  }

  return mask;
}

void helper_flip_rot_transpose(Array &array,
                               bool   do_flip,
                               int    rot_steps,
                               bool   do_transpose)
{
  if (do_flip) flip_lr(array);

  // square patches support arbitrary 90-degree step rotations and transposition
  if (array.shape.x == array.shape.y)
  {
    if (rot_steps == 1)
      rot90(array);
    else if (rot_steps == 2)
      rot180(array);
    else if (rot_steps == 3)
      rot270(array);

    if (do_transpose) transpose(array);
  }
}

Array get_random_patch(const Array          &array,
                       glm::ivec2            patch_shape,
                       std::mt19937         &gen,
                       bool                  patch_flip,
                       bool                  patch_rotate,
                       bool                  patch_transpose,
                       std::vector<Array *> *p_secondary_arrays,
                       std::vector<Array>   *p_secondary_patches)
{
  if (!validate_non_empty(array) || !validate_shape(patch_shape))
    return Array();

  // clamp patch shape to array bounds
  int pw = std::min(patch_shape.x, array.shape.x);
  int ph = std::min(patch_shape.y, array.shape.y);

  // uniform sampling of valid start coordinates across the source heightmap
  int max_i = std::max(0, array.shape.x - pw);
  int max_j = std::max(0, array.shape.y - ph);

  std::uniform_int_distribution<int> dis_i(0, max_i);
  std::uniform_int_distribution<int> dis_j(0, max_j);
  std::uniform_int_distribution<int> dis_rot(0, 3);
  std::bernoulli_distribution        dis_bool(0.5);

  int i_start = dis_i(gen);
  int j_start = dis_j(gen);
  int i_end = i_start + pw;
  int j_end = j_start + ph;

  Array patch = array.extract_slice(glm::ivec4(i_start, i_end, j_start, j_end));

  // sample transformations uniformly
  bool do_flip = patch_flip && dis_bool(gen);
  int  rot_steps = patch_rotate ? dis_rot(gen) : 0;
  bool do_transpose = patch_transpose && dis_bool(gen);

  helper_flip_rot_transpose(patch, do_flip, rot_steps, do_transpose);

  // apply the patch extraction with the same parameters to the
  // secondary arrays
  if (p_secondary_arrays && p_secondary_patches)
  {
    p_secondary_patches->clear();

    for (auto pa : *p_secondary_arrays)
    {
      if (!pa || !validate_same_shape(array, *pa)) continue;

      Array sec_patch = pa->extract_slice(
          glm::ivec4(i_start, i_end, j_start, j_end));
      helper_flip_rot_transpose(sec_patch, do_flip, rot_steps, do_transpose);
      p_secondary_patches->push_back(sec_patch);
    }
  }

  return patch;
}

} // namespace hmap
