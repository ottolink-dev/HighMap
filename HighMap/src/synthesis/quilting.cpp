/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/operator.hpp"
#include "highmap/transform.hpp"

namespace hmap
{

Array quilting(const std::vector<const Array *> &p_arrays,
               glm::ivec2                        patch_base_shape,
               glm::ivec2                        tiling,
               float                             overlap,
               std::uint32_t                     seed,
               std::vector<Array *>              secondary_arrays,
               bool                              patch_flip,
               bool                              patch_rotate,
               bool                              patch_transpose,
               float                             filter_width_ratio)
{
  if (p_arrays.empty() || !validate_shape(patch_base_shape) ||
      !validate_shape(tiling))
    return Array();

  for (const auto *p : p_arrays)
    if (!p || !validate_non_empty(*p)) return Array();

  for (const auto *p : secondary_arrays)
    if (!p || !validate_non_empty(*p)) return Array();

  std::mt19937 gen(seed);

  std::uniform_int_distribution<int> dis_a(0, (int)p_arrays.size() - 1);

  overlap = std::max(0.f, overlap);

  glm::ivec2 patch_shape = {(int)(patch_base_shape.x * (1.f + overlap)),
                            (int)(patch_base_shape.y * (1.f + overlap))};

  // keep size under control across all input arrays
  for (const auto *p : p_arrays)
  {
    patch_shape.x = std::min(patch_shape.x, p->shape.x);
    patch_shape.y = std::min(patch_shape.y, p->shape.y);
  }

  // ensure patch is at least patch_base_shape
  patch_shape.x = std::max(patch_shape.x, patch_base_shape.x);
  patch_shape.y = std::max(patch_shape.y, patch_base_shape.y);

  glm::ivec2 noverlap = {patch_shape.x - patch_base_shape.x,
                         patch_shape.y - patch_base_shape.y};

  glm::ivec2 shape_output = {tiling.x * patch_base_shape.x + noverlap.x,
                             tiling.y * patch_base_shape.y + noverlap.y};

  Array array_out = Array(shape_output);

  // reserve outputs for secondary arrays transformed using primary
  // array(s) as guide(s)
  std::vector<Array> secondary_arrays_output(secondary_arrays.size(),
                                             shape_output);

  // smoothing radius for the patch transitions
  int ir = (int)(noverlap.x * filter_width_ratio);

  for (int jt = 0; jt < tiling.y; jt++)
  {
    int   j1 = jt * patch_base_shape.y; // tile start
    Array array_strip = Array(glm::ivec2(array_out.shape.x, patch_shape.y));

    std::vector<Array> secondary_arrays_strips(secondary_arrays_output.size(),
                                               array_strip);

    // build up domain-wide horizontal strips
    for (int it = 0; it < tiling.x; it++)
    {
      int i1 = it * patch_base_shape.x;

      int array_idx = dis_a(gen);

      // used only if secondary_arrays is not empty
      std::vector<Array> secondary_patches = {};

      Array patch = get_random_patch(*p_arrays[array_idx],
                                     patch_shape,
                                     gen,
                                     patch_flip,
                                     patch_rotate,
                                     patch_transpose,
                                     &secondary_arrays,
                                     &secondary_patches);

      if (it > 0 && noverlap.x > 0)
      {
        Array error = Array(glm::ivec2(noverlap.x, patch_shape.y));
        for (int j = 0; j < patch_shape.y; j++)
          for (int i = 0; i < noverlap.x; i++)
            error(i, j) = std::abs(array_strip(i1 + i, j) - patch(i, j));

        // find cut path
        std::vector<int> cut_path_i;
        find_vertical_cut_path(error, cut_path_i);

        // define lerp factor
        Array mask = generate_mask(error.shape, cut_path_i, ir);

        // primary patches
        for (int j = 0; j < patch_shape.y; j++)
          for (int i = 0; i < noverlap.x; i++)
            patch(i, j) = lerp(array_strip(i1 + i, j), patch(i, j), mask(i, j));

        // same for secondary patches
        for (size_t k = 0; k < secondary_patches.size(); k++)
        {
          for (int j = 0; j < patch_shape.y; j++)
            for (int i = 0; i < noverlap.x; i++)
              secondary_patches[k](i, j) = lerp(
                  secondary_arrays_strips[k](i1 + i, j),
                  secondary_patches[k](i, j),
                  mask(i, j));
        }
      }

      // primary array
      for (int j = 0; j < patch_shape.y; j++)
        for (int i = 0; i < patch_shape.x; i++)
          array_strip(i1 + i, j) = patch(i, j);

      // apply the same patching to the secondary arrays
      for (size_t k = 0; k < secondary_arrays_output.size(); k++)
      {
        for (int j = 0; j < patch_shape.y; j++)
          for (int i = 0; i < patch_shape.x; i++)
            secondary_arrays_strips[k](i1 + i, j) = secondary_patches[k](i, j);
      }
    }

    // patch the horizontal stripes
    if (jt > 0 && noverlap.y > 0)
    {
      Array error = Array(glm::ivec2(shape_output.x, noverlap.y));

      for (int j = 0; j < noverlap.y; j++)
        for (int i = 0; i < shape_output.x; i++)
          error(i, j) = std::abs(array_strip(i, j) - array_out(i, j1 + j));

      Array mask = Array(error.shape);
      {
        Array            error_t = transpose(error);
        std::vector<int> cut_path_i;
        find_vertical_cut_path(error_t, cut_path_i);
        Array mask_t = generate_mask(error_t.shape,
                                     cut_path_i,
                                     (int)(noverlap.y * filter_width_ratio));
        mask = transpose(mask_t);
      }

      // primary array
      for (int i = 0; i < shape_output.x; i++)
      {
        for (int j = 0; j < noverlap.y; j++)
          array_out(i, j1 + j) = lerp(array_out(i, j1 + j),
                                      array_strip(i, j),
                                      mask(i, j));

        for (int j = noverlap.y; j < patch_shape.y; j++)
          array_out(i, j1 + j) = array_strip(i, j);
      }

      // apply the same patching to the secondary arrays
      for (size_t k = 0; k < secondary_arrays_output.size(); k++)
      {
        for (int i = 0; i < shape_output.x; i++)
        {
          for (int j = 0; j < noverlap.y; j++)
            secondary_arrays_output[k](i, j1 + j) = lerp(
                secondary_arrays_output[k](i, j1 + j),
                secondary_arrays_strips[k](i, j),
                mask(i, j));

          for (int j = noverlap.y; j < patch_shape.y; j++)
            secondary_arrays_output[k](i,
                                       j1 + j) = secondary_arrays_strips[k](i,
                                                                            j);
        }
      }
    }
    else
    {
      // primary array
      for (int j = 0; j < patch_shape.y; j++)
        for (int i = 0; i < array_out.shape.x; i++)
          array_out(i, j1 + j) = array_strip(i, j);

      // apply the same patching to the secondary arrays
      for (size_t k = 0; k < secondary_arrays_output.size(); k++)
      {
        for (int j = 0; j < patch_shape.y; j++)
          for (int i = 0; i < array_out.shape.x; i++)
            secondary_arrays_output[k](i,
                                       j1 + j) = secondary_arrays_strips[k](i,
                                                                            j);
      }
    }
  }

  // override content of input secondary arrays with output
  for (size_t k = 0; k < secondary_arrays_output.size(); k++)
    *secondary_arrays[k] = secondary_arrays_output[k];

  return array_out;
}

Array quilting_blend(const std::vector<const Array *> &p_arrays,
                     glm::ivec2                        patch_base_shape,
                     float                             overlap,
                     std::uint32_t                     seed,
                     bool                              patch_flip,
                     bool                              patch_rotate,
                     bool                              patch_transpose,
                     float                             filter_width_ratio)
{
  if (p_arrays.empty() || !validate_shape(patch_base_shape)) return Array();

  for (const auto *p : p_arrays)
    if (!p || !validate_non_empty(*p)) return Array();

  glm::ivec2 shape = p_arrays.back()->shape;

  glm::ivec2 tiling = glm::ivec2(
      (int)(std::ceil(shape.x / patch_base_shape.x)),
      (int)(std::ceil(shape.y / patch_base_shape.y)));

  Array array_out = quilting(p_arrays,
                             patch_base_shape,
                             tiling,
                             overlap,
                             seed,
                             {}, // no secondary arrays
                             patch_flip,
                             patch_rotate,
                             patch_transpose,
                             filter_width_ratio);

  // return an array with the same shape as the input
  return array_out.extract_slice(glm::ivec4(0, shape.x, 0, shape.y));
}

Array quilting_expand(const Array         &array,
                      float                expansion_ratio,
                      glm::ivec2           patch_base_shape,
                      float                overlap,
                      std::uint32_t        seed,
                      std::vector<Array *> secondary_arrays,
                      bool                 keep_input_shape,
                      bool                 patch_flip,
                      bool                 patch_rotate,
                      bool                 patch_transpose,
                      float                filter_width_ratio)
{
  if (!validate_non_empty(array) || !validate_shape(patch_base_shape))
    return Array();

  for (const auto *p : secondary_arrays)
    if (!p || !validate_non_empty(*p)) return Array();

  expansion_ratio = std::max(1.f, expansion_ratio);

  if (keep_input_shape)
  {
    // output shape is the same as the output
    glm::ivec2 work_shape = glm::ivec2(
        std::max(4, (int)(array.shape.x / expansion_ratio)),
        std::max(4, (int)(array.shape.y / expansion_ratio)));

    Array array_work = array.resample_to_shape(work_shape);

    // secondary arrays
    std::vector<Array>   secondary_arrays_storage = {};
    std::vector<Array *> secondary_arrays_work = {};

    for (auto v : secondary_arrays)
      secondary_arrays_storage.push_back(v->resample_to_shape(work_shape));

    for (auto &v : secondary_arrays_storage)
      secondary_arrays_work.push_back(&v);

    // apply
    glm::ivec2 patch_base_shape_work = glm::ivec2(
        std::max(2, (int)(patch_base_shape.x / expansion_ratio)),
        std::max(2, (int)(patch_base_shape.y / expansion_ratio)));

    // ensure patch base shape does not exceed work shape
    patch_base_shape_work.x = std::min(patch_base_shape_work.x, work_shape.x);
    patch_base_shape_work.y = std::min(patch_base_shape_work.y, work_shape.y);

    glm::ivec2 tiling = glm::ivec2(
        (int)(std::ceil((float)array.shape.x / patch_base_shape_work.x)),
        (int)(std::ceil((float)array.shape.y / patch_base_shape_work.y)));

    Array array_out = quilting({&array_work},
                               patch_base_shape_work,
                               tiling,
                               overlap,
                               seed,
                               secondary_arrays_work,
                               patch_flip,
                               patch_rotate,
                               patch_transpose,
                               filter_width_ratio);

    // override p_secondary_arrays content with output
    for (size_t k = 0; k < secondary_arrays.size(); k++)
      *secondary_arrays[k] = secondary_arrays_work[k]->extract_slice(
          glm::ivec4(0, array.shape.x, 0, array.shape.y));

    // return an array with the same shape as the input
    return array_out.extract_slice(
        glm::ivec4(0, array.shape.x, 0, array.shape.y));
  }
  else
  {
    // output shape is also expanded according to expansion factor
    glm::ivec2 expanded_shape = glm::ivec2(
        (int)(array.shape.x * expansion_ratio),
        (int)(array.shape.y * expansion_ratio));

    glm::ivec2 tiling = glm::ivec2(
        (int)(std::ceil(expanded_shape.x / patch_base_shape.x)),
        (int)(std::ceil(expanded_shape.y / patch_base_shape.y)));

    std::vector<const Array *> input = {&array};

    Array array_out = quilting(input,
                               patch_base_shape,
                               tiling,
                               overlap,
                               seed,
                               secondary_arrays,
                               patch_flip,
                               patch_rotate,
                               patch_transpose,
                               filter_width_ratio);

    // override p_secondary_arrays content with output
    for (auto v : secondary_arrays)
      *v = v->extract_slice(
          glm::ivec4(0, expanded_shape.x, 0, expanded_shape.y));

    // return an array with the expanded shape
    return array_out.extract_slice(
        glm::ivec4(0, expanded_shape.x, 0, expanded_shape.y));
  }
}

Array quilting_shuffle(const Array         &array,
                       glm::ivec2           patch_base_shape,
                       float                overlap,
                       std::uint32_t        seed,
                       std::vector<Array *> secondary_arrays,
                       bool                 patch_flip,
                       bool                 patch_rotate,
                       bool                 patch_transpose,
                       float                filter_width_ratio)
{
  if (!validate_non_empty(array) || !validate_shape(patch_base_shape))
    return Array();

  for (const auto *p : secondary_arrays)
    if (!p || !validate_non_empty(*p)) return Array();

  glm::ivec2 tiling = glm::ivec2(
      (int)(std::ceil(array.shape.x / patch_base_shape.x)),
      (int)(std::ceil(array.shape.y / patch_base_shape.y)));

  std::vector<const Array *> input = {&array};

  Array array_out = quilting(input,
                             patch_base_shape,
                             tiling,
                             overlap,
                             seed,
                             secondary_arrays,
                             patch_flip,
                             patch_rotate,
                             patch_transpose,
                             filter_width_ratio);

  // override p_secondary_arrays content with output
  for (auto v : secondary_arrays)
    *v = v->extract_slice(glm::ivec4(0, array.shape.x, 0, array.shape.y));

  // return an array with the same shape as the input
  return array_out.extract_slice(
      glm::ivec4(0, array.shape.x, 0, array.shape.y));
}

} // namespace hmap
