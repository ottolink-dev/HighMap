/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "highmap/flora/forest_growth.hpp"
#include "highmap/flora/species.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

// ============================================================================
//  InteractionMatrix Implementation
// ============================================================================

InteractionMatrix::InteractionMatrix(size_t num_species, float default_val)
    : size(num_species), values(num_species * num_species, default_val)
{
}

void InteractionMatrix::fill(float value)
{
  std::fill(values.begin(), values.end(), value);
}

float InteractionMatrix::get(size_t s1, size_t s2) const
{
  if (s1 >= size || s2 >= size) return 0.0f;
  return values[s1 * size + s2];
}

void InteractionMatrix::set(size_t s1, size_t s2, float value)
{
  if (s1 < size && s2 < size) values[s1 * size + s2] = value;
}

void InteractionMatrix::set_symmetric(size_t s1, size_t s2, float value)
{
  set(s1, s2, value);
  set(s2, s1, value);
}

InteractionMatrix InteractionMatrix::diagonal(const std::vector<float> &diag,
                                              float off_diag_val)
{
  InteractionMatrix mat(diag.size(), off_diag_val);
  for (size_t i = 0; i < diag.size(); ++i)
    mat.set(i, i, diag[i]);
  return mat;
}

InteractionMatrix InteractionMatrix::from_radii(const std::vector<float> &radii,
                                                float multiplier)
{
  size_t            num_species = radii.size();
  InteractionMatrix mat(num_species);

  for (size_t s1 = 0; s1 < num_species; ++s1)
  {
    for (size_t s2 = 0; s2 < num_species; ++s2)
    {
      mat.set(s1, s2, multiplier * (radii[s1] + radii[s2]));
    }
  }

  return mat;
}

InteractionMatrix InteractionMatrix::from_species(
    const std::vector<Species> &species,
    float                       multiplier)
{
  size_t            num_species = species.size();
  InteractionMatrix mat(num_species);

  for (size_t s1 = 0; s1 < num_species; ++s1)
  {
    for (size_t s2 = 0; s2 < num_species; ++s2)
    {
      mat.set(s1, s2, multiplier * (species[s1].radius + species[s2].radius));
    }
  }

  return mat;
}

InteractionMatrix InteractionMatrix::random(size_t   num_species,
                                            uint32_t seed,
                                            float    random_offset,
                                            bool     symmetric)
{
  InteractionMatrix                     mat(num_species);
  std::mt19937                          rng(seed);
  std::uniform_real_distribution<float> dist(1.0f - random_offset,
                                             1.0f + random_offset);

  if (symmetric)
  {
    for (size_t i = 0; i < num_species; ++i)
    {
      mat.set(i, i, dist(rng));
      for (size_t j = i + 1; j < num_species; ++j)
      {
        float val = dist(rng);
        mat.set_symmetric(i, j, val);
      }
    }
  }
  else
  {
    for (size_t i = 0; i < num_species; ++i)
    {
      for (size_t j = 0; j < num_species; ++j)
      {
        mat.set(i, j, dist(rng));
      }
    }
  }

  return mat;
}

InteractionMatrix InteractionMatrix::uniform(size_t num_species, float val)
{
  return InteractionMatrix(num_species, val);
}

} // namespace hmap
