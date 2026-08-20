/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once
#include <cstddef>
#include <cstdint>

namespace hmap
{

/**
 * @brief Generates a fast deterministic uniform float in [0,1) using a 32-bit
 * hash.
 *
 * @param  seed Base seed value.
 * @param  k    Sample index.
 * @return      Pseudo-random float in [0,1).
 * @note Faster but lower quality than splitmix64_to_unit_float().
 */
float fast_hash32_to_unit_float(unsigned int seed, size_t k);

/**
 * @brief Computes a deterministic 64-bit hash using the SplitMix64 algorithm.
 *
 * This function is used to generate reproducible pseudo-random values from
 * integer inputs without maintaining any internal state. It is suitable for
 * generating deterministic per-sample random seeds.
 *
 * @param  x Input 64-bit value.
 * @return   64-bit hashed value.
 */
uint64_t splitmix64(uint64_t x);

/**
 * @brief Generates a deterministic uniform float in [0,1) from a seed and
 * sample index.
 *
 * The function combines the base seed and sample index, hashes the result using
 * SplitMix64, and converts the generated bits into a floating-point value in
 * the range [0,1). The output is deterministic and independent for each sample
 * index.
 *
 * @param  seed Base seed value.
 * @param  k    Sample index.
 * @return      Pseudo-random float in [0,1).
 */
float splitmix64_to_unit_float(unsigned int seed, size_t k);

/**
 * @brief Converts a 64-bit hash value into a deterministic uniform float in [0,
 * 1).
 *
 * Extracts the upper 24 bits of the input hash and scales them to produce a
 * single-precision floating-point value. The 24-bit resolution matches the
 * mantissa precision of a `float`, providing a deterministic pseudo-random
 * uniform distribution.
 *
 * @param  h 64-bit hash value.
 * @return   Uniform floating-point value in the range [0, 1).
 */
float uniform01(uint64_t h);

// === PdfSampler class ===

/**
 * @brief Samples indices from a discrete probability distribution.
 */
class PdfSampler
{
public:
  /**
   * @brief Builds the sampler from a PDF and seed.
   * @param pdf  Probability weights.
   * @param seed Random generator seed.
   */
  PdfSampler(const std::vector<float> &pdf, uint32_t seed);

  /**
   * @brief Samples a float value in [0, 1[.
   * @return Sampled float.
   */
  float sample();

  /**
   * @brief Samples multiple float values.
   * @param  nb_samples Number of samples.
   * @return            Vector of sampled values.
   */
  std::vector<float> sample(size_t nb_samples);

private:
  std::vector<float>                    cdf;
  std::mt19937                          generator;
  std::uniform_real_distribution<float> distribution{0.f, 1.f};
};

} // namespace hmap