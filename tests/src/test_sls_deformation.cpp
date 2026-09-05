#include "highmap/authoring.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"

#include <gtest/gtest.h>

using namespace hmap;

namespace
{

float helper_mse(const Array &a, const Array &b)
{
  float sum = 0.f;
  for (size_t k = 0; k < a.vector.size(); ++k)
  {
    float d = a.vector[k] - b.vector[k];
    sum += d * d;
  }
  return sum / static_cast<float>(a.vector.size());
}

// fraction of vertices on the right side of the water level w.r.t. the mask
float helper_mask_agreement(const Array &z, const Array &mask, float level)
{
  int count = 0;
  for (size_t k = 0; k < z.vector.size(); ++k)
  {
    bool land = z.vector[k] > level;
    bool want = mask.vector[k] > 0.5f;
    if (land == want) ++count;
  }
  return static_cast<float>(count) / static_cast<float>(z.vector.size());
}

Array helper_disk_mask(glm::ivec2 shape, float radius_ratio)
{
  Array mask(shape, 0.f);
  float cx = 0.5f * static_cast<float>(shape.x - 1);
  float cy = 0.5f * static_cast<float>(shape.y - 1);
  float r = radius_ratio * static_cast<float>(std::min(shape.x, shape.y));
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = static_cast<float>(i) - cx;
      float dy = static_cast<float>(j) - cy;
      if (dx * dx + dy * dy < r * r) mask(i, j) = 1.f;
    }
  return mask;
}

} // namespace

TEST(SlsDeformation, GaussianPushSetsCenterAndVanishesOutsideRadius)
{
  Array z({32, 32}, 0.f);

  std::vector<GaussianPush> pushes = {{10, 10, 4, 0.5f}};
  apply_gaussian_pushes(z, pushes);

  EXPECT_FLOAT_EQ(z(10, 10), 0.5f);
  EXPECT_GT(z(10, 12), 0.f);
  EXPECT_LT(z(10, 12), 0.5f);
  EXPECT_FLOAT_EQ(z(10, 15), 0.f);
  EXPECT_FLOAT_EQ(z(15, 10), 0.f);
  EXPECT_FLOAT_EQ(z(0, 0), 0.f);
}

TEST(SlsDeformation, GaussianPushIsClippedAtArrayBorder)
{
  Array z({16, 16}, 0.f);

  std::vector<GaussianPush> pushes = {{0, 0, 8, 1.f}};
  apply_gaussian_pushes(z, pushes);

  EXPECT_FLOAT_EQ(z(0, 0), 1.f);
  EXPECT_GT(z(3, 3), 0.f);
}

TEST(SlsDeformation, MatchConstraintReducesDistanceToTarget)
{
  glm::ivec2 shape = {64, 64};

  Array z(shape, 0.f);
  Array target = gaussian_pulse(shape, 12.f);

  DeformationConstraint c;
  c.target = target;
  c.weight = Array(shape, 1.f);
  c.type = DeformationConstraintType::MATCH;

  float mse_before = helper_mse(z, target);

  Array z_out = sls_deformation(z, {c}, 1u, 100, 2, 16, 4);

  ASSERT_EQ(z_out.shape, shape);
  float mse_after = helper_mse(z_out, target);

  EXPECT_LT(mse_after, 0.25f * mse_before);
}

TEST(SlsDeformation, ShapeMaskConstraintsEnforceLandAndSea)
{
  glm::ivec2 shape = {64, 64};
  float      level = 0.5f;

  Array z = noise(NoiseType::PERLIN, shape, {4.f, 4.f}, 3u);
  remap(z, 0.f, 1.f);

  Array mask = helper_disk_mask(shape, 0.3f);
  Array sea(shape, 1.f);
  sea -= mask;

  DeformationConstraint land_c;
  land_c.target = Array(shape, level);
  land_c.weight = mask;
  land_c.type = DeformationConstraintType::ABOVE;

  DeformationConstraint sea_c;
  sea_c.target = Array(shape, level);
  sea_c.weight = sea;
  sea_c.type = DeformationConstraintType::BELOW;

  float agreement_before = helper_mask_agreement(z, mask, level);

  Array z_out = sls_deformation(z, {land_c, sea_c}, 7u, 200, 2, 16, 4);

  float agreement_after = helper_mask_agreement(z_out, mask, level);

  EXPECT_GT(agreement_after, agreement_before);
  EXPECT_GT(agreement_after, 0.9f);
}

TEST(SlsDeformation, RecordedPushesReplayExactly)
{
  glm::ivec2 shape = {48, 48};

  Array z = noise(NoiseType::PERLIN, shape, {3.f, 3.f}, 5u);
  Array target = gaussian_pulse(shape, 10.f);

  DeformationConstraint c;
  c.target = target;
  c.weight = Array(shape, 1.f);

  std::vector<GaussianPush> pushes;

  Array z_out = sls_deformation(z,
                                {c},
                                11u,
                                50,
                                2,
                                12,
                                3,
                                0.f,
                                256,
                                0.65f,
                                0.1f,
                                1e-3f,
                                &pushes);

  EXPECT_FALSE(pushes.empty());

  Array z_replay = z;
  apply_gaussian_pushes(z_replay, pushes);

  ASSERT_EQ(z_replay.vector.size(), z_out.vector.size());
  for (size_t k = 0; k < z_out.vector.size(); ++k)
    EXPECT_FLOAT_EQ(z_replay.vector[k], z_out.vector[k]);
}

TEST(SlsDeformation, SameSeedIsDeterministic)
{
  glm::ivec2 shape = {48, 48};

  Array z = noise(NoiseType::PERLIN, shape, {3.f, 3.f}, 5u);

  DeformationConstraint c;
  c.target = gaussian_pulse(shape, 10.f);
  c.weight = Array(shape, 1.f);

  Array z1 = sls_deformation(z, {c}, 42u, 30, 2, 12, 3);
  Array z2 = sls_deformation(z, {c}, 42u, 30, 2, 12, 3);

  for (size_t k = 0; k < z1.vector.size(); ++k)
    EXPECT_FLOAT_EQ(z1.vector[k], z2.vector[k]);
}

TEST(SlsDeformation, SatisfiedConstraintLeavesTerrainUnchanged)
{
  glm::ivec2 shape = {32, 32};

  Array z = noise(NoiseType::PERLIN, shape, {2.f, 2.f}, 9u);

  DeformationConstraint c;
  c.target = z;
  c.weight = Array(shape, 1.f);

  std::vector<GaussianPush> pushes;

  Array z_out = sls_deformation(z,
                                {c},
                                1u,
                                20,
                                2,
                                8,
                                3,
                                0.f,
                                128,
                                0.65f,
                                0.1f,
                                1e-3f,
                                &pushes);

  EXPECT_TRUE(pushes.empty());
  for (size_t k = 0; k < z.vector.size(); ++k)
    EXPECT_FLOAT_EQ(z_out.vector[k], z.vector[k]);
}

TEST(SlsDeformation, SlopeLimitBoundsPushAmplitude)
{
  glm::ivec2 shape = {48, 48};
  float      talus_max = 0.01f;

  Array z(shape, 0.f);

  DeformationConstraint c;
  c.target = gaussian_pulse(shape, 10.f);
  c.weight = Array(shape, 1.f);

  std::vector<GaussianPush> pushes;
  sls_deformation(z,
                  {c},
                  3u,
                  40,
                  2,
                  12,
                  3,
                  talus_max,
                  256,
                  0.65f,
                  0.1f,
                  1e-3f,
                  &pushes);

  EXPECT_FALSE(pushes.empty());
  for (auto &p : pushes)
    EXPECT_LE(std::abs(p.amplitude),
              talus_max * static_cast<float>(p.ir) + 1e-6f);
}

TEST(SlsDeformation, EmptyConstraintsReturnInput)
{
  glm::ivec2 shape = {16, 16};
  Array      z = noise(NoiseType::PERLIN, shape, {2.f, 2.f}, 1u);

  std::vector<GaussianPush> pushes;

  Array z_out = sls_deformation(z,
                                {},
                                1u,
                                10,
                                2,
                                4,
                                2,
                                0.f,
                                64,
                                0.65f,
                                0.1f,
                                1e-3f,
                                &pushes);

  EXPECT_TRUE(pushes.empty());
  for (size_t k = 0; k < z.vector.size(); ++k)
    EXPECT_FLOAT_EQ(z_out.vector[k], z.vector[k]);
}

TEST(SlsDeformation, MismatchedConstraintShapeReturnsInput)
{
  glm::ivec2 shape = {16, 16};
  Array      z = noise(NoiseType::PERLIN, shape, {2.f, 2.f}, 1u);

  DeformationConstraint c;
  c.target = Array({8, 8}, 1.f);
  c.weight = Array({8, 8}, 1.f);

  std::vector<GaussianPush> pushes;

  Array z_out = sls_deformation(z,
                                {c},
                                1u,
                                10,
                                2,
                                4,
                                2,
                                0.f,
                                64,
                                0.65f,
                                0.1f,
                                1e-3f,
                                &pushes);

  EXPECT_TRUE(pushes.empty());
  for (size_t k = 0; k < z.vector.size(); ++k)
    EXPECT_FLOAT_EQ(z_out.vector[k], z.vector[k]);
}

TEST(SlsDeformation, ConstraintScaleWeightsCompetingTerms)
{
  glm::ivec2 shape = {48, 48};

  Array z(shape, 0.f);

  // two competing full-weight targets: a bump and a flat plane
  DeformationConstraint bump;
  bump.target = gaussian_pulse(shape, 10.f);
  bump.weight = Array(shape, 1.f);

  DeformationConstraint flat;
  flat.target = Array(shape, 0.f);
  flat.weight = Array(shape, 1.f);

  // heavily favour the bump
  bump.scale = 100.f;
  flat.scale = 1.f;
  Array z_bump = sls_deformation(z, {bump, flat}, 1u, 60, 2, 12, 3);

  // heavily favour the plane
  bump.scale = 1.f;
  flat.scale = 100.f;
  Array z_flat = sls_deformation(z, {bump, flat}, 1u, 60, 2, 12, 3);

  EXPECT_LT(helper_mse(z_bump, bump.target), helper_mse(z_flat, bump.target));
}
