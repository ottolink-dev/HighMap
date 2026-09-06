#include "highmap.hpp"
#include "highmap/dbg/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(FlowSimulation, NonNegativity)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array depth_map(shape, 1.f);
  Array water_depth = gpu::flow_simulation(z, 0.05f, depth_map, 50);

  // water depth must remain non-negative everywhere
  EXPECT_GE(water_depth.min(), 0.f);
}

TEST(FlowSimulationViscous, MassConservationFlatDomain)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z(shape, 0.f); // flat terrain

  // localized water droplet in center
  Array depth_map =
      disk(shape, 0.2f, 8.f, nullptr, nullptr, nullptr, {0.5f, 0.5f});
  float water_height = 0.1f;
  Array d_init = water_height * depth_map;
  float initial_volume = d_init.sum();

  Array water_depth = gpu::flow_simulation_viscous(z,
                                                   water_height,
                                                   depth_map,
                                                   100,
                                                   0.1f,
                                                   0.f, // no dry-out
                                                   1.f,
                                                   2.5f);

  // non-negativity
  EXPECT_GE(water_depth.min(), 0.f);

  // mass conservation on flat domain with central droplet away from boundaries
  float final_volume = water_depth.sum();
  EXPECT_NEAR(initial_volume, final_volume, initial_volume * 1e-2f);
}

TEST(FlowSimulationViscous, NoSpuriousMassFromDryRidge)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z(shape, 0.f);

  // tall dry ridge on left half
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x / 2; ++i)
    {
      z(i, j) = 10.f;
    }
  }

  // water pool only on right half (depression/valley)
  Array depth_map(shape, 0.f);
  for (int j = 16; j < 48; ++j)
  {
    for (int i = 36; i < 60; ++i)
    {
      depth_map(i, j) = 1.f;
    }
  }

  float water_height = 0.05f;
  Array d_init = water_height * depth_map;
  float initial_volume = d_init.sum();

  Array water_depth = gpu::flow_simulation_viscous(z,
                                                   water_height,
                                                   depth_map,
                                                   100,
                                                   0.1f,
                                                   0.f,
                                                   1.f,
                                                   2.5f);

  // dry ridge on left must not have gained phantom fluid or caused total volume
  // explosion
  float final_volume = water_depth.sum();
  EXPECT_LE(final_volume, initial_volume * 1.05f);
  EXPECT_GE(water_depth.min(), 0.f);
}

TEST(FlowSimulationViscous, IterationParityAndDeterminism)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array       depth_map(shape, 1.f);
  const float water_height = 0.05f;

  Array d_even = gpu::flow_simulation_viscous(z, water_height, depth_map, 50);
  Array d_odd = gpu::flow_simulation_viscous(z, water_height, depth_map, 51);
  Array d_even_again = gpu::flow_simulation_viscous(z,
                                                    water_height,
                                                    depth_map,
                                                    50);

  // deterministic
  EXPECT_EQ(d_even.vector, d_even_again.vector);

  // evolution check
  const float s_even = d_even.sum();
  const float s_odd = d_odd.sum();
  EXPECT_GT(s_even, 0.f);
  EXPECT_NEAR(s_odd, s_even, 0.05f * s_even);
  EXPECT_NE(d_even.vector, (water_height * depth_map).vector);
  EXPECT_GE(d_odd.min(), 0.f);
}

TEST(FlowSimulationViscous, ZeroIterationsReturnsInitialState)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {32, 48};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 3);
  remap(z);

  Array depth_map = noise_fbm(NoiseType::PERLIN, shape, {4.f, 4.f}, 5);
  remap(depth_map);

  Array d = gpu::flow_simulation_viscous(z, 0.1f, depth_map, 0);

  EXPECT_EQ(d.vector, (0.1f * depth_map).vector);
}

TEST(FlowSimulationViscous, OutflowBoundaries)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z(shape, 0.f);

  // tilt terrain down towards right (+x)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      z(i, j) = 1.f - (float)i / (float)shape.x;

  Array depth_map(shape, 1.f);
  float water_height = 0.05f;

  Array d_closed = gpu::flow_simulation_viscous(z,
                                                water_height,
                                                depth_map,
                                                100,
                                                0.1f,
                                                0.f,
                                                1.f,
                                                2.5f,
                                                0.f,
                                                false);

  Array d_outflow = gpu::flow_simulation_viscous(z,
                                                 water_height,
                                                 depth_map,
                                                 100,
                                                 0.1f,
                                                 0.f,
                                                 1.f,
                                                 2.5f,
                                                 0.f,
                                                 true);

  EXPECT_LT(d_outflow.sum(), d_closed.sum());
}

TEST(FlowSimulationViscous, Evaporation)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z(shape, 0.f);
  Array            depth_map(shape, 1.f);
  float            water_height = 0.1f;

  Array d_no_evap = gpu::flow_simulation_viscous(z,
                                                 water_height,
                                                 depth_map,
                                                 50,
                                                 0.1f,
                                                 0.f,
                                                 1.f,
                                                 2.5f,
                                                 0.f);

  Array d_evap = gpu::flow_simulation_viscous(z,
                                              water_height,
                                              depth_map,
                                              50,
                                              0.1f,
                                              0.f,
                                              1.f,
                                              2.5f,
                                              0.1f);

  EXPECT_LT(d_evap.sum(), d_no_evap.sum());
}

TEST(FlowSimulation, VelocityExport)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array depth_map(shape, 1.f);
  Array u(shape);
  Array v(shape);

  Array water_depth = gpu::flow_simulation(z,
                                           0.05f,
                                           depth_map,
                                           20,
                                           0.5f,
                                           true,
                                           0.01f,
                                           0.f,
                                           nullptr,
                                           0.f,
                                           0.f,
                                           false,
                                           &u,
                                           &v);

  EXPECT_EQ(u.shape, shape);
  EXPECT_EQ(v.shape, shape);
  EXPECT_GE(water_depth.min(), 0.f);
}

TEST(FlowSimulation, OutflowBoundaries)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z(shape, 0.f);

  // tilt terrain down towards right (+x)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      z(i, j) = 1.f - (float)i / (float)shape.x;

  Array depth_map(shape, 1.f);
  float water_height = 0.05f;

  // run with closed boundaries
  Array d_closed = gpu::flow_simulation(z,
                                        water_height,
                                        depth_map,
                                        100,
                                        0.5f,
                                        true,
                                        0.01f,
                                        0.f,
                                        nullptr,
                                        0.f,
                                        0.f,
                                        false);

  // run with outflow boundaries
  Array d_outflow = gpu::flow_simulation(z,
                                         water_height,
                                         depth_map,
                                         100,
                                         0.5f,
                                         true,
                                         0.01f,
                                         0.f,
                                         nullptr,
                                         0.f,
                                         0.f,
                                         true);

  // outflow boundaries must allow water to leave the domain
  EXPECT_LT(d_outflow.sum(), d_closed.sum());
}

TEST(FlowSimulation, RainAndEvaporation)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z(shape, 0.f);
  Array            depth_map(shape, 0.f); // start dry

  // simulation with rain
  Array d_rain = gpu::flow_simulation(z,
                                      1.f,
                                      depth_map,
                                      50,
                                      0.1f,
                                      false,
                                      0.f,
                                      0.f,
                                      nullptr,
                                      0.01f, // rain rate
                                      0.f);

  EXPECT_GT(d_rain.sum(), 0.f);

  // simulation with evaporation on wet terrain
  Array d_wet(shape, 1.f);
  Array d_evap = gpu::flow_simulation(z,
                                      0.1f,
                                      d_wet,
                                      50,
                                      0.1f,
                                      false,
                                      0.f,
                                      0.f,
                                      nullptr,
                                      0.f,
                                      0.1f); // evap rate

  EXPECT_LT(d_evap.sum(), (0.1f * d_wet).sum());
}

TEST(FlowSimulation, IterationParityAndDeterminism)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array       depth_map(shape, 1.f);
  const float water_height = 0.05f;

  Array d_even = gpu::flow_simulation(z, water_height, depth_map, 50);
  Array d_odd = gpu::flow_simulation(z, water_height, depth_map, 51);
  Array d_even_again = gpu::flow_simulation(z, water_height, depth_map, 50);

  // deterministic
  EXPECT_EQ(d_even.vector, d_even_again.vector);

  // one more step changes little, but the field must have evolved
  const float s_even = d_even.sum();
  const float s_odd = d_odd.sum();
  EXPECT_GT(s_even, 0.f);
  EXPECT_NEAR(s_odd, s_even, 0.05f * s_even);
  EXPECT_NE(d_even.vector, (water_height * depth_map).vector);
  EXPECT_GE(d_odd.min(), 0.f);
}

TEST(FlowSimulation, ZeroIterationsReturnsInitialState)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {32, 48};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 3);
  remap(z);

  Array depth_map = noise_fbm(NoiseType::PERLIN, shape, {4.f, 4.f}, 5);
  remap(depth_map);

  Array u(shape, 1.f);
  Array v(shape, 1.f);

  Array d = gpu::flow_simulation(z,
                                 0.1f,
                                 depth_map,
                                 0,
                                 0.5f,
                                 true,
                                 0.01f,
                                 0.f,
                                 nullptr,
                                 0.f,
                                 0.f,
                                 false,
                                 &u,
                                 &v);

  EXPECT_EQ(d.vector, (0.1f * depth_map).vector);
  EXPECT_EQ(u.shape, shape);
  EXPECT_EQ(v.shape, shape);
  EXPECT_FLOAT_EQ(u.max(), 0.f);
  EXPECT_FLOAT_EQ(v.max(), 0.f);
}

TEST(FlowSimulation, RainMapOfOnesEqualsUniformRain)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array depth_map(shape, 0.f); // start dry
  Array ones(shape, 1.f);

  Array d_map = gpu::flow_simulation(z,
                                     1.f,
                                     depth_map,
                                     30,
                                     0.1f,
                                     false,
                                     0.f,
                                     0.f,
                                     &ones,
                                     0.02f,
                                     0.f,
                                     false);
  Array d_uniform = gpu::flow_simulation(z,
                                         1.f,
                                         depth_map,
                                         30,
                                         0.1f,
                                         false,
                                         0.f,
                                         0.f,
                                         nullptr,
                                         0.02f,
                                         0.f,
                                         false);

  EXPECT_GT(d_map.sum(), 0.f);
  EXPECT_EQ(d_map.vector, d_uniform.vector);
}
