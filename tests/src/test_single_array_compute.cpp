#include <cstdint>
#include <sstream>
#include <tuple>

#include "highmap/virtual_array/virtual_array.hpp"

#include <bit>
#include <gtest/gtest.h>

using namespace hmap;

namespace
{
const ComputeMode sequential{.mode = ForEachMode::VA_SEQUENTIAL};

Array pattern(glm::ivec2 shape, float offset)
{
  Array result(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      result(i, j) = offset + float((i * 17 + j * 31 + i * j * 3) % 97) / 13.f;
  return result;
}

void expect_same(const Array &actual, const Array &expected)
{
  ASSERT_EQ(actual.shape, expected.shape);
  EXPECT_EQ(actual.vector, expected.vector);
}

// XML properties allow the standalone baseline/fixed builds to compare exact
// output fingerprints in addition to the per-value legacy-contract assertions.
std::string fingerprint(const Array &array)
{
  std::uint64_t value = 14695981039346656037ULL;
  for (float sample : array.vector)
  {
    value ^= std::bit_cast<std::uint32_t>(sample);
    value *= 1099511628211ULL;
  }
  std::ostringstream result;
  result << std::hex << value;
  return result.str();
}

using Parameters = std::tuple<StorageMode, int, bool, int, bool>;
class SingleArraySnapshotTest : public ::testing::TestWithParam<Parameters>
{
};

TEST_P(SingleArraySnapshotTest, PreservesSnapshotsOptionalPortsAndOutputReads)
{
  auto [storage, halo, tiled, stride, small] = GetParam();
  const glm::ivec2 shape = small ? glm::ivec2(3, 2) : glm::ivec2(17, 13);
  const glm::ivec2 tile = tiled ? glm::ivec2(2, 2) : shape;
  VirtualArray     a(shape, tile, halo, storage);
  VirtualArray     b(shape, tile, halo, storage);
  VirtualArray     c(shape, tile, halo, storage);
  a.from_array(pattern(shape, 2.f), sequential);
  b.from_array(pattern(shape, 11.f), sequential);
  c.from_array(pattern(shape, -7.f), sequential);

  const glm::ivec2 working =
      stride > 1 ? glm::ivec2(std::max(2, int(shape.x / float(stride))),
                              std::max(2, int(shape.y / float(stride))))
                 : shape;
  const Array input_a = a.to_array(sequential).resample_to_shape(working);
  const Array input_b = b.to_array(sequential).resample_to_shape(working);
  const Array old_c = c.to_array(sequential).resample_to_shape(working);
  Array       expected_a = input_a;
  Array       expected_c = old_c;
  for (int k = 0; k < expected_a.size(); ++k)
  {
    expected_a.vector[k] += 2.f * input_a.vector[k] - input_b.vector[k];
    expected_c.vector[k] += 13.f;
  }

  // Inputs and outputs intentionally alias. Duplicate outputs must each start
  // from their original value; scattering the last slot wins as before.
  TileAccess access{{&a, nullptr, &a, &b}, {&a, &c, nullptr, &c}};
  int        calls = 0;
  for_each_tile(access,
                [&](const std::vector<const Array *> &in,
                    std::vector<Array *>             &out,
                    const UvTile                     &region)
                {
                  ++calls;
                  EXPECT_EQ(region.shape, working);
                  ASSERT_EQ(in.size(), 4u);
                  ASSERT_EQ(out.size(), 4u);
                  EXPECT_EQ(in[1], nullptr);
                  EXPECT_EQ(out[2], nullptr);
                  EXPECT_NE(in[0], in[2]);
                  EXPECT_NE(in[0], out[0]);
                  EXPECT_NE(out[1], out[3]);
                  expect_same(*in[0], input_a);
                  expect_same(*in[2], input_a);
                  expect_same(*in[3], input_b);
                  expect_same(*out[0], input_a);
                  expect_same(*out[1], old_c);
                  expect_same(*out[3], old_c);
                  for (int k = 0; k < out[0]->size(); ++k)
                  {
                    out[0]->vector[k] += 2.f * in[0]->vector[k] -
                                         in[3]->vector[k];
                    out[1]->vector[k] += 7.f;
                    out[3]->vector[k] += 13.f;
                  }
                  // Writing an aliased output must not change any input
                  // snapshot.
                  expect_same(*in[0], input_a);
                  expect_same(*in[2], input_a);
                },
                {.mode = ForEachMode::VA_SINGLE_ARRAY, .stride = stride});
  ASSERT_EQ(calls, 1);

  VirtualArray oracle_a(shape, tile, halo, storage);
  VirtualArray oracle_c(shape, tile, halo, storage);
  if (working == shape)
  {
    oracle_a.from_array(expected_a, sequential);
    oracle_c.from_array(expected_c, sequential);
  }
  else
  {
    oracle_a.from_array_bicubic(expected_a, sequential);
    oracle_c.from_array_bicubic(expected_c, sequential);
  }
  const Array actual_a = a.to_array(sequential);
  const Array actual_c = c.to_array(sequential);
  expect_same(actual_a, oracle_a.to_array(sequential));
  expect_same(actual_c, oracle_c.to_array(sequential));
  expect_same(b.to_array(sequential), pattern(shape, 11.f));
  RecordProperty("output_a_fingerprint", fingerprint(actual_a));
  RecordProperty("output_c_fingerprint", fingerprint(actual_c));
}

INSTANTIATE_TEST_SUITE_P(
    StorageGridHaloStride,
    SingleArraySnapshotTest,
    ::testing::Combine(::testing::Values(StorageMode::VA_RAM,
                                         StorageMode::VA_DISK_LRU_STREAM,
                                         StorageMode::VA_DISK_LRU_FULL),
                       ::testing::Values(0, 1),
                       ::testing::Bool(),
                       ::testing::Values(1, 2, 3),
                       ::testing::Bool()));

TEST(SingleArrayComputeTest, OutputOnlyCommunicatorKeepsExistingValues)
{
  VirtualArray output({11, 9}, {4, 3}, 1, StorageMode::VA_DISK_LRU_STREAM);
  output.from_array(pattern(output.shape, 8.f), sequential);
  Array expected = output.to_array(sequential);
  for (float &value : expected.vector)
    value += 3.f;
  int calls = 0;
  for_each_tile(output,
                [&](Array &array, const UvTile &, ICommunicator &comm)
                {
                  ++calls;
                  EXPECT_EQ(comm.size(), 1);
                  for (float &value : array.vector)
                    value += 3.f;
                },
                {.mode = ForEachMode::VA_SINGLE_ARRAY});
  ASSERT_EQ(calls, 1);
  expect_same(output.to_array(sequential), expected);
}

TEST(SingleArrayComputeTest, DirectMixedShapesRetainLegacyResampling)
{
  VirtualArray input({19, 11}, {6, 4}, 0, StorageMode::VA_RAM);
  VirtualArray output({13, 17}, {5, 6}, 1, StorageMode::VA_DISK_LRU_STREAM);
  input.from_array(pattern(input.shape, 2.f), sequential);
  output.from_array(pattern(output.shape, 5.f), sequential);
  const Array expected_input = input.to_array(sequential)
                                   .resample_to_shape(output.shape);
  const Array expected_output = output.to_array(sequential);
  int         calls = 0;
  single_array_compute(TileAccess{{&input}, {&output}},
                       [&](const std::vector<const Array *> &in,
                           std::vector<Array *>             &out,
                           const UvTile &)
                       {
                         ++calls;
                         expect_same(*in[0], expected_input);
                         expect_same(*out[0], expected_output);
                         *out[0] = *in[0];
                       });
  ASSERT_EQ(calls, 1);
  expect_same(output.to_array(sequential), expected_input);
}
} // namespace
