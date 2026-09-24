#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl(); // for hydraulic erosion

  glm::ivec2 shape = {512, 512};
  int        seed = 42;

  // input terrain with natural features and rich erosion details
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                  shape,
                                  {4.f, 4.f},
                                  seed);
  hmap::gpu::hydraulic_particle(z, 30000, seed);
  hmap::remap(z);

  // loose symmetry parameters
  int   patch_size = 16;
  int   analysis_stride = patch_size / 8;
  int   synthesis_stride = patch_size / 2;
  int   factor = 4;
  float strength = 0.8f;

  // compute loose symmetry for all symmetry types
  hmap::Array z_lr = hmap::loose_symmetry(
      z,
      hmap::SymmetryType::SYMMETRY_LEFT_TO_RIGHT,
      strength,
      factor,
      patch_size,
      analysis_stride,
      synthesis_stride);

  hmap::Array z_rl = hmap::loose_symmetry(
      z,
      hmap::SymmetryType::SYMMETRY_RIGHT_TO_LEFT,
      strength,
      factor,
      patch_size,
      analysis_stride,
      synthesis_stride);

  hmap::Array z_tb = hmap::loose_symmetry(
      z,
      hmap::SymmetryType::SYMMETRY_TOP_TO_BOTTOM,
      strength,
      factor,
      patch_size,
      analysis_stride,
      synthesis_stride);

  hmap::Array z_bt = hmap::loose_symmetry(
      z,
      hmap::SymmetryType::SYMMETRY_BOTTOM_TO_TOP,
      strength,
      factor,
      patch_size,
      analysis_stride,
      synthesis_stride);

  hmap::Array z_x = hmap::loose_symmetry(z,
                                         hmap::SymmetryType::SYMMETRY_X,
                                         strength,
                                         factor,
                                         patch_size,
                                         analysis_stride,
                                         synthesis_stride);

  hmap::Array z_y = hmap::loose_symmetry(z,
                                         hmap::SymmetryType::SYMMETRY_Y,
                                         strength,
                                         factor,
                                         patch_size,
                                         analysis_stride,
                                         synthesis_stride);

  hmap::Array z_xy = hmap::loose_symmetry(z,
                                          hmap::SymmetryType::SYMMETRY_XY,
                                          strength,
                                          factor,
                                          patch_size,
                                          analysis_stride,
                                          synthesis_stride);

  hmap::Array z_rot180 = hmap::loose_symmetry(
      z,
      hmap::SymmetryType::SYMMETRY_ROT180,
      strength,
      factor,
      patch_size,
      analysis_stride,
      synthesis_stride);

  // loose symmetry with flatten_center enabled
  hmap::Array z_lr_flat = hmap::loose_symmetry(
      z,
      hmap::SymmetryType::SYMMETRY_LEFT_TO_RIGHT,
      strength,
      factor,
      patch_size,
      analysis_stride,
      synthesis_stride,
      /* sparsity */ 1,
      /* flatten_center */ true,
      /* flatten_radius */ 0.08f);

  hmap::export_banner_png(
      "ex_loose_symmetry.png",
      {z, z_lr, z_rl, z_tb, z_bt, z_x, z_y, z_xy, z_rot180, z_lr_flat},
      hmap::Cmap::TERRAIN,
      true);
}
