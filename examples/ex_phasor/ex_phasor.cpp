#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  shape = {1024, 1024};
  float         kp = 16.f;
  std::uint32_t seed = 0;

  std::vector<hmap::PhasorProfile> profiles = {
      hmap::PhasorProfile::PP_COSINE_BULKY,
      hmap::PhasorProfile::PP_COSINE_PEAKY,
      hmap::PhasorProfile::PP_COSINE_STD,
      hmap::PhasorProfile::PP_COSINE_SQUARE,
      hmap::PhasorProfile::PP_TRIANGLE,
  };

  std::vector<hmap::Array> arrays;

  for (auto p : profiles)
  {
    hmap::Array z = hmap::gpu::phasor(p,
                                      shape,
                                      kp,
                                      seed,
                                      /* angle_shift */ 15.f);
    arrays.push_back(z);
  }

  // fbm phasor
  hmap::Array z = hmap::gpu::phasor_fbm(hmap::PhasorProfile::PP_COSINE_STD,
                                        shape,
                                        kp,
                                        seed,
                                        /* angle_shift */ 15.f);
  arrays.push_back(z);

  hmap::export_banner_png("ex_phasor.png",
                          arrays,
                          hmap::Cmap::JET,
                          false,
                          true);
}
