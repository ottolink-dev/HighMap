#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};

  // --- Random crater field

  uint        seed = 1;
  int         n_craters = 25;
  float       r_min = 0.03f;
  float       r_max = 0.18f;
  hmap::Array z1 = hmap::crater_field(shape, seed, n_craters, r_min, r_max);

  // --- Explicit centers and radii with radial noise

  float     noise_amp = 0.15f;
  glm::vec2 kw = {4.f, 4.f};
  auto      noise = noise_amp *
               hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);

  std::vector<glm::vec2> centers = {{0.25f, 0.30f},
                                    {0.65f, 0.55f},
                                    {0.45f, 0.75f},
                                    {0.80f, 0.25f}};
  std::vector<float>     radii = {0.18f, 0.22f, 0.10f, 0.12f};

  hmap::Array crater_mask;
  hmap::Array z2 = hmap::crater_field(shape,
                                      centers,
                                      radii,
                                      /* angle */ -20.f,
                                      /* inner_depth */ 0.5f,
                                      /* inner_exp */ 3.f,
                                      /* lip_height */ 0.1f,
                                      /* lip_extent */ 0.25f,
                                      /* lip_exp */ 2.f,
                                      /* asym_ratio */ 0.6f,
                                      /* central_peak_height */ 0.05f,
                                      /* central_peak_extent */ 0.4f,
                                      /* n_terraces */ 0,
                                      /* terrace_extent */ 0.1f,
                                      /* terrace_exp */ 2.f,
                                      /* terrace_persistence */ 0.5f,
                                      &noise,
                                      /* bbox */ {0.f, 1.f, 0.f, 1.f},
                                      &crater_mask);

  hmap::export_banner_png("ex_crater_field.png",
                          {z1, z2, crater_mask},
                          hmap::Cmap::TERRAIN,
                          true,
                          true);
}
