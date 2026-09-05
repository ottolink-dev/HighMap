#include <cmath>
#include <iostream>

#include "highmap.hpp"

// closed star-shaped contour: radius r0 modulated by a fixed angular profile
// (all blobs share the profile so that nested blobs never cross)
hmap::Path blob(float cx, float cy, float r0, int npts = 96)
{
  std::vector<float> x, y;
  for (int k = 0; k < npts; ++k)
  {
    float t = 2.f * M_PI * k / npts;
    float m = 1.f + 0.12f * std::sin(3.f * t + 0.7f) +
              0.08f * std::sin(5.f * t + 2.1f);
    x.push_back(cx + r0 * m * std::cos(t));
    y.push_back(cy + r0 * m * std::sin(t));
  }
  hmap::Path p(x, y);
  p.set_closed(true);
  return p;
}

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};
  int        seed = 1;

  // --- sparse contours: a main hill with 4 nested levels, a small basin in
  // its lowest ring and a second, separate hill
  std::vector<hmap::Path> contours = {blob(0.38f, 0.48f, 0.30f),
                                      blob(0.38f, 0.48f, 0.21f),
                                      blob(0.38f, 0.48f, 0.13f),
                                      blob(0.38f, 0.48f, 0.06f),
                                      blob(0.67f, 0.48f, 0.03f),
                                      blob(0.82f, 0.82f, 0.11f),
                                      blob(0.82f, 0.82f, 0.04f)};
  std::vector<float> elevations = {0.1f, 0.3f, 0.5f, 0.7f, 0.f, 0.2f, 0.4f};

  // input contours drawn at their elevation, for reference
  hmap::Array z_contours(shape);
  for (size_t k = 0; k < contours.size(); ++k)
  {
    hmap::Path p = contours[k];
    p.set_values(elevations[k]);
    p.to_array(z_contours, bbox);
  }

  // --- deterministic front propagation (weighted geodesic interpolation)
  hmap::Array z1 = hmap::elevation_from_contours(shape,
                                                 contours,
                                                 elevations,
                                                 nullptr,
                                                 0.f);

  // --- stochastic Eden growth, uniform probability
  hmap::Array z2 = hmap::elevation_from_contours(shape,
                                                 contours,
                                                 elevations,
                                                 nullptr,
                                                 1.f,
                                                 seed);

  // --- stochastic Eden growth biased by a noise probability map
  hmap::Array proba = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      {4.f, 4.f},
                                      seed,
                                      8);
  hmap::remap(proba, 0.05f, 0.95f);

  hmap::Array z3 = hmap::elevation_from_contours(shape,
                                                 contours,
                                                 elevations,
                                                 &proba,
                                                 1.f,
                                                 seed);

  // --- deterministic propagation with the same probability map
  hmap::Array z4 = hmap::elevation_from_contours(shape,
                                                 contours,
                                                 elevations,
                                                 &proba,
                                                 0.f);

  // --- evaluation: elevation error along the input contours
  for (const hmap::Array *pz : {&z1, &z2, &z3, &z4})
  {
    float err_max = 0.f;
    float err_mean = 0.f;
    int   n = 0;
    for (size_t k = 0; k < contours.size(); ++k)
    {
      std::vector<float> v = hmap::interpolate_values_from_array(contours[k],
                                                                 *pz,
                                                                 bbox);
      for (float vi : v)
      {
        float e = std::abs(vi - elevations[k]);
        err_max = std::max(err_max, e);
        err_mean += e;
        ++n;
      }
    }
    std::cout << "contour elevation error: max " << err_max << ", mean "
              << err_mean / n << "\n";
  }

  hmap::export_banner_png("ex_elevation_from_contours.png",
                          {z_contours, z1, z2, z3, z4},
                          hmap::Cmap::TERRAIN,
                          true);
}
