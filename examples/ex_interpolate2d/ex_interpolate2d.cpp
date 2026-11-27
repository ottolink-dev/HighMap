#include "highmap.hpp"
#include "highmap/internal/timer.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  // shape = {1024, 1024};
  int seed = 0;

  hmap::BBox2D bbox = {-1.f, 0.f, 0.5f, 1.5f};
  hmap::UvRect uv_rect = {-1.f, 0.f, 0.5f, 1.5f};

  hmap::Cloud cloud = hmap::Cloud(10, seed, bbox);
  cloud.snap_points_to_bounding_box(bbox);

  std::vector<float> x = cloud.get_x();
  std::vector<float> y = cloud.get_y();
  std::vector<float> values = cloud.get_values();

  // reference, pointwise values
  hmap::Array z0 = hmap::Array(shape);
  cloud.to_array(z0, bbox);

  // noise
  hmap::Array dx = 0.0f * hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                          shape,
                                          {4.f, 4.f},
                                          seed,
                                          8,
                                          0.f);

  // nearest
  hmap::Timer::Start("ITP2D_NEAREST");
  hmap::Array z1 = hmap::interpolate2d(
      shape,
      x,
      y,
      values,
      hmap::InterpolationMethod2D::ITP2D_NEAREST,
      &dx,
      &dx,
      uv_rect);
  hmap::Timer::Stop("ITP2D_NEAREST");

  hmap::Timer::Start("ITP2D_DELAUNAY");
  hmap::Array z2 = hmap::interpolate2d(
      shape,
      x,
      y,
      values,
      hmap::InterpolationMethod2D::ITP2D_DELAUNAY,
      &dx,
      &dx,
      uv_rect);
  hmap::Timer::Stop("ITP2D_DELAUNAY");

  hmap::Timer::Start("ITP2D_IDW");
  hmap::Array z3 = hmap::interpolate2d(shape,
                                       x,
                                       y,
                                       values,
                                       hmap::InterpolationMethod2D::ITP2D_IDW,
                                       &dx,
                                       &dx,
                                       uv_rect);
  hmap::Timer::Stop("ITP2D_IDW");

  hmap::Timer::Start("ITP2D_GAUSSIAN");
  hmap::Array z4 = hmap::interpolate2d(
      shape,
      x,
      y,
      values,
      hmap::InterpolationMethod2D::ITP2D_GAUSSIAN,
      &dx,
      &dx,
      uv_rect);
  hmap::Timer::Stop("ITP2D_GAUSSIAN");

  hmap::Timer::Start("ITP2D_NNI");
  hmap::Array z5 = hmap::interpolate2d(shape,
                                       x,
                                       y,
                                       values,
                                       hmap::InterpolationMethod2D::ITP2D_NNI,
                                       &dx,
                                       &dx,
                                       uv_rect);
  hmap::Timer::Stop("ITP2D_NNI");

  hmap::Timer::Start("ITP2D_CT");
  hmap::Array z6 = hmap::interpolate2d_delaunay_gradient(
      shape,
      x,
      y,
      values,
      // hmap::InterpolationMethod2D::ITP2D_NNI,
      nullptr,
      nullptr,
      uv_rect);
  hmap::Timer::Stop("ITP2D_CT");

  hmap::export_banner_png("ex_interpolate2d.png",
                          {z0, z1, z2, z3, z4, z5, z6},
                          hmap::Cmap::INFERNO);
}
