#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {64, 32};
  shape = {512, 256};

  // --- 1st frame

  hmap::CoordFrame frame1 = hmap::CoordFrame(glm::vec2(10.f, 20.f),
                                             glm::vec2(50.f, 100.f),
                                             30.f);

  // --- second frame

  hmap::CoordFrame frame2 = hmap::CoordFrame(glm::vec2(-20.f, 50.f),
                                             glm::vec2(100.f, 70.f),
                                             -30.f);

  // --- tests

  hmap::BBox2D bbox1 = frame1.compute_bounding_box();
  hmap::BBox2D bbox2 = frame2.compute_bounding_box();

  hmap::BBox2D bboxi = hmap::intersect_bounding_boxes(bbox1, bbox2);

  glm::ivec2 shape2 = {1000, 1000};

  std::vector<float> x = hmap::linspace(-200.f, 200.f, shape2.x);
  std::vector<float> y = hmap::linspace(-200.f, 200.f, shape2.y);

  hmap::Array array(shape2);

  for (int i = 0; i < shape2.x; i++)
    for (int j = 0; j < shape2.y; j++)
    {
      if (frame1.is_point_within(x[i], y[j]))
        array(i, j) = 1.f;
      else if (frame2.is_point_within(x[i], y[j]))
        array(i, j) = 0.75f;
      else if (hmap::is_point_within_bounding_box(x[i], y[j], bboxi))
        array(i, j) = 0.5f;
      else if (hmap::is_point_within_bounding_box(x[i], y[j], bbox1))
        array(i, j) = 0.1f;
      else if (hmap::is_point_within_bounding_box(x[i], y[j], bbox2))
        array(i, j) = 0.05f;
    }

  array.to_png("out.png", hmap::Cmap::JET);
}
