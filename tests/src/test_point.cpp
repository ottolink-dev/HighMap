#include <glm/glm.hpp>

#include "highmap/geometry/point.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------

static constexpr float eps = 1e-5f;

bool float_eq(float a, float b, float tol = eps)
{
  return std::abs(a - b) < tol;
}

// ------------------------------------------------------------
// Constructors
// ------------------------------------------------------------

TEST(PointTest, DefaultConstructor)
{
  Point p;

  EXPECT_TRUE(float_eq(p.x, 0.f));
  EXPECT_TRUE(float_eq(p.y, 0.f));
  EXPECT_TRUE(float_eq(p.v, 0.f));
}

TEST(PointTest, ParameterizedConstructor)
{
  Point p(1.f, 2.f, 3.f);

  EXPECT_TRUE(float_eq(p.x, 1.f));
  EXPECT_TRUE(float_eq(p.y, 2.f));
  EXPECT_TRUE(float_eq(p.v, 3.f));
}

TEST(PointTest, GlmVec2ConstructorAndConversion)
{
  glm::vec2 v(1.5f, 2.5f);
  Point     p(v, 4.f);

  EXPECT_TRUE(float_eq(p.x, 1.5f));
  EXPECT_TRUE(float_eq(p.y, 2.5f));
  EXPECT_TRUE(float_eq(p.v, 4.f));

  glm::vec2 out = p.to_vec2();
  EXPECT_TRUE(float_eq(out.x, 1.5f));
  EXPECT_TRUE(float_eq(out.y, 2.5f));
}

TEST(PointTest, GlmVec3ConstructorAndConversion)
{
  glm::vec3 v(1.5f, 2.5f, 3.5f);
  Point     p(v);

  EXPECT_TRUE(float_eq(p.x, 1.5f));
  EXPECT_TRUE(float_eq(p.y, 2.5f));
  EXPECT_TRUE(float_eq(p.v, 3.5f));

  glm::vec3 out = p.to_vec3();
  EXPECT_TRUE(float_eq(out.x, 1.5f));
  EXPECT_TRUE(float_eq(out.y, 2.5f));
  EXPECT_TRUE(float_eq(out.z, 3.5f));
}

// ------------------------------------------------------------
// Operators
// ------------------------------------------------------------

TEST(PointTest, EqualityOperator)
{
  Point a(1, 2, 3);
  Point b(1, 2, 3);

  EXPECT_TRUE(a == b);
  EXPECT_FALSE(a != b);
}

TEST(PointTest, InequalityOperator)
{
  Point a(1, 2, 3);
  Point b(1, 2, 4);

  EXPECT_TRUE(a != b);
}

TEST(PointTest, AdditionOperator)
{
  Point a(1, 2, 3);
  Point b(4, 5, 6);

  Point c = a + b;

  EXPECT_TRUE(float_eq(c.x, 5));
  EXPECT_TRUE(float_eq(c.y, 7));
  EXPECT_TRUE(float_eq(c.v, 9));
}

TEST(PointTest, SubtractionOperator)
{
  Point a(5, 7, 9);
  Point b(1, 2, 3);

  Point c = a - b;

  EXPECT_TRUE(float_eq(c.x, 4));
  EXPECT_TRUE(float_eq(c.y, 5));
  EXPECT_TRUE(float_eq(c.v, 6));
}

TEST(PointTest, ScalarMultiplicationRight)
{
  Point a(1, 2, 3);

  Point c = a * 2.f;

  EXPECT_TRUE(float_eq(c.x, 2));
  EXPECT_TRUE(float_eq(c.y, 4));
  EXPECT_TRUE(float_eq(c.v, 6));
}

TEST(PointTest, ScalarMultiplicationLeft)
{
  Point a(1, 2, 3);

  Point c = 2.f * a;

  EXPECT_TRUE(float_eq(c.x, 2));
  EXPECT_TRUE(float_eq(c.y, 4));
  EXPECT_TRUE(float_eq(c.v, 6));
}

TEST(PointTest, ScalarDivision)
{
  Point a(2, 4, 6);

  Point c = a / 2.f;

  EXPECT_TRUE(float_eq(c.x, 1));
  EXPECT_TRUE(float_eq(c.y, 2));
  EXPECT_TRUE(float_eq(c.v, 3));
}

TEST(PointTest, CompoundOperators)
{
  Point a(1, 2, 3);
  Point b(4, 5, 6);

  a += b;
  EXPECT_TRUE(float_eq(a.x, 5));
  EXPECT_TRUE(float_eq(a.y, 7));
  EXPECT_TRUE(float_eq(a.v, 9));

  a -= b;
  EXPECT_TRUE(float_eq(a.x, 1));
  EXPECT_TRUE(float_eq(a.y, 2));
  EXPECT_TRUE(float_eq(a.v, 3));

  a *= 2.f;
  EXPECT_TRUE(float_eq(a.x, 2));
  EXPECT_TRUE(float_eq(a.y, 4));
  EXPECT_TRUE(float_eq(a.v, 6));

  a /= 2.f;
  EXPECT_TRUE(float_eq(a.x, 1));
  EXPECT_TRUE(float_eq(a.y, 2));
  EXPECT_TRUE(float_eq(a.v, 3));
}

// ------------------------------------------------------------
// Geometry
// ------------------------------------------------------------

TEST(PointTest, DistanceBasic)
{
  Point a(0, 0, 0);
  Point b(3, 4, 0);

  EXPECT_TRUE(float_eq(distance(a, b), 5.f));
  EXPECT_TRUE(float_eq(distance_squared(a, b), 25.f));
}

TEST(PointTest, CrossProductOrientation)
{
  Point p0(0, 0, 0);
  Point p1(1, 0, 0);
  Point p2(0, 1, 0);

  float cross = cross_product(p0, p1, p2);

  EXPECT_GT(cross, 0.f); // CCW
}

TEST(PointTest, CrossProductColinear)
{
  Point p0(0, 0, 0);
  Point p1(1, 1, 0);
  Point p2(2, 2, 0);

  EXPECT_TRUE(float_eq(cross_product(p0, p1, p2), 0.f));
}

TEST(PointTest, AngleBetweenPoints)
{
  Point a(0, 0, 0);
  Point b(1, 0, 0);

  float ang = angle(a, b);

  EXPECT_TRUE(float_eq(ang, 0.f));
}

TEST(PointTest, AngleThreePointsRightAngle)
{
  Point p0(0, 0, 0);
  Point p1(1, 0, 0);
  Point p2(0, 1, 0);

  float ang = angle(p0, p1, p2);

  EXPECT_NEAR(std::abs(ang), M_PI / 2.f, 1e-4f);
}

// ------------------------------------------------------------
// Bounding Box
// ------------------------------------------------------------

TEST(PointTest, PointInsideBoundingBox)
{
  Point p(0.5f, 0.5f, 0.f);

  glm::vec4 bbox{0, 1, 0, 1};

  EXPECT_TRUE(is_point_within_bounding_box(p, bbox));
}

TEST(PointTest, PointOutsideBoundingBox)
{
  Point p(2.f, 0.5f, 0.f);

  glm::vec4 bbox{0, 1, 0, 1};

  EXPECT_FALSE(is_point_within_bounding_box(p, bbox));
}

// ------------------------------------------------------------
// Interpolation
// ------------------------------------------------------------

TEST(PointTest, LerpEndpoints)
{
  Point a(0, 0, 0);
  Point b(10, 10, 10);

  Point p0 = lerp(a, b, 0.f);
  Point p1 = lerp(a, b, 1.f);

  EXPECT_TRUE(p0 == a);
  EXPECT_TRUE(p1 == b);
}

TEST(PointTest, LerpMidpoint)
{
  Point a(0, 0, 0);
  Point b(10, 10, 10);

  Point m = lerp(a, b, 0.5f);

  EXPECT_TRUE(float_eq(m.x, 5.f));
  EXPECT_TRUE(float_eq(m.y, 5.f));
  EXPECT_TRUE(float_eq(m.v, 5.f));
}

TEST(PointTest, BezierEndpoints)
{
  Point p0(0, 0, 0);
  Point p1(1, 0, 0);
  Point p2(1, 1, 0);
  Point p3(0, 1, 0);

  EXPECT_TRUE(interp_bezier(p0, p1, p2, p3, 0.f) == p0);
  EXPECT_TRUE(interp_bezier(p0, p1, p2, p3, 1.f) == p3);
}

TEST(PointTest, CatmullRomEndpoints)
{
  Point p0(0, 0, 0);
  Point p1(1, 1, 1);
  Point p2(2, 2, 2);
  Point p3(3, 3, 3);

  EXPECT_TRUE(interp_catmullrom(p0, p1, p2, p3, 0.f) == p1);
  EXPECT_TRUE(interp_catmullrom(p0, p1, p2, p3, 1.f) == p2);
}

TEST(PointTest, DecasteljauEmptyAndEndpoints)
{
  std::vector<Point> empty_pts;
  EXPECT_TRUE(interp_decasteljau(empty_pts, 0.5f) == Point(0, 0, 0));

  std::vector<Point> single_pt = {Point(3, 4, 5)};
  EXPECT_TRUE(interp_decasteljau(single_pt, 0.5f) == Point(3, 4, 5));

  std::vector<Point> pts = {Point(0, 0, 0),
                            Point(1, 0, 0),
                            Point(1, 1, 0),
                            Point(0, 1, 0)};
  EXPECT_TRUE(interp_decasteljau(pts, 0.f) == pts.front());
  EXPECT_TRUE(interp_decasteljau(pts, 1.f) == pts.back());
  EXPECT_TRUE(interp_decasteljau(pts, 0.5f) ==
              interp_bezier(pts[0], pts[1], pts[2], pts[3], 0.5f));
}

// ------------------------------------------------------------
// Segment Intersection
// ------------------------------------------------------------

TEST(PointTest, SegmentIntersectionExists)
{
  Point p1(0, 0, 0), p2(1, 1, 0);
  Point q1(0, 1, 0), q2(1, 0, 0);

  auto result = segment_intersection(p1, p2, q1, q2);

  ASSERT_TRUE(result.has_value());

  EXPECT_NEAR(result->x, 0.5f, eps);
  EXPECT_NEAR(result->y, 0.5f, eps);
}

TEST(PointTest, SegmentIntersectionNone)
{
  Point p1(0, 0, 0), p2(1, 0, 0);
  Point q1(0, 1, 0), q2(1, 1, 0);

  auto result = segment_intersection(p1, p2, q1, q2);

  EXPECT_FALSE(result.has_value());
}

// ------------------------------------------------------------
// Triangle Area
// ------------------------------------------------------------

TEST(PointTest, TriangleAreaBasic)
{
  Point a(0, 0, 0);
  Point b(1, 0, 0);
  Point c(0, 1, 0);

  EXPECT_TRUE(float_eq(triangle_area(a, b, c), 0.5f));
}

TEST(PointTest, TriangleAreaColinear)
{
  Point a(0, 0, 0);
  Point b(1, 1, 0);
  Point c(2, 2, 0);

  EXPECT_TRUE(float_eq(triangle_area(a, b, c), 0.f));
}

// ------------------------------------------------------------
// Sorting
// ------------------------------------------------------------

TEST(PointTest, SortPointsLexicographic)
{
  std::vector<Point> pts = {{2, 1, 0}, {1, 2, 0}, {1, 1, 0}};

  sort_points(pts);

  EXPECT_TRUE(pts[0] == Point(1, 1, 0));
  EXPECT_TRUE(pts[1] == Point(1, 2, 0));
  EXPECT_TRUE(pts[2] == Point(2, 1, 0));
}

// ------------------------------------------------------------
// Scale
// ------------------------------------------------------------

TEST(PointTest, ScaleUniformDefaultCenter)
{
  Point p(0.f, 1.f, 5.f);
  Point scaled = scale(p, 0.8f);

  EXPECT_NEAR(scaled.x, 0.1f, eps);
  EXPECT_NEAR(scaled.y, 0.9f, eps);
  EXPECT_NEAR(scaled.v, 5.f, eps);
}

TEST(PointTest, ScaleNonUniformCustomCenter)
{
  Point p(1.f, 2.f, 3.f);
  Point scaled = scale(p, glm::vec2(2.f, 0.5f), glm::vec2(0.f, 0.f));

  EXPECT_NEAR(scaled.x, 2.f, eps);
  EXPECT_NEAR(scaled.y, 1.f, eps);
  EXPECT_NEAR(scaled.v, 3.f, eps);
}
