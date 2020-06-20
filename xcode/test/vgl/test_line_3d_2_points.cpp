// Some tests for vgl_sphere
// Ian Scott, Aug 2005.
#include <iostream>
#include <limits>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_line_3d_2_points.h>

#include <gtest/gtest.h>

TEST(plane_3d_2_points, intersection)
{
  vgl_plane_3d<double> pl1(vgl_vector_3d<double>(10,10,10), vgl_point_3d<double>(10,0,-10));

  EXPECT_EQ(pl1.d(), 0.0)<<"O is a point on pl1 \n";

  vgl_line_3d_2_points<double> l1(vgl_point_3d<double>(1,4,1),vgl_point_3d<double>(-1,-4,-1));
  EXPECT_EQ(collinear(l1, vgl_point_3d<double>(0,0,0)), true)<<"O is a point on l1\n";

  vgl_line_3d_2_points<double> l2(vgl_point_3d<double>(0,0,0),vgl_point_3d<double>(10,0,-10));
  EXPECT_EQ(collinear(l2, vgl_point_3d<double>(0,0,0)), true)<<"O is a point on l2\n";
  EXPECT_EQ(collinear(l2, vgl_point_3d<double>(10,0,-10)), true)<<"plane_pt is a point on l2\n";

  vgl_line_3d_2_points<double> l3(vgl_point_3d<double>(0,10,0),vgl_point_3d<double>(10,10,-10));
  EXPECT_EQ(collinear(l3, vgl_point_3d<double>(0,0,0)), false)<<"O is not on l3\n";
  EXPECT_EQ(dot_product(pl1.normal(), l3.direction()), 0.0)<<"plane_norm is perpendicular to l3 direction)\n";
}


TEST(plane_3d_2_points, direction_vector)
{
  vgl_point_3d<double> p1(0,0,0);
  vgl_point_3d<double> p2(1,2,3);
  vgl_line_3d_2_points<double> l1(p1, p2);
  vgl_vector_3d<double> u = p2 - p1;
  EXPECT_EQ(u, l1.direction())<<"Direction vector 1\n";
}


TEST(plane_3d_2_points, parametric_point)
{
  vgl_point_3d<double> p1(0, 0, 0);
  vgl_point_3d<double> p2(1, 2, 4);
  vgl_point_3d<double> p3(0.5, 1.0, 2.0);
  vgl_line_3d_2_points<double> l1(p1, p2);
  EXPECT_EQ(l1.point_t(0.0), p1)<<"Parametric point: t=0.0\n";
  EXPECT_EQ(l1.point_t(1.0), p2)<<"Parametric point: t=1.0\n";
  EXPECT_EQ(l1.point_t(0.5), p3)<<"Parametric point: t=0.5\n";
  EXPECT_EQ(l1.point_t(-1.0), vgl_point_3d<double>(-1,-2,-4))<<"Parametric point: t=-1.0\n";
  EXPECT_EQ(l1.point_t(2.0), vgl_point_3d<double>(2,4,8))<<"Parametric point: t=2.0\n";
}



