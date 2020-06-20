// Some tests for vgl_sphere_3d
// Ian Scott, Aug 2005.
#include <iostream>
#include <vgl/vgl_sphere_3d.h>
#include <vgl/vgl_line_3d_2_points.h>
#include <gtest/gtest.h>

TEST(sphere_3d, sphere)
{
    // Default sphere
    vgl_sphere_3d<double> s;
    // Unit sphere, centered at the origin
    vgl_sphere_3d<double> u(0, 0, 0, 1.0);
    
    EXPECT_EQ(s.is_empty(), true)<<"default is empty\n";
    EXPECT_EQ(u.is_empty(), false);
    
    EXPECT_EQ(s.contains( vgl_point_3d<double>(0,0,0) ), false);
    EXPECT_EQ(u.contains( vgl_point_3d<double>(0,0,0) ), true)<<"origin is inside unit sphere\n";
    EXPECT_EQ(u.contains( vgl_point_3d<double>(1,0,0) ), true)<<"(1,0,0) is inside unit sphere\n";
    EXPECT_EQ(u.contains( vgl_point_3d<double>(1,1,0) ), false)<<"(1,1,0) is outside unit sphere\n";
    
    // l1 is the X-axis
    vgl_line_3d_2_points<double> l1(vgl_point_3d<double>(-2,0,0),vgl_point_3d<double>(2,0,0));
    vgl_point_3d<double> p1, p2;
    EXPECT_EQ(s.clip(l1, p1, p2), false)<<"clip x-axis to empty sphere\n";
    EXPECT_EQ( u.clip(l1, p1, p2), true)<<"clip x-axis to unit sphere\n";
    EXPECT_EQ(p1, vgl_point_3d<double>(-1,0,0))<<"Intersection point 1\n";
    EXPECT_EQ(p2, vgl_point_3d<double>(1,0,0))<<"Intersection point 2\n";
    
    // l2 is the line (y=1,z=0) parallel to the X axis, touching the unit sphere in (0,1,0)
    vgl_line_3d_2_points<double> l2(vgl_point_3d<double>(-2,1,0),vgl_point_3d<double>(2,1,0));
    EXPECT_EQ(u.clip(l2, p1, p2), true)<<"clip (y=1,z=0) to unit sphere\n";
    EXPECT_EQ(p1, vgl_point_3d<double>(0,1,0))<<"Intersection point 1\n";
    EXPECT_EQ(p2, vgl_point_3d<double>(0,1,0))<<"Intersection point 2\n";
    
    // Test basic i/o
    std::cout << u << std::endl;
}

