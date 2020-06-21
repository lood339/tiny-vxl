// Some tests for vgl_line_segment_3d
// Kevin de Souza, Aug 2005.


#include <vgl/vgl_line_segment_3d.h>

#include <gtest/gtest.h>


TEST(line_segment_3d, direction) {
    vgl_point_3d<double> p1(0,0,0);
    vgl_point_3d<double> p2(1,2,3);
    vgl_line_segment_3d<double> l1(p1, p2);
    vgl_vector_3d<double> u = p2 - p1;
    EXPECT_EQ(u, l1.direction());
}

TEST(line_segment_3d, parametric_point) {
    vgl_point_3d<double> p1(0, 0, 0);
    vgl_point_3d<double> p2(1, 2, 4);
    vgl_point_3d<double> p3(0.5, 1.0, 2.0);
    vgl_line_segment_3d<double> l1(p1, p2);
    EXPECT_EQ(l1.point_t(0.0), p1)<<"Parametric point: t=0.0\n";
    EXPECT_EQ(l1.point_t(1.0), p2)<<"Parametric point: t=1.0\n";
    EXPECT_EQ(l1.point_t(0.5), p3)<<"Parametric point: t=0.5\n";
    EXPECT_EQ(l1.point_t(-1.0), vgl_point_3d<double>(-1,-2,-4))<<"Parametric point: t=-1.0\n";
    EXPECT_EQ(l1.point_t(2.0), vgl_point_3d<double>(2,4,8))<<"Parametric point: t=2.0\n";
}





