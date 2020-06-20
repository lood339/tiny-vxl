// Some tests for vgl_ray_3d
// J.L. Mundy Sept. 17, 2010

#include <iostream>
#include <vgl/vgl_ray_3d.h>
#include <vgl/vgl_closest_point.h>

#include <gtest/gtest.h>

TEST(ray_3d, constructor)
{
    vgl_vector_3d<double> t(0,0,2);
    vgl_point_3d<double> p(1,2,3);
    vgl_ray_3d<double> ray(p, t);
    vgl_point_3d<double> origin = ray.origin();
    vgl_vector_3d<double> dir = ray.direction();
    ASSERT_NEAR(origin.x()+origin.y() , 3.0, 1e-5)<<"Constructor from point and dir - compare origin\n";
    ASSERT_NEAR(dir.z_ ,1.0, 1e-5)<<"Constructor from point and dir - compare dir\n";
    vgl_point_3d<double> p1(1,2,4);
    vgl_ray_3d<double> ray1(p, p1);
    origin = ray1.origin();
    dir = ray1.direction();
    ASSERT_NEAR(origin.x()+origin.y() , 3.0, 1e-5)<<"Constructor from point-point - compare origin\n";
    ASSERT_NEAR(dir.z_ ,1.0, 1e-5)<<"Constructor from point-point - compare dir\n";
}

TEST(ray_3d, operations)
{
    vgl_vector_3d<double> t(0,0,1);
    vgl_point_3d<double> p(1,2,3), pt(1,2,2), clpt;
    vgl_ray_3d<double> ray(p, t);
    clpt = vgl_closest_point(ray, pt);
    bool con = ray.contains(clpt);
    ASSERT_EQ(con, false)<<"Contains \n";
}

