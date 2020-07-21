#include <iostream>
#include <vgl/vgl_homg_point_2d.h>
#include <vpgl/vpgl_calibration_matrix.h>

#include <gtest/gtest.h>

TEST(vpgl_calibration_matrix, calibration_matrix)
{
  double focal_length = 3;
  vgl_homg_point_2d<double> principal_point(20, 30, 2);
  double x_scale = 2;
  double y_scale = 2;
  double skew = 0;
  vpgl_calibration_matrix<double> K1(focal_length, principal_point, x_scale, y_scale, skew);

  // Test equality of constructors.
  vnl_matrix_fixed<double, 3, 3> M(0.0);
  double scale_factor = -100;
  M(0, 0) = scale_factor * x_scale * focal_length;
  M(1, 1) = scale_factor * y_scale * focal_length;
  M(2, 2) = scale_factor;
  M(0, 2) = scale_factor * 10;
  M(1, 2) = scale_factor * 15;
  vpgl_calibration_matrix<double> K1b(M);

  ASSERT_NEAR(K1.get_matrix() == K1b.get_matrix(), true, 1e-06)<<"test equality of constructors 1\n";
  ASSERT_NEAR(
            K1.focal_length() * K1.x_scale() == K1b.focal_length() * K1b.x_scale() && K1.skew() == K1b.skew(),
            true,
            1e-06)<<"test equality of constructors 2\n";

  // Test the focal length setter.
  focal_length = 5;
  vpgl_calibration_matrix<double> K2(focal_length, principal_point, x_scale, y_scale, skew);
  K1.set_focal_length(focal_length);
  ASSERT_NEAR(K1.get_matrix() == K2.get_matrix(), true, 1e-06)<<"test focal length setter\n";

  // Test the skew setter.
  skew = 2;
  vpgl_calibration_matrix<double> K3(focal_length, principal_point, x_scale, y_scale, skew);
  K1.set_skew(skew);
  ASSERT_NEAR(K1.get_matrix() == K3.get_matrix(), true, 1e-06)<<"test skew setter\n";

  // Test the x scale setter.
  x_scale = 6;
  vpgl_calibration_matrix<double> K4(focal_length, principal_point, x_scale, y_scale, skew);
  K1.set_x_scale(x_scale);
  ASSERT_NEAR(K1.get_matrix() == K4.get_matrix(), true, 1e-06)<<"test x_scale setter\n";

  // Test the principal point setter.
  principal_point.set(17, 100);
  vpgl_calibration_matrix<double> K5(focal_length, principal_point, x_scale, y_scale, skew);
  K1.set_principal_point(principal_point);
  ASSERT_NEAR(K1.get_matrix() == K5.get_matrix(), true, 1e-06)<<"test principal point setter\n";

  // Test the y scale setter.
  y_scale = 6;
  vpgl_calibration_matrix<double> K6(focal_length, principal_point, x_scale, y_scale, skew);
  K1.set_y_scale(y_scale);
  ASSERT_NEAR(K1.get_matrix() == K6.get_matrix(), true, 1e-06)<<"test y_scale setter\n";
}


