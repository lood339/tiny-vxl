// @author fsm
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

#include <vnl/vnl_least_squares_function.h>
#include <vnl/algo/vnl_levenberg_marquardt.h>

#include <gtest/gtest.h>

using vnl_double_2 = vnl_vector_fixed<double, 2>;

class vnl_rosenbrock : public vnl_least_squares_function
{
public:
  vnl_rosenbrock(bool with_grad)
    : vnl_least_squares_function(2, 2, with_grad ? use_gradient : no_gradient)
  {}
    

    void
    f(vnl_vector<double> const & x, vnl_vector<double> & y) const override
    {
        EXPECT_EQ(x.size(), 2);
        EXPECT_EQ(y.size(), 2);
        y[0] = 10 * (x[1] - x[0] * x[0]);
        y[1] = 1 - x[0];
    }

    void
    gradf(vnl_vector<double> const & x, vnl_matrix<double> & J) const override
    {
        EXPECT_EQ(x.size(), 2);
        EXPECT_EQ(J.rows() == 2 && J.cols() == 2, true);
        J[0][0] = -20 * x[0];
        J[0][1] = 10;
        J[1][0] = -1;
        J[1][1] = 0;
    }
};

class linear_est : public vnl_least_squares_function
{
public:
    linear_est(vnl_matrix<double> const &A, vnl_vector<double> b, bool with_grad)
    : vnl_least_squares_function(A.cols(), A.rows(),
                                 with_grad ? use_gradient : no_gradient)
    , A_(A)
    , b_(b)
  {
      assert(A_.rows() == b_.size());
  }
    
    void set_data(const vnl_matrix<double> &A, const vnl_vector<double>& b)
    {
        A_ = A;
        b_ = b;
        assert(A_.rows() == b_.size());
    }
    
    void
    f(vnl_vector<double> const & x, vnl_vector<double> & y) const override
    {
        y = A_ * x - b_;
    }
    
    void
    gradf(vnl_vector<double> const & x, vnl_matrix<double> & J) const override
    {
        J = A_;
    }

    vnl_matrix<double> A_;
    vnl_vector<double> b_;
};


static void
do_rosenbrock_test(bool with_grad)
{
    vnl_rosenbrock f(with_grad);

    vnl_double_2 x0(2.7, -1.3);
    std::cout << "x0 = " << x0 << std::endl;

    vnl_levenberg_marquardt<vnl_rosenbrock> lm(f);
    vnl_vector<double> x1 = x0.as_vector();
    
    if (f.has_gradient())
        lm.minimize_using_gradient(x1);
    else
        lm.minimize_without_gradient(x1);
    
    lm.diagnose_outcome(std::cout);
    std::cout << "x1 = " << x1 << std::endl;

    double err = std::abs(x1[0] - 1) + std::abs(x1[1] - 1);
    std::cout << "err = " << err << std::endl;
    ASSERT_NEAR(err, 0.0, 1e-10)<<"converged to (1, 1)\n";
}


static void
do_linear_test(bool with_grad)
{
    
  vnl_matrix<double> A(6, 2, 1.0);
  vnl_vector<double> b(6);

  A(0, 1) = 10;
  A(1, 1) = 15;
  A(2, 1) = 5.1;
  A(3, 1) = 20.2;
  A(4, 1) = -0.3;
  A(5, 1) = 25;

  b(0) = 10;
  b(1) = 15.5;
  b(2) = 4.5;
  b(3) = 21;
  b(4) = 1;
  b(5) = 24.3;

    
  linear_est f(A, b, with_grad);
  vnl_levenberg_marquardt<linear_est> lm(f);
    
  vnl_vector<double> x(2, -1000.0); // init can be far off
  // since obj function is linear
  // high precision can be achieved
  lm.set_x_tolerance(1e-12);
  lm.set_f_tolerance(1e-12);
  lm.set_g_tolerance(1e-12);

  if (f.has_gradient())
    lm.minimize_using_gradient(x);
  else
    lm.minimize_without_gradient(x);
  lm.diagnose_outcome(std::cout);
    

    vnl_vector<double> true_x(2);
    true_x[0] = 0.595607200429874;
    true_x[1] = 0.969684757298943;
    std::cout << "true      x = \n" << true_x << std::endl;
    std::cout << "estiamted x = \n" << x << std::endl;
    ASSERT_NEAR((true_x - x).two_norm(), 0, 1e-6)<<"converged to true estimate\n";
    
    /*
    // now check (inverse of) covariance approximation
    vnl_matrix<double> true_cov(2, 2);
    true_cov(0, 0) = 6;
    true_cov(1, 0) = 75;
    true_cov(0, 1) = 75;
    true_cov(1, 1) = 1384.14;

    vnl_matrix<double> covar = lm.get_JtJ();
    std::cout << "true      cov(x) =\n" << true_cov << std::endl;
    std::cout << "Estiamted cov(x) =\n" << covar << std::endl;
    ASSERT_NEAR((true_cov - covar).array_two_norm(), 0, 1e-5)<<"covariance approximation\n";
     */
}


TEST(levenberg_marquardt, rosenbrock)
{
    do_rosenbrock_test(false);
    do_rosenbrock_test(true);
}

TEST(levenberg_marquardt, linear_test)
{
    do_linear_test(true);
    do_linear_test(false);
}




