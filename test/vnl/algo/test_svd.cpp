// This is core/vnl/algo/tests/test_svd.cxx
#include <iostream>
#include <complex>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_random.h>
#include <vnl/algo/vnl_svd.h>

#include <gtest/gtest.h>

template <class T, class S>
static void
test_hilbert(T /*dummy*/, char const * type, S residual)
{
    std::cout << "----- Testing svd<" << type << ">(Hilbert_3x3) -----" << std::endl;
    using abs_t = typename vnl_numeric_traits<T>::abs_t;
    // Test inversion and recomposition of 5x5 hilbert matrix
    vnl_matrix<T> H(5, 5);
    for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      H(i, j) = T(1) / T(abs_t(i + j + 1)); // sic, because i,j are zero based

    std::cout << "H = <" << type << ">[ " << H << "]\n";

    vnl_svd<T> svd(H);
    std::cout << "rcond(H) = " << svd.well_condition() << std::endl;

    vnl_matrix<T> Hinv = svd.inverse();

    vnl_matrix<T> X = Hinv * H;

    std::cout << "H*inv(H) = " << X << std::endl;

    vnl_matrix<T> I(5, 5);
    I = 0.0;
    I.fill_diagonal(1.0);

    vnl_matrix<T> res = X - I;
    ASSERT_NEAR(res.fro_norm(), 0, residual)<<"Hilbert recomposition residual\n";
    
}


TEST(nvl_svd, test_hilbert_double)
{
    test_hilbert(double(), "double", 1.1e-10);
}

TEST(nvl_svd, test_hilbert_float)
{
    test_hilbert(float(), "float", float(0.025));
}

TEST(nvl_svd, test_hilbert_complex_double)
{
    test_hilbert(std::complex<double>(), "std::complex<double>", double(4.4e-10));
}

TEST(vnl_svd, test_hilbert_complex_float)
{
    test_hilbert(std::complex<float>(), "std::complex<float>", float(0.04));
}

//: Test recovery of parameters of least-squares parabola fit.
TEST(vnl_svd, test_least_squares_fit)
{
    std::cout << "----- Testing svd on a Least Squares problem -----" << std::endl;
    double a = 0.15;
    double b = 1.2;
    double c = 3.1;
    
    // Generate parabola design matrix
    vnl_matrix<double> D(100, 3);
    for (int n = 0; n < 100; ++n)
    {
        double x = n;
        D(n, 0) = x * x;
        D(n, 1) = x;
        D(n, 2) = 1.0;
    }
    
    // Generate Y vector
    vnl_vector<double> y(100);
    for (int n = 0; n < 100; ++n)
    {
        double x = n;
        double fx = a * x * x + b * x + c;
        // Add sawtooth "noise"
        y(n) = fx + (n % 4 - 2) / 10.0;
    }
    std::cout << "y = [" << y << "]\n";
    
    // Extract vnl_svd<double>
    vnl_svd<double> svd(D);
    
    // Solve for parameters
    vnl_vector<double> A = svd.solve(y);
    std::cout << "A = " << A << '\n';
    
    vnl_vector<double> T(3);
    T[0] = a;
    T[1] = b;
    T[2] = c;
   
    ASSERT_NEAR((A - T).squared_magnitude(), 0, 0.005)<<"Least squares residual\n";
}



//: Test nullspace extraction of rank=2 3x4 matrix.

TEST(vnl_svd, test_pmatrix)
{
    double pdata[] = {
        2, 0, 0, 0,
        3, 10, 5, 5,
        5, 12, 6, 6,
    };
    vnl_matrix<double> P(pdata, 3, 4);
    vnl_svd<double> svd(P, 1e-8);

    vnl_matrix<double> res = svd.recompose() - P;
    ASSERT_NEAR(res.fro_norm(), 0, 1e-12)<<"PMatrix recomposition residual\n";
    std::cout << " Inv = " << svd.inverse() << std::endl;

    EXPECT_EQ(svd.singularities(), 2)<<"singularities = 2\n";
    EXPECT_EQ(svd.rank(), 2)<<"rank = 2\n";
    
    vnl_matrix<double> N = svd.nullspace();
    EXPECT_EQ(N.columns(), 2)<<"nullspace dimension\n";
    std::cout << "null(P) =\n" << N << std::endl;
    
    vnl_matrix<double> PN = P * N;
    std::cout << "P * null(P) =\n" << PN << std::endl;
    ASSERT_NEAR(PN.fro_norm(), 0, 1e-12)<<"P nullspace residual\n";

    vnl_vector<double> n = svd.nullvector();
    ASSERT_NEAR((P * n).magnitude(), 0, 1e-12)<<"P nullvector residual\n";

    vnl_vector<double> l = svd.left_nullvector();
    std::cout << "left_nullvector(P) = " << l << std::endl;
    ASSERT_NEAR((l * P).magnitude(), 0, 1e-12)<<"P left nullvector residual\n";
}

TEST(vnl_svd, test_I)
{
  double Idata[] = {
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
  };
  vnl_matrix<double> P(3, 4, 12, Idata);
  vnl_svd<double> svd(P);
  std::cout << svd;

  vnl_vector_fixed<double, 4> w_expected(1, 1, 1, 0);
  ASSERT_NEAR(vnl_vector_ssd(w_expected, svd.W().diagonal()), 0, 1e-16)<<"Singular values\n";
}

template <typename T>
static
void test_util_fill_random(T * b, T * e, vnl_random & rng)
{
    for (T * p = b; p < e; ++p)
        *p = (T)rng.drand64(-1.0, +1.0);
}

template <typename T>
static
void test_util_fill_random(std::complex<T> * b, std::complex<T> * e, vnl_random & rng)
{
    for (std::complex<T> * p = b; p < e; ++p)
        *p = std::complex<T>((T)rng.drand64(-1.0, +1.0), (T)rng.drand64(-1.0, +1.0));
}

template <class T>
void
test_svd_recomposition(char const * type, double maxres, T *, vnl_random & rng)
{
  // Test inversion of 5x5 matrix of T :
  std::cout << "----- Testing vnl_svd<" << type << "> recomposition -----\n";

  vnl_matrix<T> A(5, 5);
  test_util_fill_random(A.begin(), A.end(), rng);

  std::cout << "A = [\n" << A << "]\n";
  vnl_svd<T> svd(A);

  vnl_matrix<T> B = svd.recompose();
  std::cout << "B = [\n" << B << "]\n";

  double residual = (A - B).fro_norm();
  ASSERT_NEAR(residual, 0, maxres)<<"vnl_svd<float> recomposition residual\n";
}

template <class T>
static void
test_nullvector(char const * type, double max_err, T *, vnl_random & rng)
{
  int n = 5;
  vnl_matrix<T> A(n, n + 1);
  test_util_fill_random(A.begin(), A.end(), rng);
  vnl_svd<T> svd(A);
  vnl_vector<T> x = svd.nullvector();
  vnl_vector<T> Ax = A * x;
  //std::cout << __FILE__ ": type = " << type << std::endl;
  //vnl_matlab_print(std::cout, A, "A", vnl_matlab_print_format_long);
  //std::cout << __FILE__ ": || x|| = " << x.two_norm() << std::endl
  //          << __FILE__ ": ||Ax|| = " << Ax.two_norm() << std::endl;
  ASSERT_NEAR(Ax.two_norm(), 0.0, max_err)<<"||Ax||\n";
}

TEST(vnl_svd, recomposition)
{
  vnl_random rng(9667566);
  
  test_svd_recomposition("float", 1e-5, (float *)nullptr, rng);
  test_svd_recomposition("double", 1e-10, (double *)nullptr, rng);
  test_svd_recomposition("std::complex<float>", 1e-5, (std::complex<float> *)nullptr, rng);
  test_svd_recomposition("std::complex<double>", 1e-10, (std::complex<double> *)nullptr, rng);
}

TEST(vnl_svd, nullvector)
{
    vnl_random rng(9667566);
    
    test_nullvector("float", 5e-7, (float *)nullptr, rng);
    test_nullvector("double", 5e-15, (double *)nullptr, rng);
    test_nullvector("std::complex<float>", 5e-7, (std::complex<float> *)nullptr, rng);
    test_nullvector("std::complex<double>", 5e-15, (std::complex<double> *)nullptr, rng);
}

