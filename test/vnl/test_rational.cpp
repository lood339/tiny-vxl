#include <iostream>
#include <iomanip>
#include <complex>
#include <string>
#include <utility>
#include <typeinfo>

#include <vnl/vnl_rational.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_matrix_fixed.h>
//#include <vnl/vnl_rational_traits.h>
#include <vnl/vnl_det.h>

#include <gtest/gtest.h>

inline vnl_rational
vnl_sqrt(vnl_rational x)
{
  return vnl_rational(std::sqrt(double(x)));
}

namespace
{
template <typename T>
void
test_converting_whole_number_to_rational(const T num)
{
    // Convert whole number to vnl_rational:
    const vnl_rational rat = num;

    const auto message =
    "test_converting_value_to_rational<" + std::string(typeid(T).name()) + ">(" + std::to_string(num) + ')';

    using pair_type = std::pair<vnl_rational::int_type, vnl_rational::int_type>;

    EXPECT_EQ(std::make_pair(rat.numerator(), rat.denominator()), pair_type(num, 1))<<message.c_str()<<std::endl;
}

template <typename T>
void
test_converting_floating_point_number_to_rational(T num)
{
    // Convert floating point number to vnl_rational:
    const vnl_rational rat = num;

    const auto message = "test_converting_floating_point_number_to_rational<" + std::string(typeid(T).name()) + ">(" +
                       std::to_string(num) + ')';

    using pair_type = std::pair<vnl_rational::int_type, vnl_rational::int_type>;

    EXPECT_EQ(static_cast<T>(rat.numerator()) / static_cast<T>(rat.denominator()), num)<<message.c_str()<<std::endl;
}

template <typename T>
void
test_converting_decimal_digits()
{
  for (int i{}; i < 10; ++i)
  {
    test_converting_whole_number_to_rational(static_cast<T>(i));
  }
}

template <typename T>
void
test_converting_decimal_digits_and_min_and_max()
{
  test_converting_decimal_digits<T>();

  test_converting_whole_number_to_rational(std::numeric_limits<T>::min());
  test_converting_whole_number_to_rational(std::numeric_limits<T>::max());
}

} // namespace

TEST(vnl_rational, converting_constructors)
{
  test_converting_decimal_digits_and_min_and_max<unsigned char>();
  test_converting_decimal_digits_and_min_and_max<signed char>();
  test_converting_decimal_digits_and_min_and_max<unsigned char>();
  test_converting_decimal_digits_and_min_and_max<short>();
  test_converting_decimal_digits_and_min_and_max<unsigned short>();
  test_converting_decimal_digits_and_min_and_max<int>();
  test_converting_decimal_digits_and_min_and_max<unsigned int>();
  test_converting_decimal_digits_and_min_and_max<long>();
  test_converting_decimal_digits_and_min_and_max<unsigned long>();
  test_converting_decimal_digits_and_min_and_max<long long>();
  test_converting_decimal_digits_and_min_and_max<unsigned long long>();

  test_converting_decimal_digits<float>();
  test_converting_decimal_digits<double>();

  for (float f{}; f <= 1.0f; f += 0.25f)
  {
    test_converting_floating_point_number_to_rational<float>(f);
    test_converting_floating_point_number_to_rational<double>(f);
  }
}


TEST(vnl_rational, operators)
{
  vnl_rational a(-5L), b(7, -1), c, d(3, 7), e(2, 0);
  vnl_rational z_default;
  EXPECT_EQ(z_default == 0L, true);

  vnl_rational z_int(static_cast<int>(0));
  EXPECT_EQ(z_int == 0L, true);
  vnl_rational z_uint(static_cast<unsigned int>(0));
  EXPECT_EQ(z_uint == 0L, true);

  vnl_rational z_short(static_cast<int>(0));
  EXPECT_EQ(z_short == 0L, true);
  vnl_rational z_ushort(static_cast<unsigned int>(0));
  EXPECT_EQ(z_ushort == 0L, true);

  vnl_rational z_long(static_cast<long>(0));
  EXPECT_EQ(z_long == 0L, true);
  vnl_rational z_ulong(static_cast<unsigned long>(0));
  EXPECT_EQ(z_ulong == 0L, true);

  EXPECT_EQ(a == -5L, true);
  EXPECT_EQ(5L == -a, true);
  EXPECT_EQ(b == -7, true);
  EXPECT_EQ(-7 == b, true);
  c = a + b;
  EXPECT_EQ(c, -12L);
  c = a - b;
  EXPECT_EQ(c, 2L);
  c = a * b;
  EXPECT_EQ(c, 35L);
  c = a / b;
  EXPECT_EQ(c, vnl_rational(5, 7));
  c = c % d;
  EXPECT_EQ(c, vnl_rational(2, 7));
  c = a % b;
  EXPECT_EQ(c, -5L);
  c = a % d;
  EXPECT_EQ(c, vnl_rational(-2, 7));
  c = d % a;
  EXPECT_EQ(c, d);
  c = a + 5L;
  EXPECT_EQ(c, 0L);
  c = a - 5L;
  EXPECT_EQ(c, -10L);
  c = a * 5L;
  EXPECT_EQ(c, -25L);
  c = a / 5L;
  EXPECT_EQ(c, -1L);
  c = a % 5L;
  EXPECT_EQ(c, 0L);
  c = 5L + a;
  EXPECT_EQ(c, 0L);
  c = 5L - a;
  EXPECT_EQ(c, 10L);
  c = 5L * a;
  EXPECT_EQ(c, -25L);
  c = 5L / a;
  EXPECT_EQ(c, -1L);
  c = 5L % a;
  EXPECT_EQ(c, 0L);
  c = 5 + a;
  EXPECT_EQ(c, 0L);
  c = 5 - a;
  EXPECT_EQ(c, 10L);
  c = 5 * a;
  EXPECT_EQ(c, -25L);
  c = 5 / a;
  EXPECT_EQ(c, -1L);
  c = 5 % a;
  EXPECT_EQ(c, 0L);
  c = a + 5;
  EXPECT_EQ(c, 0L);
  c = a - 5;
  EXPECT_EQ(c, -10L);
  c = a * 5;
  EXPECT_EQ(c, -25L);
  c = a / 5;
  EXPECT_EQ(c, -1L);
  EXPECT_EQ(a < d, true);
  EXPECT_EQ(a < 1L, true);
  EXPECT_EQ(a < -4.9, true);
  EXPECT_EQ( -b > d, true);
  EXPECT_EQ(b > -8, true);
  EXPECT_EQ(b > -7.1, true);
  EXPECT_EQ(c <= e, true);
  EXPECT_EQ(b >= -7L, true);
  EXPECT_EQ(2L <= e, true);
  EXPECT_EQ(1 >= d, true);
  EXPECT_EQ(truncate(1L + d), 1L);
  EXPECT_EQ(truncate(-d - 1L), -1L);
  EXPECT_EQ(round(1L + d), 1L);
  EXPECT_EQ(round(-d - 1L), -1L);
  EXPECT_EQ(round(1L - d), 1L);
  EXPECT_EQ(round(d - 1L), -1L);
  EXPECT_EQ(floor(1L + d), 1L);
  EXPECT_EQ(floor(-d - 1L), -2L);
  EXPECT_EQ(ceil(1L + d), 2L);
  EXPECT_EQ(ceil(-d - 1L), -1L);
  EXPECT_EQ(vnl_math::abs(d), d);
  EXPECT_EQ(vnl_math::abs(b), -b);
  EXPECT_EQ(vnl_math::squared_magnitude(d), vnl_rational(9, 49));
  a += b;
  a -= b;
  a *= b;
  a /= b;
  a %= b;
  std::cout << std::setprecision(20) << "a=" << a << '=' << (double)a << '\n'
            << "b=" << b << '=' << (double)b << '\n'
            << "c=" << c << '=' << (double)c << '\n'
            << "d=" << d << '=' << (double)d << '\n'
            << "e=" << e << std::endl; // (double)d ==> floating exception
  d = -7;
  d = -7L;
  std::cout << std::endl;
}

TEST(vnl_rational, infinite)
{
  vnl_rational Inf(1, 0);
  ++Inf;
  EXPECT_EQ(Inf.numerator() == 1 && Inf.denominator() == 0, true)<<"Inf+1\n";
  Inf = -Inf;
  EXPECT_EQ(Inf.numerator() == -1 && Inf.denominator() == 0, true)<<"-Inf\n";
  EXPECT_EQ(vnl_math::isfinite(Inf), false);
  EXPECT_EQ(vnl_math::isnan(Inf), false);
}


TEST(vnl_rational, frac)
{
  vnl_rational r(-15, -20), s(1234321L, -1111111L), p;
  EXPECT_EQ(vnl_math::isfinite(r), true);
  EXPECT_EQ(vnl_math::isnan(r), false);
  EXPECT_EQ(r.numerator() == 3 && r.denominator() == 4, true)<<"simplify\n";
  EXPECT_EQ(s.numerator() == -1234321L && s.denominator() == 1111111L, true)<<"sign in numerator\n";
  // All 5-digit numbers below are prime numbers, and small enough so that the multiplications in the constructors do
  // not overflow
  long p1 = 46309L, p2 = 46349L, p3 = 46327L, p4 = 46337L, p5 = 46351L;
  r = vnl_rational(p1 * p2, p3 * p4);
  s = vnl_rational(p3 * p4, p1 * p5);
  p = r * s;
  EXPECT_EQ(p.numerator() == p2 && p.denominator() == p5, true)<<"large multiplication without overflow\n";
  r = vnl_rational(p1 * p2, p3 * p4);
  s = vnl_rational(p1 * p5, p3 * p4);
  p = r * s;
  ASSERT_NEAR(p, double(r) * double(s), 1e-12)<<"large multiplication with overflow\n";
  r = vnl_rational(p1 * p2, p3 * p4);
  s = vnl_rational(p1 * p5, p3 * p4);
  p = r / s;
  EXPECT_EQ(p.numerator() == p2 && p.denominator() == p5, true)<<"large division without overflow\n";
  r = vnl_rational(p1 * p2, p3 * p4);
  s = vnl_rational(p3 * p4, p1 * p5);
  p = r / s;
  ASSERT_NEAR( p, double(r) / double(s), 1e-12)<<"large division with overflow\n";
}

TEST(vnl_rational, long_64)
{
  long l1 = 1234321234321L, l2 = 2 * l1, l3 = 123456787654321L, l4 = l3 + 1;
  // denom = 2*num
  // relatively prime 
  vnl_rational r(-l1, -l2) , s(l3, -l4) , p;
  EXPECT_EQ(vnl_math::isfinite(r), true);
  EXPECT_EQ(vnl_math::isnan(s), false);
  EXPECT_EQ(r.numerator() == 1 && r.denominator() == 2, true)<<"simplify\n";
  EXPECT_EQ(s.numerator() == -l3 && s.denominator() == l4, true)<<"sign in numerator\n";
  // The 10-digit numbers below are prime numbers, and small enough so that the multiplications in the constructors do
  // not overflow (at least, on systems where "long" is 64 bit)
  long p1 = 1999999117L, p2 = 1999999121L, p3 = 1999999151L, p4 = 1999999171L, p5 = 1999999207L;
  r = vnl_rational(p1 * p2, p3 * p4);
  s = vnl_rational(p4 * p3, p1 * p5);
  p = r * s;
  EXPECT_EQ(p.numerator() == p2 && p.denominator() == p5, true)<<"large multiplication without overflow\n";
  r = vnl_rational(p1 * p2, p3 * p4);
  s = vnl_rational(p1 * p5, p4 * p3);
  p = r * s;
  ASSERT_NEAR(p, double(r) * double(s), 1e-7)<<"large multiplication with overflow\n";
  r = vnl_rational(p1 * p2, p3 * p4);
  s = vnl_rational(p1 * p5, p4 * p3);
  p = r / s;
  EXPECT_EQ(p.numerator() == p2 && p.denominator() == p5, true)<<"large division without overflow\n";
  r = vnl_rational(p1 * p2, p3 * p4);
  s = vnl_rational(p4 * p3, p1 * p5);
  p = r / s;
  ASSERT_NEAR(p, double(r) / double(s), 1e-7)<<"large division with overflow\n";
}


TEST(vnl_rational, approx)
{
  vnl_rational d(1.0 / 3.0); // explicit constructor from double
  EXPECT_EQ(d, vnl_rational(1, 3));
  d = vnl_rational(-5.0 / 7);
  EXPECT_EQ(d, vnl_rational(-5, 7));
  d = vnl_rational(0.42857142857);
  EXPECT_EQ(d, vnl_rational(3, 7));
  d = vnl_rational(-1.23456);
  EXPECT_EQ(d, vnl_rational(-123456, 100000));
  vnl_rational pi = vnl_rational(vnl_math::pi);
  auto pi_a = double(pi);
  EXPECT_EQ(pi_a - vnl_math::pi < 1e-18 && vnl_math::pi - pi_a < 1e-18, true)<<"pi\n";
  std::cout << "Best rational approximation of pi: " << pi << " = " << pi_a << '\n'
            << "Compare this with pi in 20 decimals:                     " << vnl_math::pi << std::endl;
}

/*
TEST(vnl_rational, determinant)
{
  vnl_matrix_fixed<vnl_rational, 3, 3> m;
  m[0][0] = vnl_rational(1, 3);
  m[0][1] = vnl_rational(2, 7);
  m[0][2] = vnl_rational(2, 5);
  m[1][0] = vnl_rational(-1, 2);
  m[1][1] = vnl_rational(1, 4);
  m[1][2] = vnl_rational(6, 7);
  m[2][0] = vnl_rational(2, 3);
  m[2][1] = vnl_rational(1, 5);
  m[2][2] = vnl_rational(5, 2);
  std::cout << "rational matrix:\n" << m << "determinant = " << vnl_det(m[0], m[1], m[2]) << std::endl;
  EXPECT_EQ(vnl_det(m[0], m[1], m[2]), vnl_rational(16609, 29400))<<"determinant\n";
}
 */

TEST(vnl_rational, sqrt)
{
  vnl_rational d(16, 9);
  EXPECT_EQ(vnl_sqrt(d), vnl_rational(4, 3));
  d = vnl_sqrt(vnl_rational(2L));
  double sqrt2 = std::sqrt(2.0), sqrt_2 = double(d);
  std::cout << "Best rational approximation of sqrt(2): " << d << " = " << sqrt_2 << '\n'
            << "Compare this with sqrt(2) in 20 decimals:                     " << sqrt2 << std::endl;
  EXPECT_EQ(sqrt2 - sqrt_2 < 1e-18 && sqrt_2 - sqrt2 < 1e-18, true);
}

/*
TEST(vnl_rational, zero_one)
{
  vnl_rational n = vnl_numeric_traits<vnl_rational>::zero;
  std::cout << "zero = " << n << '\n';
  EXPECT_EQ( n, 0L);
  vnl_rational u = vnl_numeric_traits<vnl_rational>::one;
  std::cout << "one  = " << u << '\n';
  EXPECT_EQ(u, 1L);
}
 */


