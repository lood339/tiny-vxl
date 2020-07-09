// This is core/vnl/tests/test_matrix_fixed.cxx
#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <iostream>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_vector.h>
#include <gtest/gtest.h>
/*
#include <vnl/vnl_double_3x3.h>
#include <vnl/vnl_double_3.h>
#include <vnl/vnl_double_2x2.h>
#include <vnl/vnl_float_2x2.h>
#include <vnl/vnl_int_2x2.h>
*/


#undef printf // to work around a bug in libintl.h

bool verbose_malloc = false;
int malloc_count = 0;

// FIXME: Win32 will have different operator new in vnl dll from
// the one generated here, so this test fails - RWMC.
# define reset_count malloc_count = 0
#if defined(_WIN32)
# define check_count TEST("mallocs (no test)",true,true)
#else
# define check_count TEST("mallocs",malloc_count<=1,true)
#endif


// This function is used in testing later.
template< typename T, unsigned int n >
T sum_vector(const vnl_vector_fixed<T,n> &v) { return v.sum(); }

TEST(vnl_matrix_fixed, size)
{
    vnl_matrix_fixed<double,3,4> m;
    EXPECT_EQ(sizeof(m), sizeof(double[12]) )<< "memory footprint \n";
}


TEST(vnl_matrix_fixed, multiply)
{
    double data_m1[6] = {
    1, 2,
    3, 4,
    5, 6
    };
    double data_m2[8] = {
    2, 3, 4, 5,
    6, 7, 8, 9
    };
    double data_v1[2] = {
    7,
    8
    };

    vnl_matrix_fixed<double,3,2> m1( data_m1 );
    vnl_matrix_fixed<double,2,4> m2( data_m2 );
    vnl_vector_fixed<double,2> v1( data_v1 );
    
    vnl_matrix_fixed<double,3,4> mr = m1*m2;
    EXPECT_EQ(
       mr(0,0) == 14 && mr(0,1) == 17 && mr(0,2) == 20 && mr(0,3) == 23 &&
       mr(1,0) == 30 && mr(1,1) == 37 && mr(1,2) == 44 && mr(1,3) == 51 &&
       mr(2,0) == 46 && mr(2,1) == 57 && mr(2,2) == 68 && mr(2,3) == 79, true)<<"Matrix-matrix multiply\n";
    
    vnl_vector_fixed<double,3> vr = m1*v1;
    EXPECT_EQ(vr(0) == 23 && vr(1) == 53 && vr(2) == 83, true)<<"Matrix-vector multiply\n";
}

TEST(vnl_matrix_fixed, test_int)
{
    std::cout << "*********************************\n"
           << "Testing vnl_matrix_fixed<int,x,x>\n"
           << "*********************************" << std::endl;

    //////////////////
    // CONSTRUCTORS //
    //////////////////

    vnl_matrix_fixed<int,2,2> m0;
    EXPECT_EQ((m0.rows()==2 && m0.cols()==2), true);
    vnl_matrix_fixed<int,3,4> m1;
    EXPECT_EQ((m1.rows()==3 && m1.cols()==4), true);
    
    using vnl_int_2x2 = vnl_matrix_fixed<int, 2, 2>;
    vnl_int_2x2 m2(2);
    EXPECT_EQ((m2(0,0)==2 && m2(0,1)==2 && m2(1,0)==2 && m2(1,1)==2), true)<<"vnl_int_2x2 m2(2)\n";
       
    
    EXPECT_EQ(
       (m2 = vnl_int_2x2(2),
        (m2(0,0)==2 && m2(0,1)==2 && m2(1,0)==2 && m2(1,1)==2)), true);
    const vnl_int_2x2 ma = m2;
    EXPECT_EQ(
       (ma(0,0)==2 && ma(0,1)==2 && ma(1,0)==2 && ma(1,1)==2), true);
    vnl_int_2x2 mb = m2;
    EXPECT_EQ(
       (mb(0,0) = 0,
        mb(0,0)==0 && mb(0,1)==2 && mb(1,0)==2 && mb(1,1)==2), true);
    int mcvalues[4] = {1, 2, 3};
    vnl_int_2x2 mc(mcvalues);
    EXPECT_EQ((mc(0,0)==1 && mc(0,1)==2 && mc(1,0)==3 && mc(1,1)==0), true);
       
    EXPECT_EQ(
       (m0= vnl_int_2x2(2),
        (m0(0,0)==2 && m0(0,1)==2 && m0(1,0)==2 && m0(1,1)==2)), true);
    EXPECT_EQ((m0 == m2), true);
    EXPECT_EQ((m0 == m2), true);
    
  ///////////////
  // ACCESSORS //
  ///////////////
/*
#if VNL_CONFIG_CHECK_BOUNDS

  {
  // Get
  bool exceptionThrownAndCaught = false;
  try { m0.get(0,25); }  // Raise out of bounds exception.
  catch(...) { exceptionThrownAndCaught = true; }
  TEST("Out of bounds get(0,25)", exceptionThrownAndCaught, true);

  exceptionThrownAndCaught = false;
  try { m0.get(25,0); }  // Raise out of bounds exception.
  catch(...) { exceptionThrownAndCaught = true; }
  TEST("Out of bounds get(25,0)", exceptionThrownAndCaught, true);

  exceptionThrownAndCaught = false;
  try { m0.get(25,25); }  // Raise out of bounds exception.
  catch(...) { exceptionThrownAndCaught = true; }
  TEST("Out of bounds get(25,25)", exceptionThrownAndCaught, true);

  // Put
  exceptionThrownAndCaught = false;
  try { m0.put(0,25,0); }  // Raise out of bounds exception.
  catch(...) { exceptionThrownAndCaught = true; }
  TEST("Out of bounds put(0,25,0)", exceptionThrownAndCaught, true);

  exceptionThrownAndCaught = false;
  try { m0.put(25,0,0); }  // Raise out of bounds exception.
  catch(...) { exceptionThrownAndCaught = true; }
  TEST("Out of bounds put(25,0,0)", exceptionThrownAndCaught, true);

  exceptionThrownAndCaught = false;
  try { m0.put(25,25,0); }  // Raise out of bounds exception.
  catch(...) { exceptionThrownAndCaught = true; }
  TEST("Out of bounds put(25,25,0)", exceptionThrownAndCaught, true);

  }

#endif
 */
    m2(0, 0) = 2;
    m2(1, 1) = 3;
    EXPECT_EQ(m2(1,1), 3);
    
    int v2_data[] = {2,3};
    vnl_vector<int> temp(v2_data, 2);
    EXPECT_EQ(m2.get_diagonal(), temp)<<"m2.get_diagonal()\n";

    
    EXPECT_EQ((m0 == m2), false);
    EXPECT_EQ((m0 != m2), true);
    EXPECT_EQ(
       (m1.fill(3),
        (m1(0,0)==3 && m1(1,1)==3 && m1(2,2)==3 && m1(2,3)==3)), true)<<"m1.fill(3)\n";
    EXPECT_EQ(
       (m2.fill(2),
        (m2(0,0)==2 && m2(0,1)==2 && m2(1,0)==2 && m2(1,1)==2)), true)<<"m2.fill(2)\n";
    EXPECT_EQ((m0=m2, (m0==m2)), true)<<"m0=m2\n";

    
    // test additions and subtractions
    //auto tt = m2 + 3;
    
    EXPECT_EQ(
       ((m0=m2+3),
        (m0(0,0)==5 && m0(0,1)==5 && m0(1,0)==5 && m0(1,1)==5)), true)<<"m0=m2+3\n";
    EXPECT_EQ(
       ((m0=3+m2),
        (m0(0,0)==5 && m0(0,1)==5 && m0(1,0)==5 && m0(1,1)==5)), true)<<"m0=3+m2\n";
    EXPECT_EQ(
       (m0+=(-3),
        (m0(0,0)==2 && m0(0,1)==2 && m0(1,0)==2 && m0(1,1)==2)), true)<<"m0+=(-3) \n";
    EXPECT_EQ(
       (m0-=(-3),
        (m0(0,0)==5 && m0(0,1)==5 && m0(1,0)==5 && m0(1,1)==5)), true)<<"m0-=(-3) \n";
    EXPECT_EQ(
       ((m0=m2-3),
        (m0(0,0)==-1 && m0(0,1)==-1 && m0(1,0)==-1 && m0(1,1)==-1)), true)<<"m0=m2-3 \n";
    EXPECT_EQ(
       ((m0=3-m2),
        (m0(0,0)==1 && m0(0,1)==1 && m0(1,0)==1 && m0(1,1)==1)), true)<<"m0=3-m2 \n";
    EXPECT_EQ(
       (m0= -m2,
        (m0(0,0)==-2 && m0(0,1)==-2 && m0(1,0)==-2 && m0(1,1)==-2)), true)<<"m0= -m2\n";

    
    vnl_int_2x2 m5;
    m0 = m2;
    EXPECT_EQ(
       ((m5=m0+m2),
        (m5(0,0)==4 && m5(0,1)==4 && m5(1,0)==4 && m5(1,1)==4)), true)<<"m5=m0+m2\n";
    EXPECT_EQ(
       ((m5=m0-m2),
        (m5(0,0)==0 && m5(0,1)==0 && m5(1,0)==0 && m5(1,1)==0)), true)<<"m5=m0-m2\n";
    EXPECT_EQ(
       ((m0+=m2),
        (m0(0,0)==4 && m0(0,1)==4 && m0(1,0)==4 && m0(1,1)==4)), true)<<"m0+=m2\n";
    EXPECT_EQ(
       ((m0-=m2),
        (m0(0,0)==2 && m0(0,1)==2 && m0(1,0)==2 && m0(1,1)==2)), true)<<"m0-=m2\n";

   
    // test multiplications and divisions
    m2(0,0) = 1; m2(0,1) = 2; m2(1,0) = 3;
    EXPECT_EQ(
       ((m0=m2*5),
        (m0(0,0)==5 && m0(0,1)==10 && m0(1,0)==15)), true)<<"m0=m2*5\n";
    EXPECT_EQ(
       ((m0=5*m2),
        (m0(0,0)==5 && m0(0,1)==10 && m0(1,0)==15)), true)<<"m0=5*m2\n";
    EXPECT_EQ(((m2*=5), (m2== m0)), true)<<"m2*=5\n";
    EXPECT_EQ(
       ((m0=m2/5),
        (m0(0,0)==1 && m0(0,1)==2 && m0(1,0)==3)), true)<<"m0=m2/5\n";
    EXPECT_EQ(((m2/=5), (m2==m0)), true)<<"m2/=5\n";

    
    int m6values [] = {1,2,3,4};
    vnl_int_2x2 m6(m6values);
    EXPECT_EQ(m6(1,1), 4)<<"vnl_int_2x2 m6({1,2,3,4})\n";
    int m7values [] = {5,6,7,8};
    vnl_int_2x2 m7(m7values);
    EXPECT_EQ( m7(1,1), 8)<<"vnl_int_2x2 m7({5,6,7,8})\n";
    EXPECT_EQ(
       ((m5=m6*m7),
        (m5(0,0)==19 && m5(0,1)==22 && m5(1,0)==43 && m5(1,1)==50)), true)<<"m5=m6*m7\n";
    EXPECT_EQ(
       ((m6*=m7),
        (m6(0,0)==19 && m6(0,1)==22 && m6(1,0)==43 && m6(1,1)==50)), true)<<"m6*=m7\n";
    
    /////////////////////////////////////////////////////////////////
    // Test `flatten_row_major` and `flatten_column_major` Methods //
    /////////////////////////////////////////////////////////////////

    
    {
        int data[16] = { 0,  1,  2,  3,
                         4,  5,  6,  7,
                         8,  9, 10, 11,
                        12, 13, 14, 15};

        vnl_vector<int> flat(data, 16);

        vnl_matrix_fixed<int, 4, 4> sq(data);
        vnl_matrix_fixed<int, 2, 8> lg(data);
        vnl_matrix_fixed<int, 8, 2> wd(data);

        EXPECT_EQ(flat.is_equal(sq.flatten_row_major().as_vector(), 10e-6), true);
        EXPECT_EQ(flat.is_equal(lg.flatten_row_major().as_vector(), 10e-6), true);
        EXPECT_EQ(flat.is_equal(wd.flatten_row_major().as_vector(), 10e-6), true);

        EXPECT_EQ(flat.is_equal(sq.transpose().flatten_column_major().as_vector(), 10e-6), true);
        EXPECT_EQ(flat.is_equal(lg.transpose().flatten_column_major().as_vector(), 10e-6), true);
        EXPECT_EQ(flat.is_equal(wd.transpose().flatten_column_major().as_vector(), 10e-6), true);
    }
 

   
  // additional tests
  int mvalues [] = {0,-2,2,0};
  vnl_int_2x2 m(mvalues); m0 = m;
  vnl_matrix<int> m3;
  EXPECT_EQ(
       (m(0,0)==0 && m(0,1)==-2 && m(1,0)==2 && m(1,1)==0), true);
  EXPECT_EQ(m.max_value(),  2);
  EXPECT_EQ(m.min_value(), -2);
  EXPECT_EQ(m.arg_max(),   2);
  EXPECT_EQ(m.arg_min(),   1);
    
  EXPECT_EQ(
       ((m0 = m.transpose()),
        (m0(0,0)==0 && m0(0,1)==2 && m0(1,0)==-2 && m0(1,1)==0)), true);
  EXPECT_EQ(
       ((m0 = element_product(m,m)),
        (m0(0,0)==0 && m0(0,1)==4 && m0(1,0)==4 && m0(1,1)==0)), true);
  EXPECT_EQ(
       ((m2 = 2),
        (m0 = element_quotient(m,m2)),
        (m0(0,0)==0 && m0(0,1)==-1 && m0(1,0)==1 && m0(1,1)==0)), true);
  EXPECT_EQ(
       ((m3 = m.extract(1,1,1,1)),
        (m3.rows()==1 && m3.columns()==1 && m3(0,0)==m(1,1))), true);
 
    EXPECT_EQ(
       ((m3=4),
        (m.update(m3,1,1)),
        (m(0,0)==0 && m(0,1)==-2 && m(1,0)==2 && m(1,1)==4)), true)<<"m.update([4],1,1)\n";
}

TEST(vnl_matrix_fixed, test_float)
{
    std::cout << "***********************************\n"
    << "Testing vnl_matrix_fixed<float,x,x>\n"
    << "***********************************" << std::endl;
    vnl_matrix_fixed<float,2,2> d0;
    EXPECT_EQ((d0.rows()==2 && d0.columns()==2), true);
    vnl_matrix_fixed<float,3,4> d1;
    EXPECT_EQ((d1.rows()==3 && d1.columns()==4), true);
    
    
    vnl_matrix_fixed<float, 2, 2> d2(2.0f);
    //vnl_float_2x2 d2(2.0);
    EXPECT_EQ(
         (d2.get(0,0)==2.0 && d2.get(0,1)==2.0 && d2.get(1,0)==2.0 && d2.get(1,1)==2.0), true);
    EXPECT_EQ((d0=2.0,
                    (d0.get(0,0)==2.0 && d0.get(0,1)==2.0 && d0.get(1,0)==2.0 && d0.get(1,1)==2.0)), true);
    EXPECT_EQ((d0 == d2), true);
    EXPECT_EQ((d0==d2), true);
    EXPECT_EQ((d2.put(1,1,(float)3.0),d2.get(1,1)), (float)3.0);
    EXPECT_EQ(d2.get(1,1), (float)3.0);
    float v2_data[] = {2.f,3.f};
    EXPECT_EQ(d2.get_diagonal(), vnl_vector<float>(2,2,v2_data));
    EXPECT_EQ((d0 == d2), false);
    EXPECT_EQ((d0 != d2), true);
    EXPECT_EQ(
         (d1.fill(3.0),
          (d1.get(0,0)==3.0 && d1.get(1,1)==3.0 && d1.get(2,2)==3.0 && d1.get(2,3)==3.0)), true);
    EXPECT_EQ(
         (d2.fill(2.0),
          (d2.get(0,0)==2.0 && d2.get(0,1)==2.0 && d2.get(1,0)==2.0 && d2.get(1,1)==2.0)), true);
    EXPECT_EQ((d0=d2,  (d0==d2)), true);
    
    
    // test additions and subtractions
    EXPECT_EQ(
         ((d0=d2+(float)3.0),
          (d0.get(0,0)==5.0 && d0.get(0,1)==5.0 && d0.get(1,0)==5.0 && d0.get(1,1)==5.0)), true);
    EXPECT_EQ(
         (d0+=(-3.0),
          (d0.get(0,0)==2.0 && d0.get(0,1)==2.0 && d0.get(1,0)==2.0 && d0.get(1,1)==2.0)), true);
    //vnl_float_2x2 d5;
    vnl_matrix<float> d5(2, 2);
    EXPECT_EQ(
         ((d5=d0+d2),
          (d5.get(0,0)==4.0 && d5.get(0,1)==4.0 && d5.get(1,0)==4.0 && d5.get(1,1)==4.0)), true);
    EXPECT_EQ(
         ((d0+=d2),
          (d0.get(0,0)==4.0 && d0.get(0,1)==4.0 && d0.get(1,0)==4.0 && d0.get(1,1)==4.0)), true);
    
    // test multiplications and divisions
    d2(0,0) = 1; d2(0,1) = 2; d2(1,0) = 3;
    EXPECT_EQ(
         ((d0=d2*5.0f),
          (d0.get(0,0)==5 && d0.get(0,1)==10 && d0.get(1,0)==15)), true);
    EXPECT_EQ(
         ((d0=5.0f*d2),
          (d0.get(0,0)==5 && d0.get(0,1)==10 && d0.get(1,0)==15)), true);
    EXPECT_EQ(((d2*=5.0f), (d2== d0)), true);
    EXPECT_EQ(
         ((d0=d2/5.0f),
          (d0.get(0,0)==1 && d0.get(0,1)==2 && d0.get(1,0)==3)), true);
    EXPECT_EQ(((d2/=5.0f), (d2==d0)), true);
    float d6values [] = {1.0f,2.0f,
        3.0f,4.0f};
    //vnl_float_2x2 d6(d6values);
    vnl_matrix_fixed<float, 2, 2> d6(d6values);
    EXPECT_EQ(d6.get(1,1), 4.0);
    float d7values [] = {5.0,6.0,
        7.0,8.0};
    //vnl_float_2x2 d7(d7values);
    vnl_matrix<float> d7(d7values, 2, 2);
    EXPECT_EQ(d7.get(1,1), 8.0);
    EXPECT_EQ(((d5=d6*d7),
                      (d5.get(0,0)==19.0 && d5.get(0,1)==22.0 && d5.get(1,0)==43.0 && d5.get(1,1)==50.0)), true);
    EXPECT_EQ(((d6*=d7),
                    (d6.get(0,0)==19.0 && d6.get(0,1)==22.0 && d6.get(1,0)==43.0 && d6.get(1,1)==50.0)), true);
    
    // additional tests
    //vnl_float_2x2 m1, m2;
    vnl_matrix_fixed<float, 2, 2> m1, m2;
    float mvalues [] = {0,-2,2,0};
    vnl_matrix_fixed<float, 2, 2> m(mvalues);
    m1 = m; m2 = m;
    vnl_matrix<float> m3;
    EXPECT_EQ(
         (m(0,0)==0 && m(0,1)==-2 && m(1,0)==2 && m(1,1)==0), true);
    EXPECT_EQ( m.max_value(),  2);
    EXPECT_EQ( m.min_value(), -2);
    EXPECT_EQ(m.arg_max(),   2);
    EXPECT_EQ(m.arg_min(),   1);
    EXPECT_EQ(
         ((m1 = m.transpose()),
          (m1(0,0)==0 && m1(0,1)==2 && m1(1,0)==-2 && m1(1,1)==0)), true);
    EXPECT_EQ(
         ((m1 = element_product(m,m)),
          (m1(0,0)==0 && m1(0,1)==4 && m1(1,0)==4 && m1(1,1)==0)), true);
    EXPECT_EQ(
         ((m2 = 2),
          (m1 = element_quotient(m,m2)),
          (m1(0,0)==0 && m1(0,1)==-1 && m1(1,0)==1 && m1(1,1)==0)), true);
    EXPECT_EQ(
         ((m3 = m.extract(1,1,1,1)),
          (m3.rows()==1 && m3.columns()==1 && m3(0,0)==m(1,1))), true);
    EXPECT_EQ(
         ((m3=4),
          (m.update(m3,1,1)),
          (m(0,0)==0 && m(0,1)==-2 && m(1,0)==2 && m(1,1)==4)), true);
}

TEST(vnl_matrix_fixed, test_double)
{
    std::cout << "************************************\n"
    << "Testing vnl_matrix_fixed<double,x,x>\n"
    << "************************************" << std::endl;
    vnl_matrix_fixed<double,2,2> d0;
    EXPECT_EQ((d0.rows()==2 && d0.columns()==2), true);
    vnl_matrix_fixed<double,3,4> d1;
    EXPECT_EQ((d1.rows()==3 && d1.columns()==4), true);
    //vnl_double_2x2 d2(2.0);
    vnl_matrix_fixed<double, 2, 2> d2(2.0);
    EXPECT_EQ(
         (d2.get(0,0)==2.0 && d2.get(0,1)==2.0 && d2.get(1,0)==2.0 && d2.get(1,1)==2.0), true);
    EXPECT_EQ((d0=2.0,
                    (d0.get(0,0)==2.0 && d0.get(0,1)==2.0 && d0.get(1,0)==2.0 && d0.get(1,1)==2.0)), true);
    EXPECT_EQ((d0 == d2), true);
    EXPECT_EQ((d0==d2), true);
    EXPECT_EQ((d2.put(1,1,3.0),d2.get(1,1)), 3.0);
    EXPECT_EQ(d2.get(1,1), 3.0);
    double v2_data[] = {2.0,3.0};
    EXPECT_EQ(d2.get_diagonal(), vnl_vector<double>(2,2,v2_data));
    EXPECT_EQ((d0 == d2), false);
    EXPECT_EQ((d0 != d2), true);
    EXPECT_EQ(
         (d1.fill(3.0),
          (d1.get(0,0)==3.0 && d1.get(1,1)==3.0 && d1.get(2,2)==3.0 && d1.get(2,3)==3.0)), true);
    EXPECT_EQ(
         (d2.fill(2.0),
          (d2.get(0,0)==2.0 && d2.get(0,1)==2.0 && d2.get(1,0)==2.0 && d2.get(1,1)==2.0)), true);
    EXPECT_EQ((d0=d2,  (d0==d2)), true);
    
    // test additions and subtractions
    EXPECT_EQ(
         ((d0=d2+3.0),
          (d0.get(0,0)==5.0 && d0.get(0,1)==5.0 && d0.get(1,0)==5.0 && d0.get(1,1)==5.0)), true);
    EXPECT_EQ(
         (d0+=(-3.0),
          (d0.get(0,0)==2.0 && d0.get(0,1)==2.0 && d0.get(1,0)==2.0 && d0.get(1,1)==2.0)), true);
    //vnl_double_2x2 d5;
    vnl_matrix_fixed<double, 2, 2> d5;
    EXPECT_EQ(
         ((d5=d0+d2),
          (d5.get(0,0)==4.0 && d5.get(0,1)==4.0 && d5.get(1,0)==4.0 && d5.get(1,1)==4.0)), true);
    EXPECT_EQ(
         ((d0+=d2),
          (d0.get(0,0)==4.0 && d0.get(0,1)==4.0 && d0.get(1,0)==4.0 && d0.get(1,1)==4.0)), true);
    
    // test multiplications and divisions
    d2(0,0) = 1; d2(0,1) = 2; d2(1,0) = 3;
    EXPECT_EQ(
         ((d0=d2*5.0),
          (d0.get(0,0)==5 && d0.get(0,1)==10 && d0.get(1,0)==15)), true);
    EXPECT_EQ(
         ((d0=5.0*d2),
          (d0.get(0,0)==5 && d0.get(0,1)==10 && d0.get(1,0)==15)), true);
    EXPECT_EQ(((d2*=5.0), (d2== d0)), true);
    EXPECT_EQ(
         ((d0=d2/5.0),
          (d0.get(0,0)==1 && d0.get(0,1)==2 && d0.get(1,0)==3)), true);
    EXPECT_EQ(((d2/=5.0), (d2==d0)), true);
    double d6values [] = {1.0,2.0,
        3.0,4.0};
    //vnl_double_2x2 d6(d6values);
    vnl_matrix_fixed<double, 2, 2> d6(d6values);
    EXPECT_EQ(d6.get(1,1), 4.0);
    double d7values [] = {5.0,6.0,
        7.0,8.0};
    //vnl_double_2x2 d7(d7values);
    vnl_matrix_fixed<double, 2, 2> d7(d7values);
    EXPECT_EQ( d7.get(1,1), 8.0);
    EXPECT_EQ(((d5=d6*d7),
                      (d5.get(0,0)==19.0 && d5.get(0,1)==22.0 && d5.get(1,0)==43.0 && d5.get(1,1)==50.0)), true);
    EXPECT_EQ(((d6*=d7),
                    (d6.get(0,0)==19.0 && d6.get(0,1)==22.0 && d6.get(1,0)==43.0 && d6.get(1,1)==50.0)), true);
    
    // apply sqrt to every element
    double d8values [] = {0.0, 1.0, 9.0, 16.0};
    //vnl_double_2x2 d8(d8values);
    vnl_matrix_fixed<double, 2, 2> d8(d8values);
    d8 = d8.apply(std::sqrt);
    EXPECT_EQ(d8[0][0]==0 && d8[0][1]==1 && d8[1][0]==3 && d8[1][1]==4, true);
    
    {
        vnl_matrix_fixed<double,4,20> m(1.);
        vnl_vector_fixed<double,4> vr = m.apply_rowwise(sum_vector);
        for (unsigned int i = 0; i < vr.size(); ++i)
            EXPECT_EQ(vr.get(i), 20.);
        vnl_vector_fixed<double,20> vc = m.apply_columnwise(sum_vector);
        for (unsigned int i = 0; i < vc.size(); ++i)
            EXPECT_EQ(vc.get(i), 4.);
    }
    
    // normalizations
    d8.normalize_rows();
    EXPECT_EQ(d8[0][0]==0 && d8[0][1]==1, true);
    ASSERT_NEAR(d8[1][0], 0.6, 1e-12);
    ASSERT_NEAR(d8[1][1], 0.8, 1e-12);
    d8.normalize_columns();
    EXPECT_EQ(d8[0][0]==0 && d8[1][0]==1, true);
}

namespace {
    
    template<class T>
    void
    test_extract( T* )
    {
        
    }
    
} // end anonymous namespace

TEST(vnl_matrix_fixed, extract)
{
    vnl_matrix_fixed<double,2,6> m;
    m(0,0)=1; m(0,1)=2; m(0,2)=3; m(0,3)=4; m(0,4)=5; m(0,5) = 11;
    m(1,0)=6; m(1,1)=7; m(1,2)=8; m(1,3)=9; m(1,4)=0; m(1,5) = 12;
    
    vnl_matrix_fixed<double,1,3> r;
    vnl_matrix<double> rm = r.as_matrix();
    m.extract(rm, 1, 2 );
    EXPECT_EQ(rm(0,0)==8 && rm(0,1)==9 && rm(0,2)==0, true );
}


