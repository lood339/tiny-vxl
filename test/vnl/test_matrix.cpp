// This is core/vnl/tests/test_matrix.cxx
#include <iostream>
#include <cmath>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h> // necessary for tests of methods set_diagonal() and get_diagonal()
//#include <vnl/vnl_copy.h>

#include <gtest/gtest.h>

// This function is used in testing later.
template< typename T >
T sum_vector(const vnl_vector<T> &v) { return v.sum(); }


TEST(vnl_matrix, test_int)
{
    std::cout << "***********************\n"
           << "Testing vnl_matrix<int>\n"
           << "***********************" << std::endl;

    //////////////////
    // CONSTRUCTORS //
    //////////////////
    
    vnl_matrix<int> m0(2,2);
    EXPECT_EQ((m0.rows()==2 && m0.cols()==2), true)<<"vnl_matrix<int> m0(2,2)\n";
    
    
    vnl_matrix<int> m1(3,4);
    EXPECT_EQ((m1.rows()==3 && m1.cols()==4), true)<<"vnl_matrix<int> m1(3,4)\n";
    vnl_matrix<int> m2(2,2,2);
    EXPECT_EQ(
       (m2(0,0)==2 && m2(0,1)==2 && m2(1,0)==2 && m2(1,1)==2), true)<<"vnl_matrix<int> m2(2,2,2)\n";
   
    EXPECT_EQ(
       (m2 = vnl_matrix<int>(2,2, 2),
        (m2(0,0)==2 && m2(0,1)==2 && m2(1,0)==2 && m2(1,1)==2)), true)<<"m2 = vnl_matrix<int>(2,2, 2)\n";
    
    const vnl_matrix<int> ma = m2;
    EXPECT_EQ(
       (ma(0,0)==2 && ma(0,1)==2 && ma(1,0)==2 && ma(1,1)==2), true)<<"(const vnl_matrix)(i,j)\n";
    vnl_matrix<int> mb = m2;
    EXPECT_EQ(
       (mb(0,0) = 0,
        mb(0,0)==0 && mb(0,1)==2 && mb(1,0)==2 && mb(1,1)==2), true)<<"(vnl_matrix)(i,j)\n";
    
    int mcvalues[4] = {1, 2, 3};
    vnl_matrix<int> mc(2,2, 4, mcvalues);
    
    EXPECT_EQ(
       (mc(0,0)==1 && mc(0,1)==2 && mc(1,0)==3 && mc(1,1)==0), true)<<"vnl_matrix<int> mc(2,2, 4,int[])\n";
    EXPECT_EQ(
       (m0=2,
        (m0(0,0)==2 && m0(0,1)==2 && m0(1,0)==2 && m0(1,1)==2)), true)<<"m0=2\n";
    EXPECT_EQ((m0 == m2), true)<<"m0 == m2\n";
    EXPECT_EQ((m0 == m2), true)<<"(m0 == m2)\n";
    
    

  ///////////////
  // ACCESSORS //
  ///////////////

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

   
  ///////////////////////////////
  // Row and Column Operations //
  ///////////////////////////////

  {
    int some_data[2*3] = { 1, 2, 3,
                           4, 5, 6 };
    vnl_matrix<int> mrc(some_data, 2, 3);
    vnl_vector<int> v;
      
      
      // test get row and column
    EXPECT_EQ((v = mrc.row(0), (v(0)==1 && v(1)==2 && v(2)==3)), true);
    EXPECT_EQ((v = mrc.row(1),(v(0)==4 && v(1)==5 && v(2)==6)), true);
    EXPECT_EQ((v = mrc.col(0), (v(0)==1 && v(1)==4)), true);
    EXPECT_EQ((v = mrc.col(1), (v(0)==2 && v(1)==5)), true);
    EXPECT_EQ((v = mrc.col(2), (v(0)==3 && v(1)==6)), true);
      
      
     
    // Indices: {1, 0, 1}
    vnl_vector<unsigned int> indices(3, 1);
    indices.put(1,0);
    vnl_matrix<int> mi;

    EXPECT_EQ(
         (mi = mrc.get_rows(indices),
         (mi.get_row(0)==mrc.get_row(1)
       && mi.get_row(1)==mrc.get_row(0)
       && mi.get_row(2)==mrc.get_row(1))), true)<<"mi = mrc.get_rows(indices)\n";
    EXPECT_EQ(
         (mi = mrc.get_columns(indices),
         (mi.get_column(0)==mrc.get_column(1)
       && mi.get_column(1)==mrc.get_column(0)
       && mi.get_column(2)==mrc.get_column(1))), true)<<"mi = mrc.get_columns(indices)\n";
  }
 
    int v2_data[] = {2,3};
    m2(0,0) = v2_data[0];
    m2(1,1) = v2_data[1];
    EXPECT_EQ(m2.diagonal(), vnl_vector<int>(2,2,v2_data))<<"m2.diagonal()\n";
   
    EXPECT_EQ((m0 == m2), false);
    EXPECT_EQ((m0 != m2), true);
    EXPECT_EQ((m0 == m2), false);
    
    EXPECT_EQ((m1.fill(3), (m1(0,0)==3 && m1(1,1)==3 && m1(2,2)==3 && m1(2,3)==3)), true)<<"m1.fill(3)\n";
    EXPECT_EQ((m2.fill(2), (m2(0,0)==2 && m2(0,1)==2 && m2(1,0)==2 && m2(1,1)==2)), true);
    EXPECT_EQ(vnl_matrix<int>(2,2).fill(2), m2)<<"vnl_matrix<int>(2,2).fill(2)\n";
   
    EXPECT_EQ((m0.fill_diagonal(3),
               (m0(0,0)==3 && m0(1,1)==3)), true)<<"m0.fill_diagonal(3)\n";
    std::cout<<m0<<std::endl;
    int m0values [] = {7,9};
    EXPECT_EQ((m0.diagonal() = (vnl_vector<int>(2,2,m0values)),
                (m0(0,0)==7 && m0(1,1)==9 && m0(0,1)==2 && m0(1,0)==2)), true)<<"m0.diagonal(vnl_vector<int>))\n";
       
    

    
    /////////////////////////////////////////////////////////////////
    // Test `flatten_row_major` and `flatten_column_major` Methods //
    /////////////////////////////////////////////////////////////////

    {
        int data[16] = { 0,  1,  2,  3,
                         4,  5,  6,  7,
                         8,  9, 10, 11,
                        12, 13, 14, 15};

        vnl_vector<int> flat(data, 16);

        vnl_matrix<int> sq(data, 4, 4);
        vnl_matrix<int> lg(data, 2, 8);
        vnl_matrix<int> wd(data, 8, 2);

        EXPECT_EQ(flat.is_equal(sq.flatten_row_major(), 10e-6), true)<<"sq.flatten_row_major\n";
        EXPECT_EQ(flat.is_equal(lg.flatten_row_major(), 10e-6), true)<<"lg.flatten_row_major\n";
        EXPECT_EQ(flat.is_equal(wd.flatten_row_major(), 10e-6), true)<<"wd.flatten_row_major\n";

        EXPECT_EQ(flat.is_equal(sq.transpose().flatten_column_major(), 10e-6), true)<<"sq.flatten_column_major\n";
        EXPECT_EQ(flat.is_equal(lg.transpose().flatten_column_major(), 10e-6), true)<<"lg.flatten_column_major\n";
        EXPECT_EQ(flat.is_equal(wd.transpose().flatten_column_major(), 10e-6), true)<<"wd.flatten_column_major\n";
    }

    int m3values [] = {1,2,3};
    vnl_matrix<int> m3(1,3,3, m3values);
    EXPECT_EQ(
       (m3(0,0)==1 && m3(0,1)==2 && m3(0,2)==3), true)<<"m3(1,3,3,{1,2,3})\n";
    vnl_matrix<int> m4(m3);
    EXPECT_EQ((m3==m4), true)<<"vnl_matrix<int> m4(m3)\n";
    EXPECT_EQ((m0=m2, (m0==m2)), true)<<"m0=m2\n";

  // test additions and subtractions
  EXPECT_EQ(
       ((m0=m2+3),
        (m0(0,0)==5 && m0(0,1)==5 && m0(1,0)==5 && m0(1,1)==5)), true);
  EXPECT_EQ(
       ((m0=3+m2),
        (m0(0,0)==5 && m0(0,1)==5 && m0(1,0)==5 && m0(1,1)==5)), true);
  EXPECT_EQ(
       (m0+=(-3),
        (m0(0,0)==2 && m0(0,1)==2 && m0(1,0)==2 && m0(1,1)==2)), true);
  EXPECT_EQ(
       (m0-=(-3),
        (m0(0,0)==5 && m0(0,1)==5 && m0(1,0)==5 && m0(1,1)==5)), true);
  EXPECT_EQ(
       ((m0=m2-3),
        (m0(0,0)==-1 && m0(0,1)==-1 && m0(1,0)==-1 && m0(1,1)==-1)), true);
  EXPECT_EQ(
       ((m0=3-m2),
        (m0(0,0)==1 && m0(0,1)==1 && m0(1,0)==1 && m0(1,1)==1)), true);
  EXPECT_EQ(
       (m0= -m2,
        (m0(0,0)==-2 && m0(0,1)==-2 && m0(1,0)==-2 && m0(1,1)==-2)), true);

    
    vnl_matrix<int> m5(2,2);
    m0 = m2;
    EXPECT_EQ(
       ((m5=m0+m2),
        (m5(0,0)==4 && m5(0,1)==4 && m5(1,0)==4 && m5(1,1)==4)), true);
    EXPECT_EQ(
       ((m5=m0-m2),
        (m5(0,0)==0 && m5(0,1)==0 && m5(1,0)==0 && m5(1,1)==0)), true);
    EXPECT_EQ(
       ((m0+=m2),
        (m0(0,0)==4 && m0(0,1)==4 && m0(1,0)==4 && m0(1,1)==4)), true);
    EXPECT_EQ(
       ((m0-=m2),
        (m0(0,0)==2 && m0(0,1)==2 && m0(1,0)==2 && m0(1,1)==2)), true);
    

    /// test multiplications and divisions
    EXPECT_EQ(
       ((m4=m3*5),
        (m4(0,0)==5 && m4(0,1)==10 && m4(0,2)==15)), true);
    EXPECT_EQ(
       ((m4=5*m3),
        (m4(0,0)==5 && m4(0,1)==10 && m4(0,2)==15)), true);
    EXPECT_EQ(((m3*=5), (m3== m4)), true);
    EXPECT_EQ(
       ((m4=m3/5),
        (m4(0,0)==1 && m4(0,1)==2 && m4(0,2)==3)), true);
    EXPECT_EQ(((m3/=5), (m3==m4)), true);

    int m6values [] = {1,2,3,4};
    vnl_matrix<int> m6(2,2,4,m6values);
    EXPECT_EQ(m6(1,1), 4)<<"vnl_matrix<int> m6(2,2,4,{1,2,3,4})\n";
    int m7values [] = {5,6,7,8};
    vnl_matrix<int> m7(2,2,4,m7values);
    EXPECT_EQ(m7(1,1), 8)<<"vnl_matrix<int> m7(2,2,4,{5,6,7,8})\n";
    EXPECT_EQ(
       ((m5=m6*m7),
        (m5(0,0)==19 && m5(0,1)==22 && m5(1,0)==43 && m5(1,1)==50)), true)<<"m5=m6*m7\n";
    EXPECT_EQ(
       ((m6*=m7),
        (m6(0,0)==19 && m6(0,1)==22 && m6(1,0)==43 && m6(1,1)==50)), true)<<"m6*=m7\n";
    int c0values [] = {1,0};
    vnl_matrix<int> c0(2,1,2,c0values);
    vnl_matrix<int> c1;
    EXPECT_EQ(
       ((c1=m6*c0),
        c1.rows()==c0.rows() && c1.cols()==c0.cols() &&
        c1(0,0)==19 && c1(1,0)==43), true)<<"c1=m6*c0\n";
    int r0values [] = {1,0};
    vnl_matrix<int> r0(1,2,2,r0values);
    vnl_matrix<int> r1;
    EXPECT_EQ(
       ((r1=r0*m6),
        r1.rows()==r0.rows() && r1.cols()==r0.cols() &&
        r1(0,0)==19 && r1(0,1)==22), true)<<"r1=r0*m6\n";
    EXPECT_EQ(
       ((r0*=m6), r0==r1), true)<<"r0*=m6\n";
    EXPECT_EQ(
       ((m6*=c0), c1==m6), true)<<"m6*=c0\n";

   
    // additional tests
    int mvalues [] = {0,-2,2,0};
    vnl_matrix<int> m(2,2,4,mvalues);
    m0 = m; m1 = m;
    EXPECT_EQ(
       (m(0,0)==0 && m(0,1)==-2 && m(1,0)==2 && m(1,1)==0), true);
    EXPECT_EQ(
       ((m1 = m.transpose()),
        (m1(0,0)==0 && m1(0,1)==2 && m1(1,0)==-2 && m1(1,1)==0)), true);
    

    EXPECT_EQ(m.max_value(),  2);
    EXPECT_EQ(m.min_value(), -2);
    EXPECT_EQ(m.arg_max(),   2);
    EXPECT_EQ(m.arg_min(),   1);
    EXPECT_EQ(
       ((m1 = element_product(m,m)),
        (m1(0,0)==0 && m1(0,1)==4 && m1(1,0)==4 && m1(1,1)==0)), true)<<"element_product(m,m)";
    EXPECT_EQ(
       ((m2 = 2),
        (m1 = element_quotient(m,m2)),
        (m1(0,0)==0 && m1(0,1)==-1 && m1(1,0)==1 && m1(1,1)==0)), true)<<"element_quotient(m,[2])\n";
    EXPECT_EQ(
       ((m1 = m.extract(1,1,1,1)),
        (m1.rows()==1 && m1.cols()==1 && m1(0,0)==m(1,1))), true)<<"m.extract(1,1,1,1)\n";
    EXPECT_EQ(
       ((m1=4),
        (m.update(m1,1,1)),
        (m(0,0)==0 && m(0,1)==-2 && m(1,0)==2 && m(1,1)==4)), true)<<"m.update([4],1,1)\n";
    

    int vvalues[] = {1,0,0,0};
    vnl_matrix<int> v (4,1,4,vvalues);
    int v1values [] = {1,0,0};
    int v2values [] = {0,1,0};
    int v3values [] = {0,0,1};
    vnl_matrix<int> v1(3,1,3,v1values);
    vnl_matrix<int> v2(3,1,3,v2values);
    vnl_matrix<int> v3(3,1,3,v3values);
    EXPECT_EQ(
       (dot_product(v1,v2)==0 && dot_product(v1,v3)==0 && dot_product(v2,v3)==0), true)<<"dot_product(v1,v2)\n";
    v = v3;
    EXPECT_EQ((v.rows()==3 && v.cols()==1 && v==v3), true)<<"4d-v=3d-v\n";

    // Zero-size
    {
        vnl_matrix<int> m1(0,3);
        vnl_matrix<int> m2(3,4);
        vnl_matrix<int> m3(4,0);
        vnl_matrix<int> m = m1 * (m2 * m3);
        EXPECT_EQ(m.rows(), 0)<<"zero-size mult rows\n";
        EXPECT_EQ(m.cols(), 0)<<"zero-size mult cols\n";

        m = (m1 * m2) * m3;
        EXPECT_EQ(m.rows(), 0)<<"zero-size mult rows\n";
        EXPECT_EQ(m.cols(), 0)<<"zero-size mult cols\n";

        m2.clear();
        EXPECT_EQ(m2.rows(), 0)<<"zero-size after clear()\n";
        EXPECT_EQ(m2.cols(), 0)<<"zero-size after clear()";
    }
    
    {
        vnl_matrix<int> m(5,10,1);
        vnl_vector<int> vr = m.apply_rowwise(sum_vector);
        for (unsigned int i = 0; i < vr.size(); ++i)
            EXPECT_EQ(vr.get(i), 10)<<"vr.apply_rowwise(sum_vector)\n";
        vnl_vector<int> vc = m.apply_columnwise(sum_vector);
        for (unsigned int i = 0; i < vc.size(); ++i)
            EXPECT_EQ(vc.get(i), 5)<<"vc.apply_columnwise(sum_vector)\n";
    }
}

TEST(vnl_matrix, test_float)
{
    std::cout << "*************************\n"
    << "Testing vnl_matrix<float>\n"
    << "*************************" << std::endl;
    vnl_matrix<float> d0(2,2);
    EXPECT_EQ((d0.rows()==2 && d0.columns()==2), true);
    vnl_matrix<float> d1(3,4);
    EXPECT_EQ((d1.rows()==3 && d1.columns()==4), true);
    vnl_matrix<float> d2(2,2,2.0);
    EXPECT_EQ(
         (d2.get(0,0)==2.f && d2.get(0,1)==2.f && d2.get(1,0)==2.f && d2.get(1,1)==2.f), true);
    EXPECT_EQ((d0=2.f,
                    (d0.get(0,0)==2.f && d0.get(0,1)==2.f && d0.get(1,0)==2.f && d0.get(1,1)==2.f)), true);
    EXPECT_EQ((d0 == d2), true);
    EXPECT_EQ((d0==d2), true);
    EXPECT_EQ((d2.put(1,1,(float)3.0),d2.get(1,1)), (float)3.0);
    EXPECT_EQ(d2.get(1,1), (float)3.0);
    float v2_data[] = {2.f,3.f};
    EXPECT_EQ(d2.get_diagonal(), vnl_vector<float>(2,2,v2_data));
    EXPECT_EQ((d0 == d2), false);
    EXPECT_EQ((d0 != d2), true);
    EXPECT_EQ((d0==d2), false);
    EXPECT_EQ(
         (d1.fill(3.f),
          (d1.get(0,0)==3.f && d1.get(1,1)==3.f && d1.get(2,2)==3.f && d1.get(2,3)==3.f)), true);
    EXPECT_EQ(
         (d2.fill(2.f),
          (d2.get(0,0)==2.f && d2.get(0,1)==2.f && d2.get(1,0)==2.f && d2.get(1,1)==2.f)), true);
    EXPECT_EQ(vnl_matrix<float>(2,2).fill(2.f), d2);
    EXPECT_EQ(
         (d0.fill_diagonal(3.f),
          (d0.get(0,0)==3.f && d0.get(1,1)==3.f && d0.get(0,1)==2.f && d0.get(1,0)==2.f)), true);
    float d0values [] = {7.f,9.f};
    EXPECT_EQ(
         (d0.set_diagonal(vnl_vector<float>(2,2,d0values)),
          (d0.get(0,0)==7.f && d0.get(1,1)==9.f && d0.get(0,1)==2.f && d0.get(1,0)==2.f)), true);
    float d3values [] = {1.0,2.0,3.0};
    vnl_matrix<float> d3(1,3,3,d3values);
    EXPECT_EQ(
         (d3.get(0,0)==1.0 && d3.get(0,1)==2.0 && d3.get(0,2)==3.0), true);
    vnl_matrix<float> d4(d3);
    EXPECT_EQ(d3, d4);
    EXPECT_EQ((d0=d2,  (d0==d2)), true);
    EXPECT_EQ(
         ((d0=d2+(float)3.0),
          (d0.get(0,0)==5.0 && d0.get(0,1)==5.0 && d0.get(1,0)==5.0 && d0.get(1,1)==5.0)), true);
    EXPECT_EQ(
         (d0+=(-3.0),
          (d0.get(0,0)==2.0 && d0.get(0,1)==2.0 && d0.get(1,0)==2.0 && d0.get(1,1)==2.0)), true);
    vnl_matrix<float> d5(2,2);
    EXPECT_EQ(
         ((d5=d0+d2),
          (d5.get(0,0)==4.0 && d5.get(0,1)==4.0 && d5.get(1,0)==4.0 && d5.get(1,1)==4.0)), true);
    EXPECT_EQ(
         ((d0+=d2),
          (d0.get(0,0)==4.0 && d0.get(0,1)==4.0 && d0.get(1,0)==4.0 && d0.get(1,1)==4.0)), true);
    EXPECT_EQ(((d4=d3*5.0),(d4.get(0,0)==5.0 && d4.get(0,1)==10.0 && d4.get(0,2)==15.0)), true);
    EXPECT_EQ(((d3*=5.0),  (d3== d4)), true);
    float d6values [] = {1.0,2.0,
        3.0,4.0};
    vnl_matrix<float> d6(2,2,4,d6values);
    EXPECT_EQ(d6.get(1,1), 4.0);
    float d7values [] = {5.0,6.0,
        7.0,8.0};
    vnl_matrix<float> d7(2,2,4,d7values);
    EXPECT_EQ(d7.get(1,1), 8.0);
    EXPECT_EQ(((d5=d6*d7),
                      (d5.get(0,0)==19.0 && d5.get(0,1)==22.0 && d5.get(1,0)==43.0 && d5.get(1,1)==50.0)), true);
    EXPECT_EQ(((d6*=d7),
                    (d6.get(0,0)==19.0 && d6.get(0,1)==22.0 && d6.get(1,0)==43.0 && d6.get(1,1)==50.0)), true);
    
    // additional tests
    vnl_matrix<float> m0, m1, m2;
    float mvalues [] = {0,-2,2,0};
    vnl_matrix<float> m(2,2,4,mvalues);
    m0 = m; m1 = m; m2 = m;
    EXPECT_EQ(
         (m(0,0)==0 && m(0,1)==-2 && m(1,0)==2 && m(1,1)==0), true);
    EXPECT_EQ(
         ((m1 = m.transpose()),
          (m1(0,0)==0 && m1(0,1)==2 && m1(1,0)==-2 && m1(1,1)==0)), true);
    
    EXPECT_EQ(m.max_value(),  2);
    EXPECT_EQ(m.min_value(), -2);
    EXPECT_EQ(m.arg_max(),   2);
    EXPECT_EQ(m.arg_min(),   1);
    EXPECT_EQ(
         ((m1 = element_product(m,m)),
          (m1(0,0)==0 && m1(0,1)==4 && m1(1,0)==4 && m1(1,1)==0)), true);
    EXPECT_EQ(
         ((m2 = 2),
          (m1 = element_quotient(m,m2)),
          (m1(0,0)==0 && m1(0,1)==-1 && m1(1,0)==1 && m1(1,1)==0)), true);
    EXPECT_EQ(
         ((m1 = m.extract(1,1,1,1)),
          (m1.rows()==1 && m1.columns()==1 && m1(0,0)==m(1,1))), true);
    EXPECT_EQ(
         ((m1=4),
          (m.update(m1,1,1)),
          (m(0,0)==0 && m(0,1)==-2 && m(1,0)==2 && m(1,1)==4)), true);
    
    float vvalues[] = {1,0,0,0};
    vnl_matrix<float> v (4,1,4,vvalues);
    float v1values [] = {1,0,0};
    float v2values [] = {0,1,0};
    float v3values [] = {0,0,1};
    vnl_matrix<float> v1(3,1,3,v1values);
    vnl_matrix<float> v2(3,1,3,v2values);
    vnl_matrix<float> v3(3,1,3,v3values);
    EXPECT_EQ(
         (dot_product(v1,v2)==0 && dot_product(v1,v3)==0 && dot_product(v2,v3)==0), true);
    v = v3;
    EXPECT_EQ( (v.rows()==3 && v.columns()==1 && v==v3), true);
    
    v.clear();
    EXPECT_EQ(v.rows(), 0)<<"zero-size after clear()\n";
    EXPECT_EQ(v.columns(), 0)<<"zero-size after clear()\n";
    
    {
        vnl_matrix<float> m(5,10,1);
        vnl_vector<float> vr = m.apply_rowwise(sum_vector);
        for (unsigned int i = 0; i < vr.size(); ++i)
            EXPECT_EQ(vr.get(i), 10);
        vnl_vector<float> vc = m.apply_columnwise(sum_vector);
        for (unsigned int i = 0; i < vc.size(); ++i)
            EXPECT_EQ(vc.get(i), 5);
    }

    
}

TEST(vnl_matrix, test_double)
{
    std::cout << "**************************\n"
    << "Testing vnl_matrix<double>\n"
    << "**************************" << std::endl;
    vnl_matrix<double> d0(2,2);
    EXPECT_EQ((d0.rows()==2 && d0.columns()==2), true);
    vnl_matrix<double> d1(3,4);
    EXPECT_EQ((d1.rows()==3 && d1.columns()==4), true);
    vnl_matrix<double> d2(2,2,2.0);
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
    EXPECT_EQ((d0==d2), false);
    EXPECT_EQ(
         (d1.fill(3.0),
          (d1.get(0,0)==3.0 && d1.get(1,1)==3.0 && d1.get(2,2)==3.0 && d1.get(2,3)==3.0)), true);
    EXPECT_EQ(
         (d2.fill(2.0),
          (d2.get(0,0)==2.0 && d2.get(0,1)==2.0 && d2.get(1,0)==2.0 && d2.get(1,1)==2.0)), true);
    EXPECT_EQ(vnl_matrix<double>(2,2).fill(2.0), d2);
    EXPECT_EQ(
         (d0.fill_diagonal(3.0),
          (d0.get(0,0)==3.0 && d0.get(1,1)==3.0 && d0.get(0,1)==2.0 && d0.get(1,0)==2.0)), true);
    double d0values [] = {7.0,9.0};
    EXPECT_EQ(
         (d0.set_diagonal(vnl_vector<double>(2,2,d0values)),
          (d0.get(0,0)==7.0 && d0.get(1,1)==9.0 && d0.get(0,1)==2.0 && d0.get(1,0)==2.0)), true);
    double d3values [] = {1.0,2.0,3.0};
    vnl_matrix<double> d3(1,3,3,d3values);
    EXPECT_EQ(
         (d3.get(0,0)==1.0 && d3.get(0,1)==2.0 && d3.get(0,2)==3.0), true);
    vnl_matrix<double> d4(d3);
    EXPECT_EQ((d3 == d4), true);
    EXPECT_EQ((d0=d2,  (d0==d2)), true);
    EXPECT_EQ(
         ((d0=d2+3.0),
          (d0.get(0,0)==5.0 && d0.get(0,1)==5.0 && d0.get(1,0)==5.0 && d0.get(1,1)==5.0)), true);
    EXPECT_EQ(
         (d0+=(-3.0),
          (d0.get(0,0)==2.0 && d0.get(0,1)==2.0 && d0.get(1,0)==2.0 && d0.get(1,1)==2.0)), true);
    vnl_matrix<double> d5(2,2);
    EXPECT_EQ(
         ((d5=d0+d2),
          (d5.get(0,0)==4.0 && d5.get(0,1)==4.0 && d5.get(1,0)==4.0 && d5.get(1,1)==4.0)), true);
    EXPECT_EQ(
         ((d0+=d2),
          (d0.get(0,0)==4.0 && d0.get(0,1)==4.0 && d0.get(1,0)==4.0 && d0.get(1,1)==4.0)), true);
    EXPECT_EQ(((d4=d3*5.0),(d4.get(0,0)==5.0 && d4.get(0,1)==10.0 && d4.get(0,2)==15.0)), true);
    EXPECT_EQ(((d3*=5.0),  (d3== d4)), true);
    double d6values [] = {1.0,2.0,
        3.0,4.0};
    vnl_matrix<double> d6(2,2,4,d6values);
    EXPECT_EQ(d6.get(1,1), 4.0);
    double d7values [] = {5.0,6.0,
        7.0,8.0};
    vnl_matrix<double> d7(2,2,4,d7values);
    EXPECT_EQ(d7.get(1,1), 8.0);
    EXPECT_EQ(((d5=d6*d7),
                      (d5.get(0,0)==19.0 && d5.get(0,1)==22.0 && d5.get(1,0)==43.0 && d5.get(1,1)==50.0)), true);
    EXPECT_EQ(((d6*=d7),
                    (d6.get(0,0)==19.0 && d6.get(0,1)==22.0 && d6.get(1,0)==43.0 && d6.get(1,1)==50.0)), true);
    
    d0.clear();
    EXPECT_EQ(d0.rows(), 0);
    EXPECT_EQ(d0.columns(), 0);
    
    // apply sqrt to every element
    double d8values [] = {0.0, 1.0, 9.0, 16.0};
    vnl_matrix<double> d8(2,2,4,d8values);
    d8 = d8.apply(std::sqrt);
    EXPECT_EQ(d8[0][0]==0 && d8[0][1]==1 && d8[1][0]==3 && d8[1][1]==4, true);
    
    // normalizations
    d8.normalize_rows();
    EXPECT_EQ(d8[0][0]==0 && d8[0][1]==1, true);
    ASSERT_NEAR(d8[1][0], 0.6, 1e-12);
    ASSERT_NEAR(d8[1][1], 0.8, 1e-12);
    d8.normalize_columns();
    EXPECT_EQ(d8[0][0]==0 && d8[1][0]==1, true);
    
    //vnl_matrix<double> d9(2,2);
    //vnl_copy(d2, d9);
    //EXPECT_EQ("vnl_copy(S, S)", d9==d2, true);
    //vnl_matrix<float> d10(2,2);
    //vnl_copy(d2, d10);
    //d9.fill(-15.0);
    //vnl_copy(d10, d9);
    //EXPECT_EQ("vnl_copy(T, S)", d9==d2, true);
    
    {
        vnl_matrix<double> m(5,10,1);
        vnl_vector<double> vr = m.apply_rowwise(sum_vector);
        for (unsigned int i = 0; i < vr.size(); ++i)
            EXPECT_EQ(vr.get(i), 10);
        vnl_vector<double> vc = m.apply_columnwise(sum_vector);
        for (unsigned int i = 0; i < vc.size(); ++i)
            EXPECT_EQ(vc.get(i), 5);
    }
}
    


/*
#ifdef LEAK
static
void test_leak()   // use top4.1 to watch memory usage.
{
  for (;;) {       // remember to kill process.
    test_int();
    test_float();
    test_double();
  }
}
*/


namespace {
  template<class T>
  void
  test_extract( T* )
  {
    vnl_matrix<T> m( 2, 5 );
    m(0,0)=1; m(0,1)=2; m(0,2)=3; m(0,3)=4; m(0,4)=5;
    m(1,0)=6; m(1,1)=7; m(1,2)=8; m(1,3)=9; m(1,4)=0;
    std::cout << "m=\n" << m << '\n';

    vnl_matrix<T> r( 1, 3 );
    m.extract( r, 1, 2 );
    std::cout << "r=\n" << r << '\n';
    EXPECT_EQ(r(0,0)==8 && r(0,1)==9 && r(0,2)==0, true );
  }
} // end anonymous namespace

TEST(vnl_matrix, extract)
{
    test_extract( (double*)nullptr );
}


