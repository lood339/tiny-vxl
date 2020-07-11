// This is core/vnl/tests/test_vector_fixed_ref.cxx
#include <algorithm>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_cross.h>

#include <gtest/gtest.h>

TEST(vnl_vector_fixed, test_int)
{
    std::cout << "***********************\n"
    << "Testing vnl_vector_fixed<int, N>\n"
    << "***********************\n";
    
    //////////////////
    // CONSTRUCTORS //
    //////////////////
    
    
    vnl_vector_fixed<int, 2> v0;
    EXPECT_EQ(v0.size(), 2);
    
    vnl_vector_fixed<int, 2> v1;
    EXPECT_EQ(v1.size(), 2);
    
    //  vnl_vector_fixed(unsigned int len, T const& v0);
    vnl_vector_fixed<int, 2> v2(2, 2);
    EXPECT_EQ((v2.get(0) == 2 && v2.get(1) == 2), true);
    
    //  vnl_vector_fixed(unsigned int len, int n, T const values[]);
    int vcvalues[] = { 1, 0 };
    vnl_vector_fixed<int, 2> vc(vcvalues);
    EXPECT_EQ((vc(0) == 1 && vc(1) == 0), true);
    
    //  vnl_vector_fixed(T const* data_block,unsigned int n);
    vnl_vector_fixed<int, 2> vb(vcvalues);
    EXPECT_EQ((vb(0) == 1 && vb(1) == 0), true);
    
    //  vnl_vector_fixed(vnl_vector_fixed<T> const&);
    vnl_vector_fixed<int, 2> v_copy(vb);
    EXPECT_EQ((v_copy(0) == 1 && v_copy(1) == 0), true);
    
    //// test constructors, accessors
    v1[0] = 2;
    v1[1] = 2;
    EXPECT_EQ((v1 == v2), true);
    EXPECT_EQ((v2.put(1, 3), v2.get(1)), 3);
    EXPECT_EQ(v2.get(1), 3);
    EXPECT_EQ((v1.fill(3), (v1.get(0) == 3 && v1.get(1) == 3)), true);
    EXPECT_EQ((v2.fill(2), (v2.get(0) == 2 && v2.get(1) == 2)), true);
    int v3values[] = { 1, 2, 3 };
    vnl_vector_fixed<int, 3> v3(v3values);
    EXPECT_EQ((v3.get(0) == 1 && v3.get(1) == 2 && v3.get(2) == 3), true);
    
    vnl_vector_fixed<int, 3> v4(v3);
    EXPECT_EQ(v3, v4);
    
#if VNL_CONFIG_CHECK_BOUNDS
    
    {
        bool exceptionThrownAndCaught = false;
        try
        {
            v0.get(25);
        } // Raise out of bounds exception.
        catch (...)
        {
            exceptionThrownAndCaught = true;
        }
        TEST("Out of bounds get()", exceptionThrownAndCaught, true);
        
        exceptionThrownAndCaught = false;
        try
        {
            v0.put(25, 0);
        } // Raise out of bounds exception.
        catch (...)
        {
            exceptionThrownAndCaught = true;
        }
        TEST("Out of bounds put()", exceptionThrownAndCaught, true);
    }
#endif
    
    //// test additions and subtractions
    EXPECT_EQ(((v0 = v2 + 3), (v0.get(0) == 5 && v0.get(1) == 5)), true);
    EXPECT_EQ(((v0 = 3 + v2), (v0.get(0) == 5 && v0.get(1) == 5)), true);
    EXPECT_EQ((v0 += (-3), (v0.get(0) == 2 && v0.get(1) == 2)), true);
    EXPECT_EQ((v0 -= (-3), (v0.get(0) == 5 && v0.get(1) == 5)), true);
    EXPECT_EQ(((v0 = v2 - 3), (v0.get(0) == -1 && v0.get(1) == -1)), true);
    EXPECT_EQ(((v0 = 3 - v2), (v0.get(0) == 1 && v0.get(1) == 1)), true);
    EXPECT_EQ((v0 = -v2, (v0.get(0) == -2 && v0.get(1) == -2)), true);
    
    
    vnl_vector_fixed<int, 2> v5;
    v0 = v2;
    EXPECT_EQ(((v5 = v0 + v2), (v5.get(0) == 4 && v5.get(1) == 4)), true);
    EXPECT_EQ(((v5 = v0 - v2), (v5.get(0) == 0 && v5.get(1) == 0)), true);
    EXPECT_EQ(((v0 += v2), (v0.get(0) == 4 && v0.get(1) == 4)), true);
    EXPECT_EQ(((v0 -= v2), (v0.get(0) == 2 && v0.get(1) == 2)), true);
    
    //// test multiplications and divisions
    EXPECT_EQ(((v4 = v3 * 5), (v4.get(0) == 5 && v4.get(1) == 10 && v4.get(2) == 15)), true);
    EXPECT_EQ(((v4 = 5 * v3), (v4.get(0) == 5 && v4.get(1) == 10 && v4.get(2) == 15)), true);
    EXPECT_EQ(((v3 *= 5), (v3 == v4)), true);
    EXPECT_EQ(((v4 = v3 / 5), (v4.get(0) == 1 && v4.get(1) == 2 && v4.get(2) == 3)), true);
    EXPECT_EQ(((v3 /= 5), (v3 == v4)), true);
    
    //// additional tests
    int vvalues[] = { 0, -2, 2, 0 };
    vnl_vector_fixed<int, 4> v(vvalues);
   
    
    EXPECT_EQ((v(0) == 0 && v(1) == -2 && v(2) == 2 && v(3) == 0), true);
    
    EXPECT_EQ(v.max_value(), 2);
    EXPECT_EQ(v.min_value(), -2);
    EXPECT_EQ(v.arg_max(), 2);
    EXPECT_EQ(v.arg_min(), 1);
    
    
    {
        vnl_vector_fixed<int, 4> v22(v);
        vnl_vector_fixed<int, 4> v1;
        vnl_vector<int> v2;
        EXPECT_EQ(
                  ((v1 = element_product(v, v)), (v1(0) == 0 && v1(1) == 4 && v1(2) == 4 && v1(3) == 0)),
                  true);
        EXPECT_EQ(
                  ((v22.fill(2)), (v1 = element_quotient(v, v22)), (v1(0) == 0 && v1(1) == -1 && v1(2) == 1 && v1(3) == 0)),
                  true);
        EXPECT_EQ(((v2 = v.extract(1, 3)), (v2.size() == 1 && v2(0) == v(3))), true);
        EXPECT_EQ(((v2.fill(4)), (v.update(v2, 3)), (v(0) == 0 && v(1) == -2 && v(2) == 2 && v(3) == 4)), true);
        
    }
    
    { // new scope to reuse variables
        int vvalues[] = { 1, 0, 0, 0 };
        vnl_vector_fixed<int, 4> v(vvalues);
        int v1values[] = { 1, 0, 0 };
        int v2values[] = { 0, 1, 0 };
        int v3values[] = { 0, 0, 1 };
        vnl_vector_fixed<int, 3> v1(v1values);
        vnl_vector_fixed<int, 3> v2(v2values);
        vnl_vector_fixed<int, 3> v3(v3values);
        EXPECT_EQ((dot_product(v1, v2) == 0 && dot_product(v1, v3) == 0 && dot_product(v2, v3) == 0), true);
        EXPECT_EQ((dot_product(v1, v1) == 1 && dot_product(v2, v2) == 1 && dot_product(v3, v3) == 1), true);
       
        EXPECT_EQ(vnl_cross_3d(v1, v2), v3);
        EXPECT_EQ(vnl_cross_3d(v2, v3), v1);
        EXPECT_EQ(vnl_cross_3d(v1, v3), -v2);
    }
    
    
    {
        int vvalues[] = { 1, 2, 3 };
        vnl_vector_fixed<int, 3> v(vvalues);
        vnl_matrix<int> m = outer_product(v, v);
        EXPECT_EQ(
                  (m(0, 0) == 1 && m(0, 1) == 2 && m(0, 2) == 3 && m(1, 0) == 2 && m(1, 1) == 4 && m(1, 2) == 6 &&
                   m(2, 0) == 3 && m(2, 1) == 6 && m(2, 2) == 9),
                  true);
    }
   
    {
        int vvalues[] = { 1, 0, 0, 0 };
        vnl_vector_fixed<int, 4> v(vvalues);
        EXPECT_EQ((v.squared_magnitude() == 1), true);
        EXPECT_EQ((v.magnitude() == 1), true);
    }
    
    //////////
    // FLIP //
    //////////
    
    {
        int vvalues[] = { 0, 1, 2, 3 };
        vnl_vector_fixed<int, 4> v(vvalues);
        EXPECT_EQ((0 == v[0] && 1 == v[1] && 2 == v[2] && 3 == v[3]), true);
        v.flip();
        EXPECT_EQ((3 == v[0] && 2 == v[1] && 1 == v[2] && 0 == v[3]), true);
    }
}

