// This is core/vnl/tests/test_vector_fixed_ref.cxx
#include <algorithm>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_matrix_fixed.h>

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

static bool
float_equal(const float & f1, const float & f2)
{
    return std::fabs(f1 - f2) < 1.0e-6;
}

TEST(vnl_vector_fixed, test_float)
{
    std::cout << "*************************\n"
    << "Testing vnl_vector<float>\n"
    << "*************************\n";
    
    using vnl_float_3 = vnl_vector_fixed<float, 3>;
    using vnl_float_4 = vnl_vector_fixed<float, 4>;
    
    //// test constructors, accessors
    vnl_vector_fixed<float, 0> v0;
    EXPECT_EQ(v0.size(), 0);
    vnl_vector_fixed<float, 2> v1;
    EXPECT_EQ(v1.size(), 2);
    vnl_vector_fixed<float, 2> v2;
    v2.fill(2);
    EXPECT_EQ((v2.get(0) == 2 && v2.get(1) == 2), true);
    
    float vcvalues[2] = { 1 };
    vnl_vector_fixed<float, 2> vc(vcvalues);
    EXPECT_EQ((vc(0) == 1 && vc(1) == 0), true);
    EXPECT_EQ((v2.put(1, 3), v2.get(1)), 3);
    EXPECT_EQ(v2.get(1), 3);
    
    EXPECT_EQ((v1.fill(3), (v1.get(0) == 3 && v1.get(1) == 3)), true);
    EXPECT_EQ((v2.fill(2), (v2.get(0) == 2 && v2.get(1) == 2)), true);
    vnl_vector_fixed<float, 3> v3 = vnl_float_3(1.f, 2.f, 3.f);
    EXPECT_EQ((v3.get(0) == 1 && v3.get(1) == 2 && v3.get(2) == 3), true);
    vnl_vector_fixed<float, 3> v4(v3);
    EXPECT_EQ(v3, v4);
    std::cout << &v0 << " == " << v0 << std::endl;
   
    
    
    
    {
        vnl_vector_fixed<float, 2> v5;
        vnl_vector_fixed<float, 2> v6 = v2;
        EXPECT_EQ(((v5 = v6 + v2), (v5.get(0) == 4 && v5.get(1) == 4)), true);
        EXPECT_EQ(((v5 = v6 - v2), (v5.get(0) == 0 && v5.get(1) == 0)), true);
        EXPECT_EQ(((v6 += v2), (v6.get(0) == 4 && v6.get(1) == 4)), true);
        EXPECT_EQ(((v6 -= v2), (v6.get(0) == 2 && v6.get(1) == 2)), true);
    }
    
    v4 = v3;
    
    //// test multiplications and divisions
    
    EXPECT_EQ(((v4 = v3 * 5.0f), (v4.get(0) == 5 && v4.get(1) == 10 && v4.get(2) == 15)), true);
    
    EXPECT_EQ(((v4 = 5.0f * v3), (v4.get(0) == 5 && v4.get(1) == 10 && v4.get(2) == 15)), true);
    EXPECT_EQ(((v3 *= 5.0f), (v3 == v4)), true);
    EXPECT_EQ(((v4 = v3 / 5.0f), (v4.get(0) == 1 && v4.get(1) == 2 && v4.get(2) == 3)), true);
    EXPECT_EQ(((v3 /= 5.0f), (v3 == v4)), true);
    
    //// additional tests
    //  vnl_vector<float> v(4,4,0,-2,2,0); no var args with floats
    float vvalues[] = { 0, -2, 2, 0 };
    vnl_vector_fixed<float, 4> v(vvalues);
    
    v[0] = 0;
    v[1] = -2;
    v[2] = 2;
    v[3] = 0;
 
    EXPECT_EQ(v.max_value(), 2);
    EXPECT_EQ(v.min_value(), -2);
    EXPECT_EQ(v.arg_max(), 2);
    EXPECT_EQ(v.arg_min(), 1);
    
    {
        vnl_vector_fixed<float, 4> v5;
        vnl_vector<float> v1(4);
        
        EXPECT_EQ(((v1 = element_product(v, v)), (v1(0) == 0 && v1(1) == 4 && v1(2) == 4 && v1(3) == 0)),
                  true);
        EXPECT_EQ(((v5.fill(2)), (v1 = element_quotient(v, v5)), (v1(0) == 0 && v1(1) == -1 && v1(2) == 1 &&        v1(3) == 0)),
                  true);
        EXPECT_EQ(((v1 = v.extract(1, 3)), (v1.size() == 1 && v1(0) == v(3))), true);
        EXPECT_EQ(((v1.fill(4)), (v.update(v1, 3)), (v(0) == 0 && v(1) == -2 && v(2) == 2 && v(3) == 4)), true);
    }
    
    { // new scope to reuse variables
        float vvalues[] = { 1, 0, 0, 0 };
        vnl_vector_fixed<float, 4> v(vvalues);
        v[0] = 1;
        v[1] = 0;
        v[2] = 0;
        v[3] = 0;
        EXPECT_EQ(
                  (v(0) == v[0] && v[0] == 1 && v(1) == v[1] && v[1] == 0 && v(2) == v[2] && v[2] == 0 && v(3) == v[3] &&
                   v[3] == 0),
                  true);
        vnl_vector_fixed<float, 3> v1;
        v1.fill(0);
        v1[0] = 1.f;
        vnl_vector_fixed<float, 3> v2;
        v2.fill(0);
        v2[1] = 1.f;
        vnl_vector_fixed<float, 3> v3;
        v3.fill(0);
        v3[2] = 1.f;
        EXPECT_EQ((dot_product(v1, v2) == 0 && dot_product(v1, v3) == 0 && dot_product(v2, v3) == 0), true);
        EXPECT_EQ(vnl_cross_3d(v1, v2), v3);
        EXPECT_EQ(vnl_cross_3d(v2, v3), v1);
        EXPECT_EQ(vnl_cross_3d(v1, v3), -v2);
    }
    
    
    {
        vnl_float_3 v(1.f, 2.f, 3.f);
        float v2_data[4] = {1.f, 2.f, 3.f, 4.f};
        vnl_vector_fixed<float, 4> v2(v2_data);
        
        vnl_matrix_fixed<float, 3, 4> m = outer_product(v, v2.as_vector());
        EXPECT_EQ((m(0, 0) == 1 && m(0, 1) == 2 && m(0, 2) == 3 && m(0, 3) == 4 && m(1, 0) == 2 && m(1, 1) == 4 &&
                   m(1, 2) == 6 && m(1, 3) == 8 && m(2, 0) == 3 && m(2, 1) == 6 && m(2, 2) == 9 && m(2, 3) == 12),
                  true);
    }
    
    {
        vnl_float_3 v(1.f, 2.f, 3.f);
        EXPECT_EQ(v.size(), 3);
        v[0] = 1.f;
        v[1] = 2.f;
        v[2] = 3.f;
        EXPECT_EQ(v[0], 1);
        EXPECT_EQ(v[1], 2);
        EXPECT_EQ(v[2], 3);
        vnl_float_3 v1;
        v1.fill(0);
        v1[0] = 1.f;
        std::cout << "v1 = " << v1 << std::endl;
        vnl_float_3 v2;
        v2.fill(0);
        v2[1] = 1.f;
        std::cout << "v2 = " << v2 << std::endl;
        vnl_float_3 v3;
        v3.fill(0);
        v3[0] = -0.5f;
        v3[2] = 0.5f;
        std::cout << "v3 = " << v3 << std::endl << "v1 - v2 = " << v1 - v2 << std::endl;
        double ang = angle(v1, v2);
        std::cout << "angle(v1,v2) = " << ang << std::endl;
        ang *= vnl_math::deg_per_rad; // == 180/pi
        std::cout << "angle(v1,v2) in degrees = " << ang << std::endl
        << "v1.size()=" << v1.size() << std::endl
        << "v2.size()=" << v2.size() << std::endl
        << "vnl_cross_3d(v1,v2) = " << vnl_cross_3d(v1, v2) << std::endl;
        ASSERT_NEAR(ang, 90.0, 1e-5)<<"angle(v1,v2)\n";
        double ang2 = angle(v1, v3);
        std::cout << "angle(v1,v3) = " << ang << std::endl;
        ang2 *= vnl_math::deg_per_rad; // == 180/pi
        std::cout << "angle(v1,v3) in degrees = " << ang2 << std::endl;
        ASSERT_NEAR(ang2, 135.0, 1e-6)<<"angle(v1,v3)\n";
        ASSERT_NEAR(v.mean(), 2.0, 1e-6)<<"mean\n";
    }
    
    {
        double vvalues[] = { 1, 0, 0, 0 };
        vnl_vector_fixed<double, 4> vd;
        vnl_vector_fixed<double, 4> v(vvalues);
        v[0] = 1;
        v[1] = 0;
        v[2] = 0;
        v[3] = 0;
        EXPECT_EQ(v.squared_magnitude(), 1);
        EXPECT_EQ(v.magnitude(), 1);
        // Trying to track down test failure in Intel 10.0 compiler
        vd = 4.0 * v;
        std::cout << "vd.normalize() is " << vd.normalize() << " and v is " << v << "\n" << std::flush;
        std::cout << "vd.normalize() - v is " << vd.normalize() - v << "\n" << std::flush;
        EXPECT_EQ((vd = 4.0 * v, vd.normalize(), vd), v);
    }
    
    {
        float vvalues[] = { -7, -2, -3, -4 };
        vnl_vector_fixed<float, 4> v(vvalues);
        v[0] = -7;
        v[1] = -2;
        v[2] = -3;
        v[3] = -4;
        EXPECT_EQ((v(0) == -7 && v(1) == -2 && v(2) == -3 && v(3) == -4), true);
        
        EXPECT_EQ(v.max_value(), -2);
        EXPECT_EQ(v.min_value(), -7);
        EXPECT_EQ(v.arg_max(), 1);
        EXPECT_EQ(v.arg_min(), 0);
    }
    
    {
        double vvalues[] = { -7, -2, -3, -4 };
        vnl_vector_fixed<double, 4> v(vvalues);
        v[0] = -7;
        v[1] = -2;
        v[2] = -3;
        v[3] = -4;
        EXPECT_EQ((v(0) == -7 && v(1) == -2 && v(2) == -3 && v(3) == -4), true);
        
        EXPECT_EQ(v.max_value(), -2);
        EXPECT_EQ(v.min_value(), -7);
        EXPECT_EQ(v.arg_max(), 1);
        EXPECT_EQ(v.arg_min(), 0);
    }
    
    
    // We should test odd sized vector's for the special SSE2 handling of
    // different sizes.
    {
        float vvalues[] = { -7, -2, -3, -4, 5 };
        vnl_vector_fixed<float, 5> v(vvalues);
        v[0] = -7;
        v[1] = -2;
        v[2] = -3;
        v[3] = -4;
        v[4] = 5;
        EXPECT_EQ((v(0) == -7 && v(1) == -2 && v(2) == -3 && v(3) == -4 && v(4) == 5), true);
        
        EXPECT_EQ(v.max_value(), 5);
        EXPECT_EQ(v.min_value(), -7);
        EXPECT_EQ(v.arg_max(), 4);
        EXPECT_EQ(v.arg_min(), 0);
    }
    
    
    {
        float vvalues[] = { -7, -2, -3, -4, 5, 2 };
        vnl_vector_fixed<float, 6> v(vvalues);
        v[0] = -7;
        v[1] = -2;
        v[2] = -3;
        v[3] = -4;
        v[4] = 5;
        v[5] = 2;
        EXPECT_EQ((v(0) == -7 && v(1) == -2 && v(2) == -3 && v(3) == -4 && v(4) == 5 && v(5) == 2), true);
        
        EXPECT_EQ(v.max_value(), 5);
        EXPECT_EQ(v.min_value(), -7);
        EXPECT_EQ(v.arg_max(), 4);
        EXPECT_EQ(v.arg_min(), 0);
    }
    {
        double vvalues[] = { -7, -2, -3, -4, 5 };
        vnl_vector_fixed<double, 5> v(vvalues);
        v[0] = -7;
        v[1] = -2;
        v[2] = -3;
        v[3] = -4;
        v[4] = 5;
        EXPECT_EQ((v(0) == -7 && v(1) == -2 && v(2) == -3 && v(3) == -4 && v(4) == 5), true);
        
        EXPECT_EQ(v.max_value(), 5);
        EXPECT_EQ(v.min_value(), -7);
        EXPECT_EQ(v.arg_max(), 4);
        EXPECT_EQ(v.arg_min(), 0);
    }
}

TEST(vnl_vector_fixed, conversion)
{
    bool check;
    {
        // convert from a vnl_vector to a block array:
        int v1values[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
        vnl_vector_fixed<int, 12> v1(v1values);
        const int * data = v1.data_block();
        {
            check = true;
            for (int d = 0; d < 12; d++)
                if (data[d] != d + 1)
                    check = false;
        }
        EXPECT_EQ(check, true)<<"(const int*) m.data_block\n";
        
        typedef int block[12];
        const block & v2 = *((const block *)data);
        {
            check = true;
            for (int i = 0; i < 12; i++)
                if (v1(i) != v2[i])
                    check = false;
        }
        EXPECT_EQ(check, true)<<"matrix(i)==block[i]\n";
        
        // convert from a block array to a vnl_vector:
        block b1;
        for (int i = 0; i < 12; i++)
            b1[i] = i;
        data = ((const int *)b1);
        {
            check = true;
            for (int d = 0; d < 12; d++)
                if (data[d] != d)
                    check = false;
        }
        EXPECT_EQ(check, true)<<"(const int*) block\n";
        vnl_vector_fixed<int, 12> b2(data);
        {
            check = true;
            for (int i = 0; i < 12; i++)
                if (b1[i] != b2(i))
                    check = false;
        }
        EXPECT_EQ(check, true)<<"block[i]==matrix(i)\n";
    }
}
