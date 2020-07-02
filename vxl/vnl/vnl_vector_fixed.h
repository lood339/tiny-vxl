//
//  vnl_vector_fixed.h
//  vxl
//
//  Created by Jimmy Chen on 2020-06-14.
//  Copyright © 2020 Nowhere Planet. All rights reserved.
//

#ifndef vnl_vector_fixed_h
#define vnl_vector_fixed_h

#include <Eigen/Dense>
#include <vnl/vnl_error.h>

// fixed length column vector

template <typename T, unsigned int n> class vnl_vector_fixed;
template <typename T> class vnl_vector;
template <typename T> class vnl_matrix;

template <typename T, unsigned int n>
class vnl_vector_fixed : public Eigen::Matrix<T, n, 1, Eigen::ColMajor>
{
    static constexpr size_t num_bytes = n*sizeof(T);
    using base_class = Eigen::Matrix<T, n, 1, Eigen::ColMajor>;
    
public:
    //: Construct an uninitialized n-vector
    vnl_vector_fixed() = default;
    
    //: Construct a fixed-n-vector initialized from \a datablck
    //  The data \e must have enough data. No checks performed.
    explicit vnl_vector_fixed( const T* datablck )
    {
        std::memcpy( this->data(), datablck, num_bytes);
    }
    
    //: Convenience constructor for 2-D vectors
    // While this constructor is sometimes useful, consider using
    // vnl_double_2 or vnl_float_2 instead.
    vnl_vector_fixed( const T& x0, const T& x1 )
    {
        if ( n != 2 )
        {
#ifndef NDEBUG
            vnl_error_vector_dimension("vnl_vector_fixed()", 2, n);
#endif
            return;
        }
        T* data = this->data();
        data[0] = x0; data[1] = x1;
    }
    
    //: Convenience constructor for 3-D vectors
    // While this constructor is sometimes useful, consider using
    // vnl_double_3 or vnl_float_3 instead.
    vnl_vector_fixed( const T& x0, const T& x1, const T& x2 )
    {
        if ( n != 3 )
        {
#ifndef NDEBUG
            vnl_error_vector_dimension("vnl_vector_fixed()", 3, n);
#endif
            return;
        }
        T* data = this->data();
        data[0] = x0; data[1] = x1; data[2] = x2;
    }
    
    //: Convenience constructor for 4-D vectors
    vnl_vector_fixed( const T& x0, const T& x1, const T& x2, const T& x3 )
    {
        if ( n != 4 )
        {
#ifndef NDEBUG
            vnl_error_vector_dimension("vnl_vector_fixed()", 4, n);
#endif
            return;
        }
        T* data = this->data();
        data[0] = x0; data[1] = x1; data[2] = x2; data[3] = x3;
    }
    
    // This constructor allows us to construct vnl_vector_fixed from Eigen expressions
    template<typename OtherDerived>
    vnl_vector_fixed(const Eigen::MatrixBase<OtherDerived>& other):
        base_class(other)
    {}
    
    // This method allows us to assign Eigen expressions to vnl_vector_fixed
    template<typename OtherDerived>
    vnl_vector_fixed& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->base_class::operator=(other);
        return *this;
    }
    
    //: Access the contiguous block storing the elements in the vector.
    //  O(1).
    //  data_block()[0] is the first element of the vector
    T const* data_block() const
    {
        return this->data();
    }
    
    //: Access the contiguous block storing the elements in the vector.
    //  O(1).
    //  data_block()[0] is the first element of the vector
    T      * data_block()
    {
        return this->data();
    }
    
    vnl_vector<T> as_vector() const { return vnl_vector<T>(data_block(), n); }
};

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned n>
inline T dot_product( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    return dot_product( a.as_vector(), b.as_vector() );
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned n>
inline T dot_product( const vnl_vector_fixed<T,n>& a, const vnl_vector<T>& b )
{
    return dot_product( a.as_vector(), b );
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned n>
inline T dot_product( const vnl_vector<T>& a, const vnl_vector_fixed<T,n>& b )
{
    return dot_product( a, b.as_vector() );
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_matrix<T> outer_product( const vnl_vector<T>& a, const vnl_vector_fixed<T,n>& b )
{
    return outer_product( a, b.as_vector());
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_matrix<T> outer_product( const vnl_vector_fixed<T,n>& a, const vnl_vector<T>& b )
{
    return outer_product( a.as_vector(), b);
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned n>
inline T angle( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    return angle( a.as_vector(), b.as_vector() );
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned n>
inline T angle( const vnl_vector_fixed<T,n>& a, const vnl_vector<T>& b )
{
    return angle( a.as_vector(), b );
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned n>
inline T angle( const vnl_vector<T>& a, const vnl_vector_fixed<T,n>& b )
{
    return angle( a, b.as_vector() );
}


template<class T, unsigned n>
inline T vnl_vector_ssd( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    return vnl_vector_ssd( a.as_vector(), b.as_vector() );
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned n>
inline T vnl_vector_ssd( const vnl_vector_fixed<T,n>& a, const vnl_vector<T>& b )
{
    return vnl_vector_ssd( a.as_vector(), b );
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned n>
inline T vnl_vector_ssd( const vnl_vector<T>& a, const vnl_vector_fixed<T,n>& b )
{
    return vnl_vector_ssd( a, b.as_vector() );
}

#endif /* vnl_vector_fixed_h */
