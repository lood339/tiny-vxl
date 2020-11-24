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
#include <vnl/vnl_numeric_traits.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector_ref.h>

// fixed length column vector

template <typename T, unsigned int n> class vnl_vector_fixed;
template <typename T> class vnl_vector;
template <typename T> class vnl_matrix;

template <typename T, unsigned int n>
class vnl_vector_fixed : public Eigen::Matrix<T, n, 1, Eigen::ColMajor>
{
    static constexpr size_t num_bytes = n*sizeof(T);
    using base_class = Eigen::Matrix<T, n, 1, Eigen::ColMajor>;
    using this_class = vnl_vector_fixed<T, n>;
    using abs_t = typename vnl_numeric_traits<T>::abs_t;
    
public:
    //: Construct an uninitialized n-vector
    vnl_vector_fixed() = default;
    
    //: Constructs n-vector with all elements initialised to \a v
    explicit vnl_vector_fixed( const T& v ) { fill( v ); }
    
    //: Copy constructor
    //  The dimensions must match.
    vnl_vector_fixed( const vnl_vector_fixed<T,n>& rhs ) = default;
    vnl_vector_fixed( vnl_vector_fixed<T,n>&& rhs ) = default;
    //: Copy operator
    vnl_vector_fixed<T,n>& operator=( const vnl_vector_fixed<T,n>& rhs ) = default;
    vnl_vector_fixed<T,n>& operator=( vnl_vector_fixed<T,n>&& rhs ) = default;
    
    
    //: Construct a fixed-n-vector initialized from \a datablck
    //  The data \e must have enough data. No checks performed.
    explicit vnl_vector_fixed( const T* datablck ) { std::memcpy( this->data(), datablck, num_bytes);}
    
    //: Convenience constructor for 2-D vectors
    vnl_vector_fixed( const T& x0, const T& x1 );
    
    //: Convenience constructor for 3-D vectors
    vnl_vector_fixed( const T& x0, const T& x1, const T& x2 );
    
    //: Convenience constructor for 4-D vectors
    vnl_vector_fixed( const T& x0, const T& x1, const T& x2, const T& x3 );
    
   
    // This constructor allows us to construct vnl_vector_fixed from Eigen expressions
    template<typename OtherDerived>
    vnl_vector_fixed(const Eigen::MatrixBase<OtherDerived>& other):
        base_class(other)
    {assert(other.size() == n);}
     
    // This method allows us to assign Eigen expressions to vnl_vector_fixed
    template<typename OtherDerived>
    vnl_vector_fixed& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        assert(other.size() == n);
        this->base_class::operator=(other);
        return *this;
    }
    
    //: Copy data from a dynamic vector
    // The dimensions must match.
    vnl_vector_fixed<T,n>& operator=( const vnl_vector<T>& rhs) {
        assert( n == rhs.size() );
        std::memcpy( this->data(), rhs.data_block(), sizeof this->data() );
        return *this;
    }
    
    //: Length of the vector.
    // This is always \a n.
    unsigned int size() const { return n; }
    
    //: Put value at given position in vector.
    inline void put (unsigned int i, T const& v)
    {
#if VNL_CONFIG_CHECK_BOUNDS
        if (i >= this->size())           // If invalid index specified
            vnl_error_vector_index("put", i); // Raise exception
#endif
        (*this)(i) = v;
    }
    
    //: Get value at element i
    T get(unsigned int i) const {return (*this)(i);}
    
    //: Set all values to v
    vnl_vector_fixed& fill( T const& v )
    {
        this->setConstant(v);
        return *this;
    }
    
    //: Sets elements to ptr[i]
    //  Note: ptr[i] must be valid for i=0..size()-1
    vnl_vector_fixed& copy_in( T const * ptr )
    {
        for ( unsigned int i = 0; i < n; ++i )
            this->data()[i] = ptr[i];
        return *this;
    }
    
    //: Copy elements to ptr[i]
    //  Note: ptr[i] must be valid for i=0..size()-1
    void copy_out( T* ptr ) const
    {
        for ( unsigned int i = 0; i < n; ++i )
            ptr[i] = this->data()[i];
    }
    
    //: Sets elements to ptr[i]
    //  Note: ptr[i] must be valid for i=0..size()-1
    vnl_vector_fixed& set( T const *ptr ) { return copy_in(ptr); }
    
    //: Return reference to the element at specified index.
    // There are assert style boundary checks - #define NDEBUG to turn them off.
    //T       & operator() (unsigned int i);
    
    //: Return reference to the element at specified index.
    // There are assert style boundary checks - #define NDEBUG to turn them off.
    //T const & operator() (unsigned int i) const;
    
    
    //: Return the i-th element
    //T& operator[] (const size_t i);
    
    //: Return the i-th element
    //const T& operator[] (const size_t i) const;
    
    
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
    
    //----------------------------------------------------------------------
    // Conversion to vnl_vector_ref.
    
    // The const version of as_ref should return a const vnl_vector_ref
    // so that the vnl_vector_ref::non_const() cannot be used on
    // it. This prevents a const vnl_vector_fixed from being cast into a
    // non-const vnl_vector reference, giving a slight increase in type safety.
    
    //: Explicit conversion to a vnl_vector_ref or vnl_vector.
    // This is a cheap conversion for those functions that have an interface
    // for vnl_vector but not for vnl_vector_fixed. There is also a
    // conversion operator that should work most of the time.
    // \sa vnl_vector_ref::non_const
    vnl_vector_ref<T> as_ref() { return vnl_vector_ref<T>( n, data_block() ); }
    const vnl_vector_ref<T> as_ref() const { return vnl_vector_ref<T>( n, const_cast<T*>(data_block()));}
    
    vnl_vector<T> as_vector() const { return vnl_vector<T>(data_block(), n); }
    
    //: Type defs for iterators
    typedef T element_type;
    //: Type defs for iterators
    typedef T       *iterator;
    //: Iterator pointing to start of data
    iterator begin() { return this->data(); }
    
    //: Iterator pointing to element beyond end of data
    iterator end() { return this->data()+n; }
    
    //: Const iterator type
    typedef T const *const_iterator;
    //: Iterator pointing to start of data
    const_iterator begin() const { return this->data(); }
    //: Iterator pointing to element beyond end of data
    const_iterator end() const { return this->data()+n; }
    
    //: Apply f to each element.
    // Returns a new vector with the result.
    vnl_vector_fixed<T,n> apply(T (*f)(T));
    
    //: Apply f to each element.
    // Returns a new vector with the result.
    vnl_vector_fixed<T,n> apply(T (*f)(const T&));
    
    //:
    vnl_vector_fixed<T,n>& operator+=( T s ) { this_class::add( this->data(), s, this->data() ); return *this; }
    
    //:
    vnl_vector_fixed<T,n>& operator-=( T s ) { this_class::sub( this->data(), s, this->data() ); return *this; }
    
    //:
    vnl_vector_fixed<T,n>& operator*=( T s ) { this_class::mul( this->data(), s, this->data() ); return *this; }
    
    //:
    vnl_vector_fixed<T,n>& operator/=( T s ) { this_class::div( this->data(), s, this->data() ); return *this; }
    
    //:
    vnl_vector_fixed<T,n>& operator+=( const vnl_vector_fixed<T,n>& v ) {
        this_class::add( this->data(), v.data_block(), this->data() );
        return *this;
    }
    
    //:
    vnl_vector_fixed<T,n>& operator-=( const vnl_vector_fixed<T,n>& v ) { this_class::sub( this->data(), v.data_block(), this->data() );
        return *this;
    }
    
    //:
    vnl_vector_fixed<T,n>& operator+=( const vnl_vector<T>& v )
    {
        assert( v.size() == n );
        this_class::add( this->data(), v.data_block(), this->data() );
        return *this;
    }
    
    //:
    vnl_vector_fixed<T,n>& operator-=( const vnl_vector<T>& v )
    {
        assert( v.size() == n );
        this_class::sub( this->data(), v.data_block(), this->data() );
        return *this;
    }
    
    //:
    vnl_vector_fixed<T,n> operator-() const
    {
        vnl_vector_fixed<T,n> result;
        this_class::sub( (T)0, this->data(), result.data() );
        return result;
    }
    
    //: Returns a subvector specified by the start index and length. O(n).
    vnl_vector<T> extract (unsigned int len, unsigned int start=0) const
    {
        assert( start < n && start + len <= n );
        return vnl_vector<T>( this->data() + start, len );
    }
    
    //: Replaces elements with index beginning at start, by values of v. O(n).
    vnl_vector_fixed& update(vnl_vector<T> const&, unsigned int start=0);
    
    //: Return sum of squares of elements
    abs_t squared_magnitude() const { return this->squaredNorm(); }
    
    //: Return magnitude (length) of vector
    abs_t magnitude() const { return this->norm(); }
    
    //: Return sum of absolute values of the elements
    abs_t one_norm() const { return this->template lpNorm<1>(); }
    
    //: Return sqrt of sum of squares of values of elements
    abs_t two_norm() const { return this->norm(); }
    
    //: Return largest absolute element value
    abs_t inf_norm() const { return this->template lpNorm<Eigen::Infinity>(); }
    
    //: Normalise by dividing through by the magnitude
    vnl_vector_fixed<T,n>& normalize() {
        base_class::normalize();
        return *this;
    }
    
    // These next 6 functions are should really be helper functions since they aren't
    // really proper functions on a vector in a philosophical sense.
    
    //: Root Mean Squares of values
    abs_t rms     () const { return std::sqrt((abs_t)base_class::squaredNorm()/size());; }
    
    //: Smallest value
    T min_value () const;
    
    //: Largest value
    T max_value () const;
    
    //: Location of smallest value
    unsigned arg_min() const;
    
    //: Location of largest value
    unsigned arg_max() const;
    
    //: Mean of values in vector
    T mean() const;
    
    //: Sum of values in a vector
    T sum() const;
    
    //: Reverse the order of the elements
    //  Element i swaps with element size()-1-i
    vnl_vector_fixed& flip();
    
    //: Return true if it's finite
    bool is_finite() const;
    
    //: Return true iff all the entries are zero.
    bool is_zero() const;
    
    //: Return true iff the size is zero.
    bool empty() const { return n==0; }
    
    //: Return true if *this == v
    bool operator_eq (vnl_vector_fixed<T,n> const& v) const;
    //: Return true if *this == v
    bool operator_eq (vnl_vector<T> const& v) const;
    
    // help functions
    inline static void add( const T* a, const T* b, T* r )
    {
        for ( unsigned int i=0; i < n; ++i,++r,++a,++b )
            *r = *a + *b;
    }
    
    inline static void add( const T* a, T b, T* r )
    {
        for ( unsigned int i=0; i < n; ++i,++r,++a )
            *r = *a + b;
    }
    
    inline static void sub( const T* a, const T* b, T* r )
    {
        for ( unsigned int i=0; i < n; ++i,++r,++a,++b )
            *r = *a - *b;
    }
    
    inline static void sub( const T* a, T b, T* r )
    {
        for ( unsigned int i=0; i < n; ++i,++r,++a )
            *r = *a - b;
    }
    
    inline static void sub( T a, const T* b, T* r )
    {
        for ( unsigned int i=0; i < n; ++i,++r,++b )
            *r = a - *b;
    }
    
    inline static void mul( const T* a, const T* b, T* r )
    {
        for ( unsigned int i=0; i < n; ++i,++r,++a,++b )
            *r = *a * *b;
    }
    
    inline static void mul( const T* a, T b, T* r )
    {
        for ( unsigned int i=0; i < n; ++i,++r,++a )
            *r = *a * b;
    }
    
    inline static void div( const T* a, const T* b, T* r )
    {
        for ( unsigned int i=0; i < n; ++i,++r,++a,++b )
            *r = *a / *b;
    }
    
    inline static void div( const T* a, T b, T* r )
    {
        for ( unsigned int i=0; i < n; ++i,++r,++a )
            *r = *a / b;
    }
};

template<class T, unsigned int n>
vnl_vector_fixed<T,n>::vnl_vector_fixed( const T& x0, const T& x1 )
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
template<class T, unsigned int n>
vnl_vector_fixed<T,n>::vnl_vector_fixed( const T& x0, const T& x1, const T& x2 )
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
template<class T, unsigned int n>
vnl_vector_fixed<T,n>::vnl_vector_fixed( const T& x0, const T& x1, const T& x2, const T& x3 )
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

template<class T, unsigned int n>
vnl_vector_fixed<T,n>
vnl_vector_fixed<T,n>::apply( T (*f)(T) )
{
    vnl_vector_fixed<T,n> ret;
    for ( unsigned int i = 0; i < n; ++i )
        ret[i] = f( this->data()[i] );
    return ret;
}

template<class T, unsigned int n>
vnl_vector_fixed<T,n>
vnl_vector_fixed<T,n>::apply( T (*f)(const T&) )
{
    vnl_vector_fixed<T,n> ret;
    for ( unsigned int i = 0; i < n; ++i )
        ret[i] = f( this->data()[i] );
    return ret;
}

template<class T, unsigned int n>
vnl_vector_fixed<T,n>&
vnl_vector_fixed<T,n>::update( const vnl_vector<T>& v, unsigned int start )
{
    unsigned stop = start + v.size();
    assert( stop <= n );
    for (unsigned i = start; i < stop; i++)
        (*this)(i) = v[i-start];
    return *this;
}

//: Smallest value
template<class T, unsigned int n>
T vnl_vector_fixed<T,n>::min_value() const
{
    return this->minCoeff();
}

//: Largest value
template<class T, unsigned int n>
T vnl_vector_fixed<T,n>::max_value() const
{
    return this->maxCoeff();
}

//: Location of smallest value
template<class T, unsigned int n>
unsigned vnl_vector_fixed<T,n>::arg_min() const
{
    unsigned int index = 0;
    this->minCoeff(&index);
    return index;
}

//: Location of largest value
template<class T, unsigned int n>
unsigned vnl_vector_fixed<T,n>::arg_max() const
{
    unsigned int index = 0;
    this->maxCoeff(&index);
    return index;
}

//: Mean of values in vector
template<class T, unsigned int n>
T vnl_vector_fixed<T,n>::mean() const
{
    if(this->size() == 0) return T{0};
    T sum = T{0};
    for(int i = 0; i<this->size(); ++i) {
        sum += (*this)[i];
    }
    return sum/this->size();
}

//: Sum of values in a vector
template<class T, unsigned int n>
T vnl_vector_fixed<T,n>::sum() const
{
    T sum = T{0};
    for(int i = 0; i<this->size(); ++i) {
        sum += (*this)[i];
    }
    return sum;
}

template <class T, unsigned int n>
vnl_vector_fixed<T,n>&
vnl_vector_fixed<T,n>::flip()
{
    for ( unsigned int i=0; 2*i+1 < n; ++i )
        std::swap( this->data()[i], this->data()[n-1-i] );
    return *this;
}

template <class T, unsigned int n>
bool
vnl_vector_fixed<T,n>::is_finite() const
{
    for ( unsigned i = 0; i < this->size(); ++i )
        if ( !vnl_math::isfinite( (*this)[i] ) )
            return false;
    
    return true;
}


template <class T, unsigned int n>
bool
vnl_vector_fixed<T,n>::is_zero() const
{
    T const zero(0);
    for ( unsigned i = 0; i < this->size(); ++i )
        if ( !( (*this)[i] == zero) )
            return false;
    
    return true;
}

template <class T, unsigned int n>
bool vnl_vector_fixed<T,n>::operator_eq (vnl_vector_fixed<T,n> const& v) const
{
    for ( unsigned i = 0; i < n; ++i )
        if ( (*this)[i] != v[i] )
            return false;
    return true;
}

//: Return true if *this == v
template <class T, unsigned int n>
bool vnl_vector_fixed<T,n>::operator_eq (vnl_vector<T> const& v) const
{
    assert( v.size() == n );
    for ( unsigned i = 0; i < n; ++i )
        if ( (*this)[i] != v[i] )
            return false;
    return true;
}

// --- Vector-scalar operators ----------------------------------------

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> operator+( const vnl_vector_fixed<T,n>& v, T s )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::add( v.data_block(), s, r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> operator+( const T& s,
                                       const vnl_vector_fixed<T,n>& v )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::add( v.data_block(), s, r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> operator-( const vnl_vector_fixed<T,n>& v, T s )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::sub( v.data_block(), s, r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> operator-( const T& s,
                                       const vnl_vector_fixed<T,n>& v )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::sub( s, v.data_block(), r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> operator*( const vnl_vector_fixed<T,n>& v, T s )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::mul( v.data_block(), s, r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> operator*( const T& s,
                                       const vnl_vector_fixed<T,n>& v )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::mul( v.data_block(), s, r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> operator/( const vnl_vector_fixed<T,n>& v, T s )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::div( v.data_block(), s, r.data_block() );
    return r;
}

// --- Vector-vector operators ----------------------------------------
//
// Includes overloads for the common case of mixing a fixed with a
// non-fixed. Because the operators are templated, the fixed will not
// be automatically converted to a non-fixed-ref. These do it for you.

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> operator+( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::add( a.data_block(), b.data_block(), r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector<T> operator+( const vnl_vector_fixed<T,n>& a, const vnl_vector<T>& b )
{
    return a.as_vector() + b;
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector<T> operator+( const vnl_vector<T>& a, const vnl_vector_fixed<T,n>& b )
{
    return a + b.as_vector();
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> operator-( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::sub( a.data_block(), b.data_block(), r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector<T> operator-( const vnl_vector_fixed<T,n>& a, const vnl_vector<T>& b )
{
    return a.as_vector() - b;
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector<T> operator-( const vnl_vector<T>& a, const vnl_vector_fixed<T,n>& b )
{
    return a - b.as_vector();
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> element_product( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::mul( a.data_block(), b.data_block(), r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector<T> element_product( const vnl_vector_fixed<T,n>& a, const vnl_vector<T>& b )
{
    assert( b.size() == n );
    vnl_vector<T> r(n);
    vnl_vector_fixed<T,n>::mul( a.data_block(), b.data_block(), r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector<T> element_product( const vnl_vector<T>& a, const vnl_vector_fixed<T,n>& b )
{
    assert( a.size() == n );
    vnl_vector<T> r(n);
    vnl_vector_fixed<T,n>::mul( a.data_block(), b.data_block(), r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector_fixed<T,n> element_quotient( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    vnl_vector_fixed<T,n> r;
    vnl_vector_fixed<T,n>::div( a.data_block(), b.data_block(), r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector<T> element_quotient( const vnl_vector_fixed<T,n>& a, const vnl_vector<T>& b )
{
    assert( b.size() == n );
    vnl_vector<T> r(n);
    vnl_vector_fixed<T,n>::div( a.data_block(), b.data_block(), r.data_block() );
    return r;
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline vnl_vector<T> element_quotient( const vnl_vector<T>& a, const vnl_vector_fixed<T,n>& b )
{
    assert( a.size() == n );
    vnl_vector<T> r(n);
    vnl_vector_fixed<T,n>::div( a.data_block(), b.data_block(), r.data_block() );
    return r;
}


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
template<class T, unsigned int n>
inline vnl_matrix<T> outer_product( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    return outer_product( a.as_vector(), b.as_vector());
}

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

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline bool operator==( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    return a.operator_eq(b);
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline bool operator==( vnl_vector_fixed<T,n> const& a, vnl_vector<T> const& b )
{
    return a.operator_eq(b);
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline bool operator==( vnl_vector<T> const& a, vnl_vector_fixed<T,n> const& b )
{
    return b.operator_eq(a);
}

//:
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline bool operator!=( const vnl_vector_fixed<T,n>& a, const vnl_vector_fixed<T,n>& b )
{
    return ! a.operator_eq(b);
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline bool operator!=( vnl_vector_fixed<T,n> const& a, vnl_vector<T> const& b )
{
    return ! a.operator_eq(b);
}

//:
// \relatesalso vnl_vector
// \relatesalso vnl_vector_fixed
template<class T, unsigned int n>
inline bool operator!=( vnl_vector<T> const& a, vnl_vector_fixed<T,n> const& b )
{
    return ! b.operator_eq(a);
}


#endif /* vnl_vector_fixed_h */
