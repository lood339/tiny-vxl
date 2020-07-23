//
//  vnl_vector.h
//  vxl
//
//  Created by Jimmy Chen on 2020-06-14.
//  Copyright © 2020 Nowhere Planet. All rights reserved.
//

#ifndef vnl_vector_h
#define vnl_vector_h

#include <Eigen/Dense>
#include <vnl/vnl_numeric_traits.h>
#include <vnl/vnl_error.h>
#include <vnl/vnl_math.h>

template <typename T> class vnl_matrix;
template <typename T> class vnl_vector;

template <typename T>
class vnl_vector: public Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>
{
    using base_class = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;
    using abs_t = typename vnl_numeric_traits<T>::abs_t;
public:
    vnl_vector()=default;
    
    //: Creates a vector containing n uninitialized elements.
    explicit vnl_vector(size_t len):base_class(len){}
   
    //: Creates a vector containing n elements, all set to v0.
    vnl_vector(size_t len, T const& v0):base_class(len) {
        this->setConstant(v0);
    }
    
    // Creates a vector of specified length and initialize first n elements with values. O(n).
    vnl_vector(size_t len, size_t n, const T* data_block):base_class(len){
        n = std::min(len, n);
        std::copy(data_block, data_block+n, this->data());
    }
    
    //: Construct a fixed-n-vector initialized from \a datablck
    //  The data \e must have enough data. No checks performed.
    explicit vnl_vector( const T* datablck, size_t n ){
        *this = base_class::Zero(n);
        std::memcpy(this->data(), datablck, n*sizeof(T));
    }
    
    //: Copy constructor.
    //vnl_vector(vnl_vector<T> const&);
    
    // NOTE: move-assignment must be allowed to throw an exception, because we need to maintain
    //       backwards compatibility and the move-construction & move-aasignment
    //       operators fall back to the copy-assignment operator behavior in
    //       cases when the memory is externally managed.
    //: Move-constructor.
    //vnl_vector(vnl_vector<T> &&);
    //: Move-assignment operator
    //vnl_vector<T>& operator=(vnl_vector<T>&& rhs);
    
    //: Destructor
    /** This destructor *must* be virtual to ensure that the vnl_vector_ref subclass destructor
     * is called and memory is not accidently de-allocated. */
    //virtual ~vnl_vector();
    
    template<typename OtherDrived>
    vnl_vector(const Eigen::MatrixBase<OtherDrived>& other):
    Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>(other){}
   
    
    template<typename OtherDrived>
    vnl_vector & operator=(const Eigen::MatrixBase<OtherDrived>& other){
        this->Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>::operator=(other);
        return *this;
    }
    
    //: Return the length, number of elements, dimension of this vector.
    //size_t size() const { return this->num_elmts; }
    
    //: Put value at given position in vector.
    inline void put(size_t i, T const& v)
    {
        assert(i<this->size());
        this->data()[i] = v;
    }
    
    //: Get value at element i
    inline T get(size_t  i) const
    {
        assert(i<this->size());
        return this->data()[i];
    }
    
    //: Set all values to v
    vnl_vector& fill(T const& v)
    {
        this->setConstant(v);
        return *this;
    }
    
    //: Sets elements to ptr[i]
    //  Note: ptr[i] must be valid for i=0..size()-1
    vnl_vector& copy_in(T const * ptr);
    
    //: Copy elements to ptr[i]
    //  Note: ptr[i] must be valid for i=0..size()-1
    void copy_out(T *) const; // from vector to array[].
    
    //: Sets elements to ptr[i]
    //  Note: ptr[i] must be valid for i=0..size()-1
    vnl_vector& set(T const *ptr) { return copy_in(ptr); }
    
    
    //: Return reference to the element at specified index.
    // There are assert style boundary checks - #define NDEBUG to turn them off.
    T       & operator()(size_t i)
    {
        assert(i<this->size());
        return this->data()[i];
    }
    //: Return reference to the element at specified index. No range checking.
    // There are assert style boundary checks - #define NDEBUG to turn them off.
    T const & operator()(size_t i) const
    {
        assert(i<this->size());
        return this->data()[i];
    }
    
    //: Return reference to the element at specified index. No range checking.
    T       & operator[](size_t i) { return this->data()[i]; }
    //: Return reference to the element at specified index. No range checking.
    T const & operator[](size_t i) const { return this->data()[i]; }
    
    //: Set all elements to value v
    vnl_vector<T>& operator=(T const&v) { fill(v); return *this; }
    
    //: Copy operator
    //vnl_vector<T>& operator=(vnl_vector<T> const& rhs);
    
    //: Add scalar value to all elements
    vnl_vector<T>& operator+=(T value);
    
    //: Subtract scalar value from all elements
    vnl_vector<T>& operator-=(T value);
    
    //: Multiply all elements by scalar
    vnl_vector<T>& operator*=(T );
    
    //: Divide all elements by scalar
    vnl_vector<T>& operator/=(T );
    
    //: Add rhs to this and return *this
    vnl_vector<T>& operator+=(vnl_vector<T> const& rhs);
    
    //: Subtract rhs from this and return *this
    vnl_vector<T>& operator-=(vnl_vector<T> const& rhs);
    
    //: *this = M*(*this) where M is a suitable matrix.
    //  this is treated as a column vector
    vnl_vector<T>& pre_multiply(vnl_matrix<T> const& M);
    
    //: *this = (*this)*M where M is a suitable matrix.
    //  this is treated as a row vector
    vnl_vector<T>& post_multiply(vnl_matrix<T> const& M);
    
    //: *this = (*this)*M where M is a suitable matrix.
    //  this is treated as a row vector
    vnl_vector<T>& operator*=(vnl_matrix<T> const& m) { return this->post_multiply(m); }
    
    //: Unary plus operator
    // Return new vector = (*this)
    vnl_vector<T> operator+() const { return *this; }
    
    //: Unary minus operator
    // Return new vector = -1*(*this)
    vnl_vector<T> operator-() const;
    
    vnl_vector<T> operator+(T v) const;
    
    vnl_vector<T> operator-(T v) const;
    
    vnl_vector<T> operator*(T v) const;
    
    vnl_vector<T> operator/(T v) const;
    
    vnl_vector<T> operator+(vnl_vector<T> const &v) const;
    
    vnl_vector<T> operator-(vnl_vector<T> const &v) const;
    
    // as a row vectror
    vnl_vector<T> operator*(vnl_matrix<T> const &M) const;
    //--------------------------------------------------------------------------------
    
    //: Access the contiguous block storing the elements in the vector. O(1).
    //  data_block()[0] is the first element of the vector
    T const* data_block() const { return this->data(); }
    
    //: Access the contiguous block storing the elements in the vector. O(1).
    //  data_block()[0] is the first element of the vector
    T      * data_block() { return this->data(); }
    
    //: Type defs for iterators
    typedef T       *iterator;
    //: Iterator pointing to start of data
    iterator begin() { return this->data(); }
    
    //: Iterator pointing to element beyond end of data
    iterator end() { return this->data()+this->size(); }
    
    //: Const iterator type
    typedef T const *const_iterator;
    //: Iterator pointing to start of data
    const_iterator begin() const { return this->data(); }
    //: Iterator pointing to element beyond end of data
    const_iterator end() const { return this->data()+this->size(); }
    
    //: Applies function to elements
    vnl_vector<T> apply(T (*f)(T)) const;
    //: Applies function to elements
    vnl_vector<T> apply(T (*f)(T const&)) const;
    
    //: Returns a subvector specified by the start index and length. O(n).
    vnl_vector<T> extract(size_t len, size_t start=0) const;
    
    //: Replaces elements with index beginning at start, by values of v. O(n).
    vnl_vector<T>& update(vnl_vector<T> const&, size_t start=0);
    
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
    vnl_vector<T>& normalize() {
        base_class::normalize();
        return *this;
    }
    
    // These next 6 functions are should really be helper functions since they aren't
    // really proper functions on a vector in a philosophical sense.
    
    //: Root Mean Squares of values
    abs_t rms() const { return std::sqrt(base_class::squaredNorm()/this->size()); }
    
    //: Smallest value
    T min_value() const;
    
    //: Largest value
    T max_value() const;
    
    //: Location of smallest value
    size_t arg_min() const;
    
    //: Location of largest value
    size_t arg_max() const;
    
    //: Mean of values in vector
    T mean() const {return base_class::mean();}
    
    //: Sum of values in a vector
    T sum() const {return base_class::sum();}
    
    //: Reverse the order of the elements
    //  Element i swaps with element size()-1-i
    vnl_vector<T>& flip();
    
    //: Reverse the order of the elements from index b to 1-e, inclusive.
    //  When b = 0 and e = size(), this is equivalent to flip();
    vnl_vector<T>& flip(const size_t &b, const size_t &e);
    
    //: Roll the vector forward by the specified shift.
    //  The shift is cyclical, such that the elements which
    //  are displaced from the end reappear at the beginning.
    //  Negative shifts and shifts >= the length of the array are supported.
    //  A new vector is returned; the underlying data is unchanged.
    vnl_vector<T> roll(const int &shift) const;
    
    //: Roll the vector forward by the specified shift.
    //  The shift is cyclical, such that the elements which
    //  are displaced from the end reappear at the beginning.
    //  Negative shifts and shifts >= the length of the array are supported.
    //
    vnl_vector& roll_inplace(const int &shift);
    
    //: Set this to that and that to this
    //void swap(vnl_vector<T> & that) noexcept;
    
    //: Return true if it's finite
    bool is_finite() const;
    
    //: Return true iff all the entries are zero.
    bool is_zero() const;
    
    //: Return true iff the size is zero.
    bool empty() const { return base_class::size() == 0; }
    
    //:  Return true if all elements of vectors are equal, within given tolerance
    bool is_equal(vnl_vector<T> const& rhs, double tol) const;
    
    //: Return true if *this == v
    bool operator_eq(vnl_vector<T> const& v) const;
    
    //: Equality test
    bool operator==(vnl_vector<T> const &that) const { return  this->operator_eq(that); }
    
    //: Inequality test
    bool operator!=(vnl_vector<T> const &that) const { return !this->operator_eq(that); }
    
    //: Resize to n elements.
    // This is a destructive resize, in that the old data is lost if size() != \a n before the call.
    // If size() is already \a n, this is a null operation.
    bool set_size(size_t n);
    
    //: Make the vector as if it had been default-constructed.
    void clear();
    
};

template<class T>
vnl_vector<T>&
vnl_vector<T>::copy_in (T const *ptr)
{
    std::copy( ptr, ptr + this->size(), this->data() );
    return *this;
}

//: Sets elements of an array to those in vector. O(n).

template<class T>
void vnl_vector<T>::copy_out (T *ptr) const
{
    std::copy( this->data(), this->data() + this->size(), ptr );
}

//: Add scalar value to all elements
template<typename T>
vnl_vector<T>& vnl_vector<T>::operator+=(T value)
{
    this->array() += value;
    return *this;
}

//: Subtract scalar value from all elements
template<typename T>
vnl_vector<T>& vnl_vector<T>::operator-=(T value)
{
    this->array() -= value;
    return *this;
}

//: Multiply all elements by scalar
template<typename T>
vnl_vector<T>& vnl_vector<T>::operator*=(T value)
{
    this->array() *= value;
    return *this;
}

//: Divide all elements by scalar
template<typename T>
vnl_vector<T>& vnl_vector<T>::operator/=(T value)
{
    this->array() /= value;
    return *this;
}

//: Add rhs to this and return *this
template<typename T>
vnl_vector<T>& vnl_vector<T>::operator+=(vnl_vector<T> const& rhs)
{
    if (this->size() != rhs.size()) {
        vnl_error_vector_dimension("operator+=", this->size(), rhs.size());
    }
    const unsigned int n = this->size();
    T *a = this->data();
    T const *b = rhs.data();
    for(unsigned int i = 0; i<n; ++i) {
        a[i] = T(a[i] + b[i]);
    }
    return *this;
}

//: Subtract rhs from this and return *this
template<typename T>
vnl_vector<T>& vnl_vector<T>::operator-=(vnl_vector<T> const& rhs)
{
    if (this->size() != rhs.size()) {
        vnl_error_vector_dimension("operator+=", this->size(), rhs.size());
    }
    const unsigned int n = this->size();
    T *a = this->data();
    T const *b = rhs.data();
    for(unsigned int i = 0; i<n; ++i) {
        a[i] = T(a[i] - b[i]);
    }
    return *this;
}

//: *this = M*(*this) where M is a suitable matrix.
//  this is treated as a column vector
template<typename T>
vnl_vector<T>& vnl_vector<T>::pre_multiply(vnl_matrix<T> const& m)
{
    if(m.cols() != this->size()) {
        vnl_error_vector_dimension("pre_multiply", m.cols(), this->size());
    }
    
    vnl_vector<T> v_copy = *this;
    *this = base_class::Zero(m.rows());
    for(int i = 0; i<m.rows(); ++i) {
        T s = T{0};
        for(int j = 0; j<m.cols(); ++j) {
            s += m(i,j) * v_copy[j];
        }
        (*this)[i] = s;
    }
    return *this;
}

//: *this = (*this)*M where M is a suitable matrix.
//  this is treated as a row vector
template<typename T>
vnl_vector<T>& vnl_vector<T>::post_multiply(vnl_matrix<T> const& m)
{
    if(m.rows() != this->size()) {
        vnl_error_vector_dimension("post_multiply", m.rows(), this->size());
    }
    
    vnl_vector<T> v_copy = *this;
    *this = base_class::Zero(m.cols());
    for(int j = 0; j<m.cols(); ++j) {
        T s = T{0};
        for(int i = 0; i<m.rows(); ++i) {
            s += v_copy[i] * m(i, j);
        }
        (*this)[j] = s;
    }
    return *this;
}

//: Creates new vector containing the negation of THIS vector. O(n).

template<typename T>
vnl_vector<T> vnl_vector<T>::operator- () const
{
    vnl_vector<T> result(this->size());
    for (size_t i = 0; i < this->size(); i++)
        result[i] = - (*this)[i];           // negate element
    return result;
}

template<typename T>
vnl_vector<T> vnl_vector<T>::operator+(T v) const {
    vnl_vector<T> result(this->size());
    std::transform(this->begin(), this->end(), result.begin(),
                   [v](T d) -> T { return d + v; });
    return result;
}

template<typename T>
vnl_vector<T> vnl_vector<T>::operator-(T v) const {
    vnl_vector<T> result(this->size());
    std::transform(this->begin(), this->end(), result.begin(),
                   [v](T d) -> T { return d - v; });
    return result;
}

template<typename T>
vnl_vector<T> vnl_vector<T>::operator*(T v) const {
    vnl_vector<T> result(this->size());
    std::transform(this->begin(), this->end(), result.begin(),
                   [v](T d) -> T { return d * v; });
    return result;
}

template<typename T>
vnl_vector<T> vnl_vector<T>::operator/(T v) const {
    vnl_vector<T> result(this->size());
    std::transform(this->begin(), this->end(), result.begin(),
                   [v](T d) -> T { return d / v; });
    return result;
}

template<typename T>
vnl_vector<T> vnl_vector<T>::operator+(vnl_vector<T> const &v) const {
    assert(this->size() == v.size());
    vnl_vector<T> result(this->size());
    for(int i = 0; i<this->size(); ++i) {
        result[i] = (*this)[i] + v[i];
    }
    return result;
}

template<typename T>
vnl_vector<T> vnl_vector<T>::operator-(vnl_vector<T> const &v) const {
    assert(this->size() == v.size());
    vnl_vector<T> result(this->size());
    for(int i = 0; i<this->size(); ++i) {
        result[i] = (*this)[i] - v[i];
    }
    return result;
}

// as a row vectror
template<typename T>
vnl_vector<T> vnl_vector<T>::operator*(vnl_matrix<T> const &M) const {
    
    vnl_vector<T> result(M.cols());
    if (this->size() != M.rows())
        vnl_error_vector_dimension("vnl_vector<>::operator*(M)", this->size(),
                                   M.rows());
    assert(this->size() == M.rows());
    
    for(int i = 0; i<result.size(); ++i) {
        T v = T{0};
        for(int j = 0; j<this->size(); ++j) {
            v += (*this)(j) * M(j, i);
        }
        result[i] = v;
    }
    return result;
}

//: Returns a subvector specified by the start index and length. O(n).
template<typename T>
vnl_vector<T> vnl_vector<T>::extract(size_t len, size_t start) const
{
    assert(start + len <= this->size());
    vnl_vector<T> result(len);
    for(int i = 0; i<len; ++i) {
        result[i] = (*this)[i+start];
    }
    return result;
}

template<typename T>
vnl_vector<T>& vnl_vector<T>::update (vnl_vector<T> const& v, size_t start)
{
    size_t stop = start + v.size();
    assert(stop <= this->size());
    for (size_t i = start; i < stop; i++)
        (*this)[i] = v[i-start];
    return *this;
}

//: Smallest value
template<typename T>
T vnl_vector<T>::min_value() const
{
    return this->minCoeff();
}

//: Largest value
template<typename T>
T vnl_vector<T>::max_value() const
{
    return this->maxCoeff();
}

//: Location of smallest value
template<typename T>
size_t vnl_vector<T>::arg_min() const
{
    unsigned int index = 0;
    this->minCoeff(&index);
    return index;
}

//: Location of largest value
template<typename T>
size_t vnl_vector<T>::arg_max() const
{
    unsigned int index = 0;
    this->maxCoeff(&index);
    return index;
}

//: Reverse the order of the elements
//  Element i swaps with element size()-1-i
template<typename T>
vnl_vector<T>& vnl_vector<T>::flip()
{
    size_t b = 0;
    size_t e = this->size() - 1;
    while(b < e) {
        T temp = (*this)[b];
        (*this)[b] = (*this)[e];
        (*this)[e] = temp;
        b++;
        e--;
    }
    return *this;
}

//: Reverse the order of the elements from index b to 1-e, inclusive.
//  When b = 0 and e = size(), this is equivalent to flip();
template<typename T>
vnl_vector<T>& vnl_vector<T>::flip(const size_t &b, const size_t &e)
{
    assert (!(b > this->size() || e > this->size() || b > e));
    for (size_t i=b;i<(e-b)/2+b;++i) {
        T tmp=(*this)[i];
        const size_t endIndex = e - 1 - (i-b);
        (*this)[i]=(*this)[endIndex];
        (*this)[endIndex]=tmp;
    }
    return *this;
    
}

//: Roll the vector forward by the specified shift.
//  The shift is cyclical, such that the elements which
//  are displaced from the end reappear at the beginning.
//  Negative shifts and shifts >= the length of the array are supported.
//  A new vector is returned; the underlying data is unchanged.
template<typename T>
vnl_vector<T> vnl_vector<T>::roll(const int &shift) const
{
    vnl_vector<T> v(this->size());
    const size_t wrapped_shift = shift % this->size();
    for (size_t i = 0; i < this->size(); ++i)
    {
        v[(i + wrapped_shift)%this->size()] = this->data()[i];
    }
    return v;
}

//: Roll the vector forward by the specified shift.
//  The shift is cyclical, such that the elements which
//  are displaced from the end reappear at the beginning.
//  Negative shifts and shifts >= the length of the array are supported.
template<typename T>
vnl_vector<T>& vnl_vector<T>::roll_inplace(const int &shift)
{
    const size_t wrapped_shift = (shift % (int)this->size() + this->size())%this->size();
    if (0 == wrapped_shift)
        return *this;
    return this->flip().flip(0,wrapped_shift).flip(wrapped_shift,this->size());
}

template <class T>
bool vnl_vector<T>::is_finite() const
{
    for (size_t i = 0; i < this->size();++i)
        if (!vnl_math::isfinite( (*this)[i] ))
            return false;
    
    return true;
}

template <class T>
bool vnl_vector<T>::is_zero() const
{
    T const zero(0);
    for (size_t i = 0; i < this->size();++i)
        if ( !( (*this)[i] == zero) )
            return false;
    
    return true;
}

template<class T>
bool vnl_vector<T>::operator_eq (vnl_vector<T> const& rhs) const
{
    if (this == &rhs)                               // same object => equal.
        return true;
    
    if (this->size() != rhs.size())                 // Size different ?
        return false;                                 // Then not equal.
    for (size_t i = 0; i < this->size(); i++)           // For each index
        if (!((*this)(i) == rhs[i]))          // Element different ?
            return false;                               // Then not equal.
    
    return true;                                    // Else same; return true.
}

//:  Return true if all elements of vectors are equal, within given tolerance
template <typename T>
bool vnl_vector<T>::is_equal(vnl_vector<T> const& rhs, double tol) const
{
    if (this == &rhs)                                         //Same object ? => equal.
        return true;
    
    if (this->size() != rhs.size())                           //Size different ?
        return false;
    for (size_t i = 0; i < this->size(); i++)
        if (std::abs(this->data()[i] - rhs.data()[i]) > tol)    //Element different ?
            return false;
    return true;
}

template<class T>
void vnl_vector<T>::clear()
{
    base_class::resize(0);
}

template<class T>
bool vnl_vector<T>::set_size(size_t n)
{
    if (!this->empty()) {
        // if no change in size, do not reallocate.
        if (this->size() == n)
            return false;
        base_class::resize(0);
    }
    else {
        // this happens if the vector is default constructed.
        base_class::resize(0);
    }
    return true;
}

template<typename T>
T dot_product (vnl_vector<T> const& v1, vnl_vector<T> const& v2)
{
    if (v1.size() != v2.size())
        vnl_error_vector_dimension ("dot_product", v1.size(), v2.size());

    return v1.dot(v2);
}

//: Hermitian inner product. O(n)
template<typename T>
T inner_product (vnl_vector<T> const& v1, vnl_vector<T> const& v2)
{
    if (v1.size() != v2.size())
        vnl_error_vector_dimension ("inner_product", v1.size(), v2.size());
    return v1.dot(v2.conjugate());
}

//: Returns the 'matrix element' <u|A|v> = u^t * A * v. O(mn).
template<typename T>
T bracket(vnl_vector<T> const &u, vnl_matrix<T> const &A, vnl_vector<T> const &v)
{
#ifndef NDEBUG
    if (u.size() != A.rows())
        vnl_error_vector_dimension("bracket",u.size(),A.rows());
    if (A.columns() != v.size())
        vnl_error_vector_dimension("bracket",A.columns(),v.size());
#endif
    T brak(0);
    for (size_t i=0; i<u.size(); ++i)
        for (size_t j=0; j<v.size(); ++j)
            brak += u[i]*A(i,j)*v[j];
    return brak;
}

// fsm : cos_angle should return a T, or a double-precision extension
// of T. "double" is wrong since it won't work if T is complex.
template <class T>
T cos_angle(vnl_vector<T> const& a, vnl_vector<T> const& b)
{
    typedef typename vnl_numeric_traits<T>::real_t real_t;
    typedef typename vnl_numeric_traits<T>::abs_t abs_t;
    typedef typename vnl_numeric_traits<abs_t>::real_t abs_r;
    
    real_t ab = inner_product(a,b);
    real_t a_b = static_cast<real_t>(
                                     std::sqrt( abs_r(a.squared_magnitude() * b.squared_magnitude()) ));
    return T( ab / a_b);
}

//: Returns smallest angle between two non-zero n-dimensional vectors. O(n).

template<class T>
double angle (vnl_vector<T> const& a, vnl_vector<T> const& b)
{
    typedef typename vnl_numeric_traits<T>::abs_t abs_t;
    typedef typename vnl_numeric_traits<abs_t>::real_t abs_r;
    const abs_r c = abs_r( cos_angle(a, b) );
    // IMS: sometimes cos_angle returns 1+eps, which can mess up std::acos.
    if (c >= 1.0) return 0;
    if (c <= -1.0) return vnl_math::pi;
    return std::acos( c );
}

//: Returns the nxn outer product of two nd-vectors, or [v1]^T*[v2]. O(n).

template<typename T>
vnl_matrix<T> outer_product (vnl_vector<T> const& v1,
                             vnl_vector<T> const& v2) {
    vnl_matrix<T> out(v1.size(), v2.size());
    for (size_t i = 0; i < out.rows(); i++)             // v1.column() * v2.row()
        for (size_t j = 0; j < out.columns(); j++)
            out(i, j) = v1[i] * v2[j];
    return out;
}

//: Returns new vector whose elements are the products v1[i]*v2[i]. O(n).
template<class T>
vnl_vector<T> element_product (vnl_vector<T> const& v1, vnl_vector<T> const& v2)
{
    assert(v1.size() == v2.size());

    vnl_vector<T> result(v1.size());
    for(int i = 0; i<v1.size(); ++i) {
        result[i] = v1[i] * v2[i];
    }
    
    return result;
}

template<class T>
vnl_vector<T> element_quotient (vnl_vector<T> const& v1, vnl_vector<T> const& v2)
{
    assert(v1.size() == v2.size());
    vnl_vector<T> result(v1.size());
    for(int i = 0; i<v1.size(); ++i) {
        result[i] = v1[i]/v2[i];
    }
    return result;
}

//: add scalar and vector. O(n).
// \relatesalso vnl_vector
template<typename T>
inline vnl_vector<T> operator+(T s, vnl_vector<T> const& v)
{
    return v.operator+(s);
}

//: subtract vector from scalar. O(n).
// \relatesalso vnl_vector
template<typename T>
inline vnl_vector<T> operator-(T s, vnl_vector<T> const& v)
{
    vnl_vector<T> result(v.size());
    for(size_t i=0; i< result.size(); ++i)
        result[i]= s - v[i];
    return result;
}

//: multiply scalar and vector. O(n).
// \relatesalso vnl_vector
template<typename T>
inline vnl_vector<T> operator*(T s, vnl_vector<T> const& v)
{
    return v*s;
}

//: Euclidean Distance between two vectors.
// Sum of Differences squared.
// \relatesalso vnl_vector
template<class T>
inline T vnl_vector_ssd(vnl_vector<T> const& v1, vnl_vector<T> const& v2)
{
#ifndef NDEBUG
    if (v1.size() != v2.size())
        vnl_error_vector_dimension("vnl_vector_ssd", v1.size(), v2.size());
#endif
    vnl_vector<T> dif = v1 - v2;
    return dif.squared_magnitude();
}


//: multiply matrix and (column) vector. O(m*n).
// \relatesalso vnl_vector
// \relatesalso vnl_matrix
template<typename T>
inline vnl_vector<T> operator*(vnl_matrix<T> const& M, vnl_vector<T> const& v)
{
    assert(M.cols() == v.size());
    vnl_vector<T> result(M.rows());
    for(int i = 0; i<result.size(); ++i) {
        T s = T{0};
        for(int j = 0; j<M.cols(); ++j) {
            s += M(i, j) * v(j);
        }
        result[i] = s;
    }
    //vnl_sse<T>::matrix_x_vector(M.begin(), v.begin(), result.begin(), M.rows(), M.cols());
    return result;
}





#endif /* vnl_vector_h */
