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
    
    //: Construct a fixed-n-vector initialized from \a datablck
    //  The data \e must have enough data. No checks performed.
    explicit vnl_vector( const T* datablck, size_t n )
    {
        *this = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>::Zero(n);
        std::memcpy(this->data(), datablck, n*sizeof(T));
    }
    
    // Creates a vector of specified length and initialize first n elements with values. O(n).
    vnl_vector(size_t len, size_t n, const T* data_block):base_class(len)
    {
        n = std::min(len, n);
        std::copy(data_block, data_block+n, this->data());
    }
    
    template<typename OtherDrived>
    vnl_vector(const Eigen::MatrixBase<OtherDrived>& other):
        Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>(other)
    {}
    
    template<typename OtherDrived>
    vnl_vector & operator=(const Eigen::MatrixBase<OtherDrived>& other)
    {
        this->Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>::operator=(other);
        return *this;
    }
    
    
    //:  Return true if all elements of vectors are equal, within given tolerance
    bool is_equal(vnl_vector<T> const& rhs, double tol) const
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
    
    //: Add scalar value to all elements
    vnl_vector<T>& operator+=(T );
    
    //: Subtract scalar value from all elements
    vnl_vector<T>& operator-=(T value) { return *this += T(-value); }
    
    //: Multiply all elements by scalar
    vnl_vector<T>& operator*=(T );
    
    //: Divide all elements by scalar
    vnl_vector<T>& operator/=(T );
    
    //: Add rhs to this and return *this
    vnl_vector<T>& operator+=(vnl_vector<T> const& rhs);
    
    //: Subtract rhs from this and return *this
    vnl_vector<T>& operator-=(vnl_vector<T> const& rhs);
    
    vnl_vector<T> operator+(vnl_vector<T> const &v) const {
        assert(this->size() == v.size());
        vnl_vector<T> result(this->size());
        for(int i = 0; i<this->size(); ++i) {
            result[i] = (*this)[i] + v[i];
        }
        return result;
    }
    
    vnl_vector<T> operator-(vnl_vector<T> const &v) const {
        assert(this->size() == v.size());
        vnl_vector<T> result(this->size());
        for(int i = 0; i<this->size(); ++i) {
            result[i] = (*this)[i] - v[i];
        }
        return result;
    }
    
    vnl_vector<T> operator*(vnl_matrix<T> const &M) const {
        
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
    //abs_t one_norm() const { return vnl_c_vector<T>::one_norm(begin(), size()); }
    
    //: Return sqrt of sum of squares of values of elements
    abs_t two_norm() const { return this->norm(); }
    
    //: Return largest absolute element value
    //abs_t inf_norm() const { return vnl_c_vector<T>::inf_norm(begin(), size()); }
    
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
    
    
};

template<class T>
vnl_vector<T>& vnl_vector<T>::update (vnl_vector<T> const& v, size_t start)
{
    size_t stop = start + v.size();
    assert(stop <= this->size());
    for (size_t i = start; i < stop; i++)
        (*this)[i] = v[i-start];
    return *this;
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
    
    //vnl_sse<T>::element_product(v1.begin(), v2.begin(), result.begin(), v1.size());
    
    return result;
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
