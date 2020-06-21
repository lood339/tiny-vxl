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

template <typename T> class vnl_matrix;
template <typename T> class vnl_vector;

template <typename T>
class vnl_vector: public Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>
{
    using base_class = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;
public:
    vnl_vector()=default;
    
    //: Creates a vector containing n uninitialized elements.
    explicit vnl_vector(size_t len):base_class(len)
    {
        
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
    
    
};

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

/*
//: multiply matrix and (column) vector. O(m*n).
// \relatesalso vnl_vector
// \relatesalso vnl_matrix
template<class T>
inline vnl_vector<T> operator*(vnl_matrix<T> const& M, vnl_vector<T> const& v)
{
    vnl_vector<T> result(M.rows());
#ifndef NDEBUG
    if (M.cols() != v.size())
        vnl_error_vector_dimension ("vnl_vector<>::operator*(M, v)", M.cols(), v.size());
#endif
    //vnl_sse<T>::matrix_x_vector(M.begin(), v.begin(), result.begin(), M.rows(), M.cols());
    return result;
}
 */




#endif /* vnl_vector_h */
