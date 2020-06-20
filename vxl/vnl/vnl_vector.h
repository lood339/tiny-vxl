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




#endif /* vnl_vector_h */
