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

// fixed length column vector

template <typename T, unsigned int n> class vnl_vector_fixed;

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
};


#endif /* vnl_vector_fixed_h */
