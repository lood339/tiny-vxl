//
//  vnl_diag_matrix.h
//  vxl
//
//  Created by Jimmy Chen on 2020-06-20.
//  Copyright © 2020 Nowhere Planet. All rights reserved.
//

#ifndef vnl_diag_matrix_h
#define vnl_diag_matrix_h

#include <Eigen/Dense>

template<typename T> class vnl_matrix;
template<typename T> class vnl_vector;
template<typename T, unsigned int n> class vnl_vector_fixed;
template<typename T> class vnl_diag_matrix;

template<typename T>
class vnl_diag_matrix: public Eigen::DiagonalMatrix<T, Eigen::Dynamic, Eigen::Dynamic>
{
    using base_class = Eigen::DiagonalMatrix<T, Eigen::Dynamic, Eigen::Dynamic>;
public:
    vnl_diag_matrix()=default;
    
    //: Construct an empty diagonal matrix.
    vnl_diag_matrix(unsigned nn) : base_class(nn) {}
    
    //: Construct a diagonal matrix with diagonal elements equal to value.
    vnl_diag_matrix(unsigned nn, T const& value) : base_class(nn, nn) {
        this->setConstant(value);
    }
    
    template<unsigned int n>
    vnl_diag_matrix(const vnl_vector_fixed<T, n>& v) = delete;
    
    //: Construct a diagonal matrix from a vnl_vector.
    //  The vector elements become the diagonal elements.
    //vnl_diag_matrix(vnl_vector<T> const& that): diagonal_(that) {}
    ~vnl_diag_matrix() = default;
    
    // This constructor allows us to construct from Eigen expressions
    // https://eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html
    template<typename OtherDerived>
    vnl_diag_matrix(const Eigen::MatrixBase<OtherDerived>& other):
    base_class(other)
    {}
    
    // This method allows us to assign Eigen expressions
    template<typename OtherDerived>
    vnl_diag_matrix<T> & operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->base_class::operator=(other);
        return *this;
    }
    
    // Data Access---------------------------------------------------------------
    
    inline T operator () (unsigned i, unsigned j) const {
        return (i != j) ? T(0) : base_class::diagonal()[i];
    }
    
    inline T& operator () (unsigned i, unsigned j) {
        assert(i == j); (void)j;
        return base_class::diagonal()[i];
    }
    
    //: Return the total number of elements stored by the matrix.
    // Since vnl_diag_matrix only stores values on the diagonal
    // and must be square, size() == rows() == cols().
    inline unsigned int size() const { return base_class::rows(); }
    
    //: Return the number of rows.
    inline unsigned int rows() const { return base_class::rows(); }
    
    //: Return the number of columns.
    // A synonym for columns().
    inline unsigned int cols() const { return base_class::rows(); }
    
    
    //: Return diagonal elements as a vector
    inline vnl_vector<T> diagonal() const
    {
        vnl_vector<T> diag(this->size());
        for(int i =0 ; i<this->size(); ++i) {
            diag[i] = base_class::diagonal()[i];
        }
        return diag;
    }
    
    
};

template <typename T>
std::ostream& operator<< (std::ostream& s, vnl_diag_matrix<T> const& D)
{
    s << "diag([ ";
    for (unsigned i=0; i<D.rows(); ++i)
        s << D(i,i) << ' ';
    return s << "])";
}


//: Add two vnl_diag_matrices.  Just add the diag elements - n flops
// \relatesalso vnl_diag_matrix
template <class T>
inline vnl_diag_matrix<T> operator+ (vnl_diag_matrix<T> const& A, vnl_diag_matrix<T> const& B)
{
    assert(A.size() == B.size());
    vnl_diag_matrix<T> ret = A;
    std::cout<<"Debug A.size()"<<A.size()<<std::endl;
    for (unsigned i = 0; i < A.size(); ++i)
        ret(i,i) += B(i,i);
    return ret;
}

//: Add a vnl_diag_matrix to a vnl_matrix.  n adds, mn copies.
// \relatesalso vnl_diag_matrix
// \relatesalso vnl_matrix
template <class T>
inline vnl_matrix<T> operator+ (vnl_matrix<T> const& A, vnl_diag_matrix<T> const& D)
{
    const unsigned n = D.size();
    assert(A.rows() == n); assert(A.columns() == n);
    vnl_matrix<T> ret(A);
    T const* d = D.data();
    for (unsigned j = 0; j < n; ++j)
        ret(j,j) += d[j];
    return ret;
}

//: Add a vnl_matrix to a vnl_diag_matrix.  n adds, mn copies.
// \relatesalso vnl_diag_matrix
// \relatesalso vnl_matrix
template <class T>
inline vnl_matrix<T> operator+ (vnl_diag_matrix<T> const& D, vnl_matrix<T> const& A)
{
    return A + D;
}

//: Subtract two vnl_diag_matrices.  Just subtract the diag elements - n flops
// \relatesalso vnl_diag_matrix
template <class T>
inline vnl_diag_matrix<T> operator- (vnl_diag_matrix<T> const& A, vnl_diag_matrix<T> const& B)
{
    assert(A.size() == B.size());
    vnl_diag_matrix<T> ret = A;
    for (unsigned i = 0; i < A.size(); ++i)
        ret(i,i) -= B(i,i);
    return ret;
}

//: Subtract a vnl_diag_matrix from a vnl_matrix.  n adds, mn copies.
// \relatesalso vnl_diag_matrix
// \relatesalso vnl_matrix
template <class T>
inline vnl_matrix<T> operator- (vnl_matrix<T> const& A, vnl_diag_matrix<T> const& D)
{
    const unsigned n = D.size();
    assert(A.rows() == n); assert(A.cols() == n);
    vnl_matrix<T> ret(A);
    T const* d = D.data();
    for (unsigned j = 0; j < n; ++j)
        ret(j,j) -= d[j];
    return ret;
}

//: Subtract a vnl_matrix from a vnl_diag_matrix.  n adds, mn copies.
// \relatesalso vnl_diag_matrix
// \relatesalso vnl_matrix
template <class T>
inline vnl_matrix<T> operator- (vnl_diag_matrix<T> const& D, vnl_matrix<T> const& A)
{
    const unsigned n = D.size();
    assert(A.rows() == n); assert(A.cols() == n);
    vnl_matrix<T> ret(n, n);
    T const* d = D.data();
    for (unsigned i = 0; i < n; ++i)
    {
        for (unsigned j = 0; j < i; ++j)
            ret(i,j) = -A(i,j);
        for (unsigned j = i+1; j < n; ++j)
            ret(i,j) = -A(i,j);
        ret(i,i) = d[i] - A(i,i);
    }
    return ret;
}


//: Multiply two vnl_diag_matrices.  Just multiply the diag elements - n flops
// \relatesalso vnl_diag_matrix
template <class T>
inline vnl_diag_matrix<T> operator* (vnl_diag_matrix<T> const& A, vnl_diag_matrix<T> const& B)
{
    assert(A.size() == B.size());
    vnl_diag_matrix<T> ret = A;
    for (unsigned i = 0; i < A.size(); ++i)
        ret(i,i) *= B(i,i);
    return ret;
}


//: Multiply a vnl_matrix by a vnl_diag_matrix.  Just scales the columns - mn flops
// \relatesalso vnl_diag_matrix
// \relatesalso vnl_matrix
template <class T>
inline vnl_matrix<T> operator* (vnl_matrix<T> const& A, vnl_diag_matrix<T> const& D)
{
    assert(A.columns() == D.size());
    vnl_matrix<T> ret(A.rows(), A.cols());
    for (unsigned i = 0; i < A.rows(); ++i)
        for (unsigned j = 0; j < A.cols(); ++j)
            ret(i,j) = A(i,j) * D(j,j);
    return ret;
}

//: Multiply a vnl_diag_matrix by a vnl_matrix.  Just scales the rows - mn flops
// \relatesalso vnl_diag_matrix
// \relatesalso vnl_matrix
template <class T>
inline vnl_matrix<T> operator* (vnl_diag_matrix<T> const& D, vnl_matrix<T> const& A)
{
    assert(A.rows() == D.size());
    vnl_matrix<T> ret(A.rows(), A.cols());
    T const* d = D.data_block();
    for (unsigned i = 0; i < A.rows(); ++i)
        for (unsigned j = 0; j < A.cols(); ++j)
            ret(i,j) = A(i,j) * d[i];
    return ret;
}


//: Multiply a vnl_diag_matrix by a vnl_vector.  n flops.
// \relatesalso vnl_diag_matrix
// \relatesalso vnl_vector
template <class T>
inline vnl_vector<T> operator* (vnl_diag_matrix<T> const& D, vnl_vector<T> const& A)
{
    assert(A.size() == D.size());
    
    return element_product(D.diagonal(), A);
}


/*
//: Multiply a vnl_vector by a vnl_diag_matrix.  n flops.
// \relatesalso vnl_diag_matrix
// \relatesalso vnl_vector
template <class T>
inline vnl_vector<T> operator* (vnl_vector<T> const& A, vnl_diag_matrix<T> const& D)
{
    assert(A.size() == D.size());
    return element_product(D.diagonal(), A);
}
 */


#endif /* vnl_diag_matrix_h */
