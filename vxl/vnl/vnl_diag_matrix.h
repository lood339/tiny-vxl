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
    vnl_diag_matrix() = default;
    vnl_diag_matrix(const vnl_diag_matrix<T> &)  = default;
    vnl_diag_matrix(vnl_diag_matrix<T> &&)  = default;
    vnl_diag_matrix& operator=(const vnl_diag_matrix<T> &)  = default;
    vnl_diag_matrix& operator=(vnl_diag_matrix<T> &&)  = default;
    ~vnl_diag_matrix() = default;
   
    //: Construct an empty diagonal matrix.
    vnl_diag_matrix(unsigned nn) : base_class(nn) {}
    
    //: Construct a diagonal matrix with diagonal elements equal to value.
    vnl_diag_matrix(unsigned nn, T const& value) : base_class(nn) {
        for(int i = 0; i<base_class::diagonal().size(); ++i) {
            base_class::diagonal()[i] = value;
        }
    }
    
    //: Construct a diagonal matrix from a vnl_vector.
    //  The vector elements become the diagonal elements.
    vnl_diag_matrix(const vnl_vector<T>& v):base_class(v.size()){
        std::copy(v.data(), v.data()+v.size(), this->data());
    }
    
    template<unsigned int n>
    vnl_diag_matrix(const vnl_vector_fixed<T, n>& v) = delete;
    
    //: Construct a diagonal matrix from a vnl_vector.
    //  The vector elements become the diagonal elements.
    //vnl_diag_matrix(vnl_vector<T> const& that): diagonal_(that) {}
    
    
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
    
    // Operations----------------------------------------------------------------
    
    //: In-place arithmetic operation
    inline vnl_diag_matrix<T>& operator*=(T v) { this->array() *= v; return *this; }
    //: In-place arithmetic operation
    inline vnl_diag_matrix<T>& operator/=(T v) { this->array() /= v; return *this; }

    // Computations--------------------------------------------------------------
    
    vnl_diag_matrix& invert_in_place();
    T determinant() const;
    vnl_vector<T> solve(vnl_vector<T> const& b) const;
    void solve(vnl_vector<T> const& b, vnl_vector<T>* out) const;

    
    // Data Access---------------------------------------------------------------
    
    inline T operator () (unsigned i, unsigned j) const {
        return (i != j) ? T(0) : base_class::diagonal()[i];
    }
    
    inline T& operator () (unsigned i, unsigned j) {
        assert(i == j); (void)j;
        return base_class::diagonal()[i];
    }
    
    inline T& operator() (unsigned i) { return base_class::diagonal()[i]; }
    inline T const& operator() (unsigned i) const { return base_class::diagonal()[i]; }
    
    inline T& operator[] (unsigned i) { return base_class::diagonal()[i]; }
    inline T const& operator[] (unsigned i) const { return base_class::diagonal()[i]; }
    
    //: Return a vector (copy) with the content of the (main) diagonal
    inline vnl_vector<T> get_diagonal() const {
        vnl_vector<T> diag(this->size());
        for(int i =0 ; i<this->size(); ++i) {
            diag[i] = base_class::diagonal()[i];
        }
        return diag;
    }
    //: Return diagonal elements as a vector
    // oritinal return type: vnl_vector<T> const&
    vnl_vector<T> diagonal() const {
        vnl_vector<T> diag(this->size());
        for(int i =0 ; i<this->size(); ++i) {
            diag[i] = base_class::diagonal()[i];
        }
        return diag;
    }
    
    //: Set all diagonal elements of matrix to specified value.
    inline vnl_diag_matrix& fill_diagonal (T const& v) { this->setConstant(v); return *this; }
    
    //: Sets the diagonal elements of this matrix to the specified list of values.
    inline vnl_diag_matrix& set_diagonal(vnl_vector<T> const& v) {
        this->resize(v.size());
        std::copy(v.data(), v.data()+v.size(), this->data());
        return *this;
    }
    
    /*
    // iterators
    typedef typename vnl_vector<T>::iterator iterator;
    inline iterator begin() { return diagonal_.begin(); }
    inline iterator end() { return diagonal_.end(); }
    typedef typename vnl_vector<T>::const_iterator const_iterator;
    inline const_iterator begin() const { return diagonal_.begin(); }
    inline const_iterator end() const { return diagonal_.end(); }
     */
    
    
    //: Return the total number of elements stored by the matrix.
    // Since vnl_diag_matrix only stores values on the diagonal
    // and must be square, size() == rows() == cols().
    inline unsigned int size() const { return base_class::rows(); }
    
    //: Return the number of rows.
    inline unsigned int rows() const { return base_class::rows(); }
    
    //: Return the number of columns.
    // A synonym for columns().
    inline unsigned int cols() const { return base_class::rows(); }
    
    //: Return the number of columns.
    // A synonym for cols().
    inline unsigned int columns() const { return base_class::rows(); }
    
    //: set element with boundary checks.
    inline void put (unsigned r, unsigned c, T const& v) {
        assert(r == c);
        assert(r < this->rows());
        this->data[r] = v;
    }
    
    //: get element with boundary checks.
    inline T get (unsigned r, unsigned c) const {
        assert(r == c);
        assert(r < this->rows());
        return this->data[r];
    }
    
    /*
#if VXL_LEGACY_FUTURE_REMOVE
    VXL_DEPRECATED_MSG("Deprecated inconsistent return type.\nWARNING: .as_ref returns a vnl_matrix, not a vnl_matrix_ref, use .as_matrix() directly")
#endif
    vnl_matrix<T> as_ref() const { return as_matrix(); }
     */
    
    // Need this until we add a vnl_diag_matrix ctor to vnl_matrix;
    vnl_matrix<T> as_matrix() const;
    
    inline void set_size(int n) { this->resize(n); }
    
    inline void clear() { this->resize(0); }
    inline vnl_diag_matrix& fill(T const &x) {
        for(int i = 0; i<base_class::diagonal().size(); ++i) {
            base_class::diagonal()[i] = x;
        }
        return *this;
    }
    
    //: Return pointer to the diagonal elements as a contiguous 1D C array;
    inline T*       data_block()       { return base_class::diagonal().data(); }
    inline T const* data_block() const { return base_class::diagonal().data(); }
    
    
    //: Set diagonal elements using vector
    inline vnl_diag_matrix& set(vnl_vector<T> const& v)  {
        this->resize(v.size());
        std::copy(v.data(), v.data()+v.size(), this->data());
        return *this;
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

//: Convert a vnl_diag_matrix to a Matrix.
template <class T>
vnl_matrix<T> vnl_diag_matrix<T>::as_matrix() const
{
    unsigned len = this->size();
    vnl_matrix<T> ret(len, len);
    for (unsigned i = 0; i < len; ++i)
    {
        unsigned j;
        for (j = 0; j < i; ++j)
            ret(i,j) = T(0);
        for (j = i+1; j < len; ++j)
            ret(i,j) = T(0);
        ret(i,i) = base_class::diagonal()[i];
    }
    return ret;
}

//: Invert a vnl_diag_matrix in-situ.
// Just replaces each element with its reciprocal.
template <class T>
inline vnl_diag_matrix<T>& vnl_diag_matrix<T>::invert_in_place()
{
    unsigned len = this->size();
    T* d = this->data();
    T one = T(1);
    for (unsigned i = 0; i < len; ++i)
        d[i] = one / d[i];
    return *this;
}

//: Return determinant as product of diagonal values.
template <class T>
inline T vnl_diag_matrix<T>::determinant() const
{
    T det = T(1);
    T const* d = this->data();
    unsigned len = this->size();
    for (unsigned i = 0; i < len; ++i)
        det *= d[i];
    return det;
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
    assert(A.rows() == n); assert(A.cols() == n);
    vnl_matrix<T> ret(A);
    T const* d = D.data_block();
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
    T const* d = D.data_block();
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
    T const* d = D.data_block();
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
    assert(A.cols() == D.size());
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

//: Multiply a vnl_vector by a vnl_diag_matrix.  n flops.
// \relatesalso vnl_diag_matrix
// \relatesalso vnl_vector
template <class T>
inline vnl_vector<T> operator* (vnl_vector<T> const& A, vnl_diag_matrix<T> const& D)
{
    assert(A.size() == D.size());
    return element_product(D.diagonal(), A);
}


#endif /* vnl_diag_matrix_h */
