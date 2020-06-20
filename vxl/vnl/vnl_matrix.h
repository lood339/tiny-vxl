//
//  vnl_matrix.h
//  vxl
//
//  Created by Jimmy Chen on 2020-06-14.
//  Copyright © 2020 Nowhere Planet. All rights reserved.
//

#ifndef vnl_matrix_h
#define vnl_matrix_h

#include <Eigen/Dense>

template <typename T> class vnl_matrix;
template <typename T> class vnl_vector;

template <typename T>
class vnl_matrix: public Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
    using base_class = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
public:
    vnl_matrix()=default;
    
    //: Construct a matrix of size r rows by c columns
    // Contents are unspecified.
    vnl_matrix(unsigned int r, unsigned int c):base_class(r, c)
    {
    }
    
    //: Construct a matrix of size r rows by c columns, and all elements equal to v0
    // Complexity $O(r.c)$
    vnl_matrix(unsigned int r, unsigned int c, T const& v0):base_class(r, c)
    {
        this->setConstant(v0);
    }
    
    //: Construct a matrix of size r rows by c columns, initialised by an automatic array
    // The first n elements, are initialised row-wise, to values.
    // Complexity $O(n)$
    vnl_matrix(unsigned int r, unsigned int c, unsigned int n, T const values[]):base_class(r, c) // use automatic arrays.
    {
        if (n > this->rows()*this->cols()) {
            n = static_cast<unsigned int>(this->rows()*this->cols());
        }
        std::copy(values, values + n, this->data());
    }
    
    //: Construct a matrix of size r rows by c columns, initialised by a memory block
    // The values are initialise row wise from the data.
    // Complexity $O(r.c)$
    vnl_matrix(T const* data_block, unsigned int r, unsigned int c):base_class(r, c)      // fill row-wise.
    {
        std::copy(data_block, data_block + r*c, this->data());
    }
    
    // This constructor allows us to construct from Eigen expressions
    // https://eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html
    template<typename OtherDerived>
    vnl_matrix(const Eigen::MatrixBase<OtherDerived>& other):
        base_class(other)
    {}
    
    // This method allows us to assign Eigen expressions
    template<typename OtherDerived>
    vnl_matrix<T> & operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->base_class::operator=(other);
        return *this;
    }
    
    //: Set all elements to value v
    // Complexity $O(r.c)$
    vnl_matrix<T>& operator=(T const&v)
    {
        this->setConstant(v);
        return *this;
    }
    
    // ----------------------- Arithmetic --------------------------------
    // note that these functions should not pass scalar as a const&.
    // Look what would happen to A /= A(0,0).
    
    //: Add rhs to each element of lhs matrix in situ
    vnl_matrix<T>& operator+=(T value)
    {
        this->array() += value;
        return *this;
    }
    
    //: Subtract rhs from each element of lhs matrix in situ
    vnl_matrix<T>& operator-=(T value)
    {
        this->array() -= value;
        return *this;
    }
    
    //: Scalar multiplication in situ of lhs matrix  by rhs
    vnl_matrix<T>& operator*=(T value)
    {
        this->array() *= value;
        return *this;
    }
    
    //: Scalar division of lhs matrix  in situ by rhs
    vnl_matrix<T>& operator/=(T value)
    {
        this->array() /= value;
        return *this;
    }
    
    //: Add rhs to lhs  matrix in situ
    vnl_matrix<T>& operator+=(vnl_matrix<T> const& rhs)
    {
        assert(rhs.rows() == this->rows());
        assert(rhs.cols() == this->cols());
        
        const unsigned int n = this->rows() * this->cols();
        T *a = this->data();
        T const *b = rhs.data();
        for(unsigned int i = 0; i<n; ++i) {
            a[i] = T(a[i] + b[i]);
        }
        return *this;
    }
    //: Subtract rhs from lhs matrix in situ
    vnl_matrix<T>& operator-=(vnl_matrix<T> const& rhs)
    {
        assert(rhs.rows() == this->rows());
        assert(rhs.cols() == this->cols());
        
        const unsigned int n = this->rows() * this->cols();
        T *a = this->data();
        T const *b = rhs.data();
        for(unsigned int i = 0; i<n; ++i) {
            a[i] = T(a[i] - b[i]);
        }
        return *this;
    }
    //: Multiply lhs matrix in situ by rhs
    vnl_matrix<T>& operator*=(vnl_matrix<T> const&rhs) { return *this = (*this) * rhs; }
    
    //: Negate all elements of matrix
    vnl_matrix<T> operator-() const
    {
        vnl_matrix<T> result(this->rows(), this->cols());
        for (unsigned int i = 0; i < this->rows(); i++)
            for (unsigned int j = 0; j < this->cols(); j++)
                result(i, j) = - (*this)(i, j);
        return result;
    }
    
    //: Add rhs to each element of lhs matrix and return result in new matrix
    vnl_matrix<T> operator+(T const& v) const {
        vnl_matrix<T> result(this->rows(), this->cols());
        const unsigned int n = this->rows() * this->cols();
        T const *m = this->data();
        T *dst = result.data();
        
        for (unsigned int i = 0; i < n; ++i)
            dst[i] = T(m[i] + v);
        return result;
    }
    
    //: Subtract rhs from each element of lhs matrix and return result in new matrix
    vnl_matrix<T> operator-(T const& v) const {
        vnl_matrix<T> result(this->rows(), this->cols());
        const unsigned int n = this->rows() * this->cols();
        T const *m = this->data();
        T *dst = result.data();
        
        for (unsigned int i = 0; i < n; ++i)
            dst[i] = T(m[i] - v);
        return result;
    }
    
    //: Scalar multiplication of lhs matrix by rhs  and return result in new matrix
    vnl_matrix<T> operator*(T const& v) const {
        vnl_matrix<T> result(this->rows(), this->cols());
        const unsigned int n = this->rows() * this->cols();
        T const *m = this->data();
        T *dst = result.data();
        
        for (unsigned int i = 0; i < n; ++i)
            dst[i] = T(m[i] * v);
        return result;
    }
    
    //: Scalar division of lhs matrix by rhs and return result in new matrix
    vnl_matrix<T> operator/(T const& v) const {
        vnl_matrix<T> result(this->rows(), this->cols());
        const unsigned int n = this->rows() * this->cols();
        T const *m = this->data();
        T *dst = result.data();
        
        for (unsigned int i = 0; i < n; ++i)
            dst[i] = T(m[i] / v);
        return result;
    }
    
    //: Matrix add rhs to lhs matrix and return result in new matrix
    vnl_matrix<T> operator+(vnl_matrix<T> const& rhs) const
    {
        assert(rhs.rows() == this->rows());
        assert(rhs.cols() == this->cols());
        vnl_matrix<T> result(this->rows(), this->cols());
        const unsigned int n = this->rows() * this->cols();
        T const *a = this->data();
        T const *b = rhs.data();
        T *dst = result.data();
        for(unsigned int i = 0; i<n; ++i) {
            dst[i] = T(a[i] + b[i]);
        }
        return result;
    }
    
    //: Matrix subtract rhs from lhs and return result in new matrix
    vnl_matrix<T> operator-(vnl_matrix<T> const& rhs) const
    {
        assert(rhs.rows() == this->rows());
        assert(rhs.cols() == this->cols());
        vnl_matrix<T> result(this->rows(), this->cols());
        const unsigned int n = this->rows() * this->cols();
        T const *a = this->data();
        T const *b = rhs.data();
        T *dst = result.data();
        for(unsigned int i = 0; i<n; ++i) {
            dst[i] = T(a[i] - b[i]);
        }
        return result;
    }
    
    //: Matrix multiply lhs by rhs matrix and return result in new matrix
    vnl_matrix<T> operator*(vnl_matrix<T> const& rhs) const
    {
        assert(this->cols() == rhs.rows());
        vnl_matrix<T> result(this->rows(), rhs.cols());
        const int l = this->rows();
        const int m = this->cols();
        const int n = rhs.cols();
        
        for(int i =0; i<l; ++i) {
            for(int k = 0; k<n; ++k) {
                T sum{0};
                for(int j =0; j<m; ++j) {
                    sum += (*this)(i, j) * rhs(j, k);
                }
                result(i, k) = sum;
            }
        }
        return result;
    }

    
    
    // fill and copy
    vnl_matrix& fill(const T v)
    {
        this->setConstant(v);
        return *this;
    }
    
    vnl_matrix& fill_diagonal(const T v)
    {
        for(int i = 0; i<this->rows() && i<this->cols(); ++i) {
            this->diagonal()[i] = v;
        }
        return *this;
    }
    
    
    //: Flatten row-major (C-style)
    vnl_vector<T> flatten_row_major() const
    {
        const unsigned int len = this->rows() * this->cols();
        vnl_vector<T> v(len);
        std::copy(this->data(), this->data()+len, v.data());
        return v;
    }
    
    //: Flatten column-major (Fortran-style)
    vnl_vector<T> flatten_column_major() const
    {
        vnl_vector<T> v(this->rows() * this->cols());
        for (unsigned int c = 0; c < this->cols(); ++c)
            for (unsigned int r = 0; r < this->rows(); ++r)
                v[c*this->rows()+r] = (*this)(r, c);
        return v;
    }
    
    /*
    //: Make a vector by applying a function across rows.
    vnl_vector<T> apply_rowwise(T (*f)(vnl_vector<T> const&)) const
    {
        
    }
    
    //: Make a vector by applying a function across columns.
    vnl_vector<T> apply_columnwise(T (*f)(vnl_vector<T> const&)) const
    {
        
    }
     */
    
    
    //: Return transpose
    vnl_matrix<T> transpose() const;
    
    //: Return conjugate transpose
    //vnl_matrix<T> conjugate_transpose() const;
    
    //: Return minimum value of elements
    T min_value() const {
        unsigned int r, c;
        T min_v = this->minCoeff(&r, &c);
        return min_v;
    }
    
    //: Return maximum value of elements
    T max_value() const {
        unsigned int r, c;
        T max_v = this->maxCoeff(&r, &c);
        return max_v;
    }
    
    //: Return location of minimum value of elements
    unsigned int arg_min() const {
        unsigned int r, c;
        this->minCoeff(&r, &c);
        return r * this->cols() + c;
    }
    
    //: Return location of maximum value of elements
    unsigned int arg_max() const {
        unsigned int r, c;
        this->maxCoeff(&r, &c);
        return r * this->cols() + c;
    }
    
    //: Set values of this matrix to those of M, starting at [top,left]
    vnl_matrix<T>& update(vnl_matrix<T> const& other, unsigned top=0, unsigned left=0)
    {
        assert(top + other.rows() <= this->rows());
        assert(left + other.cols() <= this->cols());
        
        for(int i = 0; i<other.rows(); ++i) {
            for(int j = 0; j<other.cols(); ++j) {
                (*this)(i+top, j+left) = other(i, j);
            }
        }
        return *this;
    }
    
    //: Extract a sub-matrix of size r x c, starting at (top,left)
    //  Thus it contains elements  [top,top+r-1][left,left+c-1]
    vnl_matrix<T> extract(unsigned r, unsigned c,
                          unsigned top=0, unsigned left=0) const
    {
        vnl_matrix<T> result = this->block(top, left, r, c);
        return result;
    }
    
    //: Make the matrix as if it had been default-constructed.
    void clear()
    {
        *this = vnl_matrix<T>();
    }
    
    
};


// Implementation

////--------------------------- Additions------------------------------------

//: Returns new matrix with rows and columns transposed.
// O(m*n).

template <class T>
vnl_matrix<T> vnl_matrix<T>::transpose() const
{
    vnl_matrix<T> result(this->cols(), this->rows());
    for (unsigned int i = 0; i < this->cols(); i++)
        for (unsigned int j = 0; j < this->rows(); j++)
            result(i, j) = (*this)(j, i);
    return result;
}

//:
// \relatesalso vnl_matrix
template<class T>
inline vnl_matrix<T> operator+(T const& value, vnl_matrix<T> const& m)
{
    return m+value;
}

template <class T>
vnl_matrix<T> operator- (T const& value, vnl_matrix<T> const& m)
{
    vnl_matrix<T> result(m.rows(),m.cols());
    for (unsigned int i = 0; i < m.rows(); i++)  // For each row
        for (unsigned int j = 0; j < m.cols(); j++) // For each element in column
            result(i, j) = T(value - m(i,j));    // subtract from value element.
    return result;
}


//template <class T> VNL_EXPORT m operator+(T const&, m const&);
//template <class T> VNL_EXPORT m operator-(T const&, m const&);
//template <class T> VNL_EXPORT m operator*(T const&, m const&);
template <class T>
vnl_matrix<T> element_product(vnl_matrix<T> const& a, vnl_matrix<T> const& b)
{
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    
    vnl_matrix<T> c(a.rows(), a.cols());
    for(int i = 0; i<a.rows(); ++i) {
        for(int j = 0; j<a.cols(); ++j) {
            c(i, j) = a(i, j) * b(i, j);
        }
    }
    return c;
}
template <class T>
vnl_matrix<T> element_quotient(vnl_matrix<T> const& a, vnl_matrix<T> const& b)
{
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    
    vnl_matrix<T> c(a.rows(), a.cols());
    for(int i = 0; i<a.rows(); ++i) {
        for(int j = 0; j<a.cols(); ++j) {
            c(i, j) = a(i, j) / b(i, j);
        }
    }
    return c;
}
template <class T>
T dot_product(vnl_matrix<T> const& a, vnl_matrix<T> const& b)
{
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    
    T sum{0};
    for(int i = 0; i<a.rows(); ++i) {
        for(int j = 0; j<a.cols(); ++j) {
            sum = a(i, j) * b(i, j);
        }
    }
    return sum;
    
}

//template <class T> VNL_EXPORT T inner_product(m const&, m const&);
//template <class T> VNL_EXPORT T cos_angle(m const&, m const& );
//template <class T> VNL_EXPORT std::ostream& operator<<(std::ostream&, m const&);
//template <class T> VNL_EXPORT std::istream& operator>>(std::istream&, m&);

#endif /* vnl_matrix_h */
