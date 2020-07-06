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
#include <vnl/vnl_numeric_traits.h>
#include <vnl/vnl_error.h>
#include <vnl/vnl_vector.h>

template <typename T> class vnl_matrix;
template <typename T> class vnl_vector;


template <typename T>
class vnl_matrix: public Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
{
    using base_class = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using abs_t = typename vnl_numeric_traits<T>::abs_t;
    //using abs_t = typename Eigen::NumTraits<T>::Real;
    //typedef typename vnl_numeric_traits<T>::abs_t abs_t;
public:
    vnl_matrix()=default;
    
    //: Construct a matrix of size r rows by c columns
    // Contents are unspecified.
    vnl_matrix(unsigned int r, unsigned int c):base_class(r, c){}
    
    
    //: Construct a matrix of size r rows by c columns, and all elements equal to v0
    // Complexity $O(r.c)$
    vnl_matrix(unsigned int r, unsigned int c, T const& v0):base_class(r, c){
        this->setConstant(v0);
    }
    
    //: Construct a matrix of size r rows by c columns, initialised by an automatic array
    // The first n elements, are initialised row-wise, to values.
    // Complexity $O(n)$
    vnl_matrix(unsigned int r, unsigned int c, unsigned int n, T const values[]):base_class(r, c)
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
    
    //: Copy construct a matrix
    // Complexity $O(r.c)$
    vnl_matrix(vnl_matrix<T> const& other) = default;    // from another matrix.
    
    /*
    // NOTE: move-assignment must be allowed to throw an exception, because we need to maintain
    //       backwards compatibility and the move-construction & move-aasignment
    //       operators fall back to the copy-assignment operator behavior in
    //       cases when the memory is externally managed.
    //: Move-constructor.
    vnl_matrix(vnl_matrix<T> &&);
    //: Move-assignment operator
    vnl_matrix<T>& operator=(vnl_matrix<T>&& rhs);
     */
    
    //: Matrix destructor
    virtual ~vnl_matrix() = default;
    
    // Basic 2D-Array functionality-------------------------------------------
    
    //: Return the total number of elements stored by the matrix.
    // This equals rows() * cols()
    inline unsigned int size() const { return base_class::size(); }
    
    //: Return the number of rows.
    //inline unsigned int rows() const { return this->rows(); }
    
    //: Return the number of columns.
    // A synonym for columns().
    //inline unsigned int cols() const { return this->cols; }
    
    //: Return the number of columns.
    // A synonym for cols().
    inline unsigned int columns() const { return base_class::cols(); }
    
    //: set element with boundary checks if error checking is on.
    inline void put(unsigned r, unsigned c, T const& v) {
        assert(r <= this->rows());
        assert( c <= this->cols());
        (*this)(r, c) = v;
    }
    
    //: get element with boundary checks if error checking is on.
    inline T get(unsigned r, unsigned c) const {
        assert(r <= this->rows());
        assert( c <= this->cols());
        return (*this)(r, c);
    }

   
    //: return pointer to given row
    // No boundary checking here.
    T       * operator[](unsigned r) { return this->row(r).data(); }
    
    //: return pointer to given row
    // No boundary checking here.
    T const * operator[](unsigned r) const { return this->row(r).data(); }
    
    
    //: Access an element for reading or writing
    // There are assert style boundary checks - #define NDEBUG to turn them off.
    //T       & operator()(unsigned r, unsigned c){return (*this)(i, j);}
    
    //: Access an element for reading
    // There are assert style boundary checks - #define NDEBUG to turn them off.
    //T const & operator()(unsigned r, unsigned c) const;
    
    // ----------------------- Filling and copying -----------------------
    
    //: Sets all elements of matrix to specified value, and returns "*this".
    //  Complexity $O(r.c)$
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to set a matrix to a column-normalized all-elements-equal matrix, say
    //  \code
    //     M.fill(1).normalize_columns();
    //  \endcode
    //  Returning "*this" also allows passing such a matrix as argument
    //  to a function f, without having to name the constructed matrix:
    //  \code
    //     f(vnl_matrix<double>(5,5,1.0).normalize_columns());
    //  \endcode
    vnl_matrix& fill(T const&);
    
    //: Sets all diagonal elements of matrix to specified value; returns "*this".
    //  Complexity $O(\min(r,c))$
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to set a 3x3 matrix to [5 0 0][0 10 0][0 0 15], just say
    //  \code
    //     M.fill_diagonal(5).scale_row(1,2).scale_column(2,3);
    //  \endcode
    //  Returning "*this" also allows passing a diagonal-filled matrix as argument
    //  to a function f, without having to name the constructed matrix:
    //  \code
    //     f(vnl_matrix<double>(3,3).fill_diagonal(5));
    //  \endcode
    vnl_matrix& fill_diagonal(T const&);
    
    //: Sets the diagonal elements of this matrix to the specified list of values.
    //  Returning "*this" allows "chaining" two or more operations: see the
    //  reasoning (and the examples) in the documentation for method
    //  fill_diagonal().
    vnl_matrix& set_diagonal(vnl_vector<T> const&);
    
    //: Fills (laminates) this matrix with the given data, then returns it.
    //  We assume that the argument points to a contiguous rows*cols array, stored rowwise.
    //  No bounds checking on the array.
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to fill a square matrix column-wise, fill it rowwise then transpose:
    //  \code
    //     M.copy_in(array).inplace_transpose();
    //  \endcode
    //  Returning "*this" also allows passing a filled-in matrix as argument
    //  to a function f, without having to name the constructed matrix:
    //  \code
    //     f(vnl_matrix<double>(3,3).copy_in(array));
    //  \endcode
    //vnl_matrix& copy_in(T const *);
    
    //: Fills (laminates) this matrix with the given data, then returns it.
    // A synonym for copy_in()
    //vnl_matrix& set(T const *d) { return copy_in(d); }
    
    //: Fills the given array with this matrix.
    //  We assume that the argument points to a contiguous rows*cols array, stored rowwise.
    // No bounds checking on the array.
    //void copy_out(T *) const;
    
    //: Set all elements to value v
    // Complexity $O(r.c)$
    vnl_matrix<T>& operator=(T const&v)
    {
        this->setConstant(v);
        return *this;
    }
    
    //: Copies all elements of rhs matrix into lhs matrix.
    // Complexity $O(\min(r,c))$
    // vnl_matrix<T>& operator=(vnl_matrix<T> const&);
    
    // ----------------------- Arithmetic --------------------------------
    // note that these functions should not pass scalar as a const&.
    // Look what would happen to A /= A(0,0).
    
    //: Add rhs to each element of lhs matrix in situ
    vnl_matrix<T>& operator+=(T value);
    
    //: Subtract rhs from each element of lhs matrix in situ
    vnl_matrix<T>& operator-=(T value);
    
    //: Scalar multiplication in situ of lhs matrix  by rhs
    vnl_matrix<T>& operator*=(T value);
   
    //: Scalar division of lhs matrix  in situ by rhs
    vnl_matrix<T>& operator/=(T value);
    
    //: Add rhs to lhs  matrix in situ
    vnl_matrix<T>& operator+=(vnl_matrix<T> const& rhs);
    
    //: Subtract rhs from lhs matrix in situ
    vnl_matrix<T>& operator-=(vnl_matrix<T> const& rhs);
    
    //: Multiply lhs matrix in situ by rhs
    vnl_matrix<T>& operator*=(vnl_matrix<T> const&rhs);
    
    //: Negate all elements of matrix
    vnl_matrix<T> operator-() const;
    
    //: Add rhs to each element of lhs matrix and return result in new matrix
    vnl_matrix<T> operator+(T const& v) const;
    
    //: Subtract rhs from each element of lhs matrix and return result in new matrix
    vnl_matrix<T> operator-(T const& v) const;
    
    //: Scalar multiplication of lhs matrix by rhs  and return result in new matrix
    vnl_matrix<T> operator*(T const& v) const;
    
    //: Scalar division of lhs matrix by rhs and return result in new matrix
    vnl_matrix<T> operator/(T const& v) const;
    
    //: Matrix add rhs to lhs matrix and return result in new matrix
    vnl_matrix<T> operator+(vnl_matrix<T> const& rhs) const;
    
    //: Matrix subtract rhs from lhs and return result in new matrix
    vnl_matrix<T> operator-(vnl_matrix<T> const& rhs) const;
    
    //: Matrix multiply lhs by rhs matrix and return result in new matrix
    vnl_matrix<T> operator*(vnl_matrix<T> const& rhs) const;
    
    ////--------------------------- Additions ----------------------------
    
    //: Make a new matrix by applying function to each element.
    vnl_matrix<T> apply(T (*f)(T)) const;
    
    //: Make a new matrix by applying function to each element.
    vnl_matrix<T> apply(T (*f)(T const&)) const;
    
    //: Make a vector by applying a function across rows.
    vnl_vector<T> apply_rowwise(T (*f)(vnl_vector<T> const&)) const;
    
    //: Make a vector by applying a function across columns.
    vnl_vector<T> apply_columnwise(T (*f)(vnl_vector<T> const&)) const;
    
    //: Return transpose
    vnl_matrix<T> transpose() const;
    
    //: Return conjugate transpose
    vnl_matrix<T> conjugate_transpose() const;
    
    //: Set values of this matrix to those of M, starting at [top,left]
    vnl_matrix<T>& update(vnl_matrix<T> const&, unsigned top=0, unsigned left=0);
    
    //: Set the elements of the i'th column to v[i]  (No bounds checking)
    vnl_matrix& set_column(unsigned i, T const * v);
    
    //: Set the elements of the i'th column to value, then return *this.
    vnl_matrix& set_column(unsigned i, T value );
    
    //: Set j-th column to v, then return *this.
    vnl_matrix& set_column(unsigned j, vnl_vector<T> const& v);
    
    //: Set columns to those in M, starting at starting_column, then return *this.
    vnl_matrix& set_columns(unsigned starting_column, vnl_matrix<T> const& M);
    
    //: Set the elements of the i'th row to v[i]  (No bounds checking)
    vnl_matrix& set_row(unsigned i, T const * v);
    
    //: Set the elements of the i'th row to value, then return *this.
    vnl_matrix& set_row(unsigned i, T value );
    
    //: Set the i-th row
    vnl_matrix& set_row(unsigned i, vnl_vector<T> const&);
    
    //: Extract a sub-matrix of size r x c, starting at (top,left)
    //  Thus it contains elements  [top,top+r-1][left,left+c-1]
    vnl_matrix<T> extract(unsigned r, unsigned c,
                          unsigned top=0, unsigned left=0) const;
    
    //: Extract a sub-matrix starting at (top,left)
    //
    //  The output is stored in \a sub_matrix, and it should have the
    //  required size on entry.  Thus the result will contain elements
    //  [top,top+sub_matrix.rows()-1][left,left+sub_matrix.cols()-1]
    void extract ( vnl_matrix<T>& sub_matrix,
                  unsigned top=0, unsigned left=0) const;
    
    //: Get a vector equal to the given row
    vnl_vector<T> get_row(unsigned r) const;
    
    //: Get a vector equal to the given column
    vnl_vector<T> get_column(unsigned c) const;
    
    //: Get a matrix composed of rows from the indices specified in the supplied vector.
    vnl_matrix<T> get_rows(const vnl_vector<unsigned int>& i) const;
    
    //: Get a matrix composed of columns from the indices specified in the supplied vector.
    vnl_matrix<T> get_columns(const vnl_vector<unsigned int>& i) const;
    
    //: Get n rows beginning at rowstart
    vnl_matrix<T> get_n_rows(unsigned rowstart, unsigned n) const;
    
    //: Get n columns beginning at colstart
    vnl_matrix<T> get_n_columns(unsigned colstart, unsigned n) const;
    
    //: Return a vector with the content of the (main) diagonal
    vnl_vector<T> get_diagonal() const;
    
    //: Flatten row-major (C-style)
    vnl_vector<T> flatten_row_major() const;
    
    //: Flatten column-major (Fortran-style)
    vnl_vector<T> flatten_column_major() const;
    
    // ==== mutators ====
    
    //: Sets this matrix to an identity matrix, then returns "*this".
    //  Returning "*this" allows e.g. passing an identity matrix as argument to
    //  a function f, without having to name the constructed matrix:
    //  \code
    //     f(vnl_matrix<double>(5,5).set_identity());
    //  \endcode
    //  Returning "*this" also allows "chaining" two or more operations:
    //  e.g., to set a 3x3 matrix to [3 0 0][0 2 0][0 0 1], one could say
    //  \code
    //     M.set_identity().scale_row(0,3).scale_column(1,2);
    //  \endcode
    //  If the matrix is not square, anyhow set main diagonal to 1, the rest to 0.
    vnl_matrix& set_identity();
    
    //: Transposes this matrix efficiently, and returns it.
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to fill a square matrix column-wise, fill it rowwise then transpose:
    //  \code
    //     M.copy_in(array).inplace_transpose();
    //  \endcode
    vnl_matrix& inplace_transpose();
    
    //: Reverses the order of rows, and returns "*this".
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to flip both up-down and left-right, one could just say
    //  \code
    //     M.flipud().fliplr();
    //  \endcode
    vnl_matrix& flipud();
    
    //: Reverses the order of columns, and returns "*this".
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to flip both up-down and left-right, one could just say
    //  \code
    //     M.flipud().fliplr();
    //  \endcode
    vnl_matrix& fliplr();
    
    //: Normalizes each row so it is a unit vector, and returns "*this".
    //  Zero rows are not modified
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to set a matrix to a row-normalized all-elements-equal matrix, say
    //  \code
    //     M.fill(1).normalize_rows();
    //  \endcode
    //  Returning "*this" also allows passing such a matrix as argument
    //  to a function f, without having to name the constructed matrix:
    //  \code
    //     f(vnl_matrix<double>(5,5,1.0).normalize_rows());
    //  \endcode
    vnl_matrix& normalize_rows();
    
    //: Normalizes each column so it is a unit vector, and returns "*this".
    //  Zero columns are not modified
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to set a matrix to a column-normalized all-elements-equal matrix, say
    //  \code
    //     M.fill(1).normalize_columns();
    //  \endcode
    //  Returning "*this" also allows passing such a matrix as argument
    //  to a function f, without having to name the constructed matrix:
    //  \code
    //     f(vnl_matrix<double>(5,5,1.0).normalize_columns());
    //  \endcode
    vnl_matrix& normalize_columns();
    
    //: Scales elements in given row by a factor T, and returns "*this".
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to set a 3x3 matrix to [3 0 0][0 2 0][0 0 1], one could say
    //  \code
    //     M.set_identity().scale_row(0,3).scale_column(1,2);
    //  \endcode
    vnl_matrix& scale_row(unsigned row, T value);
    
    //: Scales elements in given column by a factor T, and returns "*this".
    //  Returning "*this" allows "chaining" two or more operations:
    //  e.g., to set a 3x3 matrix to [3 0 0][0 2 0][0 0 1], one could say
    //  \code
    //     M.set_identity().scale_row(0,3).scale_column(1,2);
    //  \endcode
    vnl_matrix& scale_column(unsigned col, T value);
    
    //: Swap this matrix with that matrix
    void swap(vnl_matrix<T> & that) noexcept;
    
    
    //: Type def for norms.
    //typedef typename vnl_c_vector<T>::abs_t abs_t;
    /*
    //: Return sum of absolute values of elements
    abs_t array_one_norm() const { return vnl_c_vector<T>::one_norm(begin(), size()); }
    
    //: Return square root of sum of squared absolute element values
    abs_t array_two_norm() const { return vnl_c_vector<T>::two_norm(begin(), size()); }
    
    //: Return largest absolute element value
    abs_t array_inf_norm() const { return vnl_c_vector<T>::inf_norm(begin(), size()); }
    
    //: Return sum of absolute values of elements
    abs_t absolute_value_sum() const { return array_one_norm(); }
    
    //: Return largest absolute value
    abs_t absolute_value_max() const { return array_inf_norm(); }
    
    // $ || M ||_1 := \max_j \sum_i | M_{ij} | $
    abs_t operator_one_norm() const;
    
    // $ || M ||_\inf := \max_i \sum_j | M_{ij} | $
    abs_t operator_inf_norm() const;
    */
    //: Return Frobenius norm of matrix (sqrt of sum of squares of its elements)
    abs_t frobenius_norm() const { return this->norm(); }
    
    //: Return Frobenius norm of matrix (sqrt of sum of squares of its elements)
    abs_t fro_norm() const { return frobenius_norm(); }
    /*
    //: Return RMS of all elements
    abs_t rms() const { return vnl_c_vector<T>::rms_norm(begin(), size()); }
    */
    //: Return minimum value of elements
    T min_value() const;
    
    //: Return maximum value of elements
    T max_value() const;
    
    //: Return location of minimum value of elements
    unsigned arg_min() const;
    //{ return vnl_c_vector<T>::arg_min(begin(), size()); }
    
    //: Return location of maximum value of elements
    unsigned arg_max() const;
    //{ return vnl_c_vector<T>::arg_max(begin(), size()); }
    
    //: Return mean of all matrix elements
    T mean() const;
    //{ return vnl_c_vector<T>::mean(begin(), size()); }
    
    // predicates
    
    //: Return true iff the size is zero.
    bool empty() const;
    //{ return !data || !num_rows || !num_cols; }
    
    //:  Return true if all elements equal to identity.
    bool is_identity() const;
    
    //:  Return true if all elements equal to identity, within given tolerance
    bool is_identity(double tol) const;
    
    //: Return true if all elements equal to zero.
    bool is_zero() const;
    
    //: Return true if all elements equal to zero, within given tolerance
    bool is_zero(double tol) const;
    
    //:  Return true if all elements of both matrices are equal, within given tolerance
    bool is_equal(vnl_matrix<T> const& rhs, double tol) const;
    
    //: Return true if finite
    bool is_finite() const;
    
    //: Return true if matrix contains NaNs
    bool has_nans() const;
    
    ////----------------------- Input/Output ----------------------------
    
    //: Read a vnl_matrix from an ascii std::istream, automatically determining file size if the input matrix has zero size.
    static vnl_matrix<T> read(std::istream& s);
    
    // : Read a vnl_matrix from an ascii std::istream, automatically determining file size if the input matrix has zero size.
    bool read_ascii(std::istream& s);
    
    //--------------------------------------------------------------------------------
    
    //: Access the contiguous block storing the elements in the matrix row-wise. O(1).
    // 1d array, row-major order.
    T const* data_block() const;// { return data[0]; }
    
    //: Access the contiguous block storing the elements in the matrix row-wise. O(1).
    // 1d array, row-major order.
    T      * data_block();// { return data[0]; }
    
    //: Access the 2D array, so that elements can be accessed with array[row][col] directly.
    //  2d array, [row][column].
    T const* const* data_array() const;// { return data; }
    
    //: Access the 2D array, so that elements can be accessed with array[row][col] directly.
    //  2d array, [row][column].
    T      *      * data_array();// { return data; }
    
    
    //typedef T element_type;
    
    
    //: Iterators
    typedef T       *iterator;
    //: Iterator pointing to start of data
    iterator       begin() { return this->data(); }
    //: Iterator pointing to element beyond end of data
    iterator       end() { return this->data()+this->size(); }
    
    //: Const iterators
    typedef T const *const_iterator;
    //: Iterator pointing to start of data
    const_iterator begin() const { return this->data(); }
    //: Iterator pointing to element beyond end of data
    const_iterator end() const { return this->data()+this->size(); }
    
     /*
    //: Return a reference to this.
    // Useful in code which would prefer not to know if its argument
    // is a matrix, matrix_ref or a matrix_fixed.  Note that it doesn't
    // return a matrix_ref, so it's only useful in templates or macros.
    vnl_matrix<T> const& as_ref() const { return *this; }
    
    //: Return a reference to this.
    vnl_matrix<T>&       as_ref()       { return *this; }
    */
    //--------------------------------------------------------------------------------
    
    //: Return true if *this == rhs
    bool operator_eq(vnl_matrix<T> const & rhs) const;
    
    //: Equality operator
    //bool operator==(vnl_matrix<T> const &that) const; //{ return  this->operator_eq(that); }
    
    //: Inequality operator
    //bool operator!=(vnl_matrix<T> const &that) const;// { return !this->operator_eq(that); }
    
    //: Print matrix to os in some hopefully sensible format
    void print(std::ostream& os) const;
    
    //: Make the matrix as if it had been default-constructed.
    void clear();
    
    //: Resize to r rows by c columns. Old data lost.
    // Returns true if size changed.
    bool set_size(unsigned r, unsigned c);
};


// Implementation
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::fill(const T& v)
{
    this->setConstant(v);
    return *this;
}

template <typename T>
vnl_matrix<T>& vnl_matrix<T>::fill_diagonal(const T& v)
{
    for(int i = 0; i<this->rows() && i<this->cols(); ++i) {
        this->diagonal()[i] = v;
    }
    return *this;
}

//: Sets the diagonal elements of this matrix to the specified list of values.

template <class T>
vnl_matrix<T>& vnl_matrix<T>::set_diagonal(vnl_vector<T> const& diag)
{
    assert(diag.size() >= this->rows() ||
           diag.size() >= this->cols());
    // The length of the diagonal of a non-square matrix is the minimum of
    // the matrix's width & height; that explains the "||" in the assert,
    // and the "&&" in the upper bound for the "for".
    for (unsigned int i = 0; i < this->rows() && i < this->cols(); ++i) {
        (*this)(i, i) = diag[i];
    }
    return *this;
}

//: Add rhs to each element of lhs matrix in situ
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::operator+=(T value)
{
    this->array() += value;
    return *this;
}

//: Subtract rhs from each element of lhs matrix in situ
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::operator-=(T value)
{
    this->array() -= value;
    return *this;
}

//: Scalar multiplication in situ of lhs matrix  by rhs
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::operator*=(T value)
{
    this->array() *= value;
    return *this;
}

//: Scalar division of lhs matrix  in situ by rhs
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::operator/=(T value)
{
    this->array() /= value;
    return *this;
}

//: Add rhs to lhs  matrix in situ
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::operator+=(vnl_matrix<T> const& rhs)
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
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::operator-=(vnl_matrix<T> const& rhs)
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
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::operator*=(vnl_matrix<T> const&rhs)
{
    return *this = (*this) * rhs;
}

//: Negate all elements of matrix
template <typename T>
vnl_matrix<T> vnl_matrix<T>::operator-() const
{
    vnl_matrix<T> result(this->rows(), this->cols());
    for (unsigned int i = 0; i < this->rows(); i++)
        for (unsigned int j = 0; j < this->cols(); j++)
            result(i, j) = - (*this)(i, j);
    return result;
}

//: Add rhs to each element of lhs matrix and return result in new matrix
template <typename T>
vnl_matrix<T> vnl_matrix<T>::operator+(T const& v) const {
    vnl_matrix<T> result(this->rows(), this->cols());
    const unsigned int n = this->rows() * this->cols();
    T const *m = this->data();
    T *dst = result.data();
    
    for (unsigned int i = 0; i < n; ++i)
        dst[i] = T(m[i] + v);
    return result;
}

//: Subtract rhs from each element of lhs matrix and return result in new matrix
template <typename T>
vnl_matrix<T> vnl_matrix<T>::operator-(T const& v) const {
    vnl_matrix<T> result(this->rows(), this->cols());
    const unsigned int n = this->rows() * this->cols();
    T const *m = this->data();
    T *dst = result.data();
    
    for (unsigned int i = 0; i < n; ++i)
        dst[i] = T(m[i] - v);
    return result;
}

//: Scalar multiplication of lhs matrix by rhs  and return result in new matrix
template <typename T>
vnl_matrix<T> vnl_matrix<T>::operator*(T const& v) const {
    vnl_matrix<T> result(this->rows(), this->cols());
    const unsigned int n = this->rows() * this->cols();
    T const *m = this->data();
    T *dst = result.data();
    
    for (unsigned int i = 0; i < n; ++i)
        dst[i] = T(m[i] * v);
    return result;
}

//: Scalar division of lhs matrix by rhs and return result in new matrix
template <typename T>
vnl_matrix<T> vnl_matrix<T>::operator/(T const& v) const {
    vnl_matrix<T> result(this->rows(), this->cols());
    const unsigned int n = this->rows() * this->cols();
    T const *m = this->data();
    T *dst = result.data();
    
    for (unsigned int i = 0; i < n; ++i)
        dst[i] = T(m[i] / v);
    return result;
}

//: Matrix add rhs to lhs matrix and return result in new matrix
template <typename T>
vnl_matrix<T> vnl_matrix<T>::operator+(vnl_matrix<T> const& rhs) const
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
template <typename T>
vnl_matrix<T> vnl_matrix<T>::operator-(vnl_matrix<T> const& rhs) const
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
template <typename T>
vnl_matrix<T> vnl_matrix<T>::operator*(vnl_matrix<T> const& rhs) const
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

//: Return the matrix made by applying "f" to each element.
template <typename T>
vnl_matrix<T> vnl_matrix<T>::apply(T (*f)(T const&)) const
{
    vnl_matrix<T> ret(this->rows(), this->cols());
    for(int i = 0; i<this->rows(); ++i) {
        for(int j = 0; j<this->cols(); ++j) {
            ret(i, j) = f(this(i, j));
        }
    }
    //vnl_c_vector<T>::apply(this->data[0], num_rows * num_cols, f, ret.data_block());
    return ret;
}

//: Return the matrix made by applying "f" to each element.
template <typename T>
vnl_matrix<T> vnl_matrix<T>::apply(T (*f)(T)) const
{
    vnl_matrix<T> ret(this->rows(), this->cols());
    for(int i = 0; i<this->rows(); ++i) {
        for(int j = 0; j<this->cols(); ++j) {
            ret(i, j) = f((*this)(i, j));
        }
    }
    return ret;
}

//: Make a vector by applying a function across rows.
template <typename T>
vnl_vector<T> vnl_matrix<T>::apply_rowwise(T (*f)(vnl_vector<T> const&)) const
{
    vnl_vector<T> v(this->rows());
    for (unsigned int i = 0; i < this->rows(); ++i)
        v.put(i,f(this->get_row(i)));
    return v;
}

//: Make a vector by applying a function across columns.
template <typename T>
vnl_vector<T> vnl_matrix<T>::apply_columnwise(T (*f)(vnl_vector<T> const&)) const
{
    vnl_vector<T> v(this->cols());
    for (unsigned int i = 0; i < this->cols(); ++i)
        v.put(i,f(this->get_column(i)));
    return v;
}

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

//: Return conjugate transpose
template <typename T>
vnl_matrix<T> vnl_matrix<T>::conjugate_transpose() const
{
    // matrix.adjoint() (or matrix.conjugate().transpose()
    //https://stackoverflow.com/questions/30451040/complex-number-matrix-multiplication-eigen-vs-matlab
    return this->conjugate().transpose();
}

//: Set values of this matrix to those of M, starting at [top,left]
template <class T>
vnl_matrix<T>& vnl_matrix<T>::update(vnl_matrix<T> const& other, unsigned top, unsigned left)
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
//--------------------------------------------------------------------------------

//: Set row[row_index] to data at given address. No bounds check.
template <class T>
vnl_matrix<T>& vnl_matrix<T>::set_row(unsigned row_index, T const *v)
{
    for (unsigned int j = 0; j < this->cols(); j++)    // For each element in row
        (*this)(row_index, j) = v[j];
    return *this;
}

//: Set row[row_index] to given vector.
template <class T>
vnl_matrix<T>& vnl_matrix<T>::set_row(unsigned row_index, vnl_vector<T> const &v)
{
#ifndef NDEBUG
    if (v.size() != this->cols())
        vnl_error_vector_dimension ("vnl_matrix::set_row", v.size(), this->cols());
#endif
    set_row(row_index,v.data_block());
    return *this;
}

//: Set row[row_index] to given value.
template <class T>
vnl_matrix<T>& vnl_matrix<T>::set_row(unsigned row_index, T v)
{
    for (unsigned int j = 0; j < this->num_cols; j++)    // For each element in row
        this->data[row_index][j] = v;
    return *this;
}

//--------------------------------------------------------------------------------

//: Set column[column_index] to data at given address.
template <class T>
vnl_matrix<T>& vnl_matrix<T>::set_column(unsigned column_index, T const *v)
{
    for (unsigned int i = 0; i < this->rows(); i++)    // For each element in row
        (*this)(i, column_index) = v[i];
    return *this;
}

//: Set column[column_index] to given vector.
template <class T>
vnl_matrix<T>& vnl_matrix<T>::set_column(unsigned column_index, vnl_vector<T> const &v)
{
#ifndef NDEBUG
    if (v.size() != this->rows())
        vnl_error_vector_dimension ("vnl_matrix::set_column", v.size(), this->rows());
#endif
    set_column(column_index,v.data_block());
    return *this;
}

//: Set column[column_index] to given value.
template <class T>
vnl_matrix<T>& vnl_matrix<T>::set_column(unsigned column_index, T v)
{
    for (unsigned int j = 0; j < this->num_rows; j++)    // For each element in row
        (*this)(j, column_index) = v;
    return *this;
}


//: Set columns starting at starting_column to given matrix
template <class T>
vnl_matrix<T>& vnl_matrix<T>::set_columns(unsigned starting_column, vnl_matrix<T> const& m)
{
#ifndef NDEBUG
    if (this->num_rows != m.num_rows ||
        this->num_cols < m.num_cols + starting_column)           // Size match?
        vnl_error_matrix_dimension ("set_columns",
                                    this->num_rows, this->num_cols,
                                    m.num_rows, m.num_cols);
#endif
    
    for (unsigned int j = 0; j < m.num_cols; ++j)
        for (unsigned int i = 0; i < this->num_rows; i++)    // For each element in row
            (*this)(i, starting_column + j) = m(i, j);
    return *this;
}

//: Extract a sub-matrix of size r x c, starting at (top,left)
//  Thus it contains elements  [top,top+r-1][left,left+c-1]
template<typename T>
vnl_matrix<T> vnl_matrix<T>::extract(unsigned r, unsigned c,
                      unsigned top, unsigned left) const
{
    vnl_matrix<T> result = this->block(top, left, r, c);
    return result;
}

template <class T>
void vnl_matrix<T>::extract( vnl_matrix<T>& submatrix,
                            unsigned top, unsigned left) const {
    unsigned const rowz = submatrix.rows();
    unsigned const colz = submatrix.cols();
#ifndef NDEBUG
    unsigned int bottom = top + rowz;
    unsigned int right = left + colz;
    if ((this->rows() < bottom) || (this->cols() < right))
        vnl_error_matrix_dimension ("extract",
                                    this->rows(), this->cols(), bottom, right);
#endif
    for (unsigned int i = 0; i < rowz; i++)      // actual copy of all elements
        for (unsigned int j = 0; j < colz; j++)    // in submatrix
            submatrix(i, j) = (*this)(top+i, left+j);
}

//: Returns a copy of n rows, starting from "row"
template <class T>
vnl_matrix<T> vnl_matrix<T>::get_n_rows (unsigned row, unsigned n) const
{
#ifndef NDEBUG
    if (row + n > this->num_rows)
        vnl_error_matrix_row_index ("get_n_rows", row);
#endif
    unsigned int i = row;
    unsigned int j = 0;
    unsigned int p = n;
    unsigned int q = this->cols();
    
    vnl_matrix<T> result = this->block(i, j, p, q);
    return result;
    
    // Extract data rowwise.
    //return vnl_matrix<T>(data[row], n, this->num_cols);
}

//: Returns a copy of n columns, starting from "column".
template <class T>
vnl_matrix<T> vnl_matrix<T>::get_n_columns (unsigned column, unsigned n) const
{
#ifndef NDEBUG
    if (column + n > this->num_cols)
        vnl_error_matrix_col_index ("get_n_columns", column);
#endif
    unsigned int i = 0;
    unsigned int j = column;
    unsigned int p = this->rows();
    unsigned int q = n;
    
    vnl_matrix<T> result = this->block(i, j, p, q);
    //for (unsigned int c = 0; c < n; ++c)
    //    for (unsigned int r = 0; r < this->num_rows; ++r)
    //        result(r, c) = data[r][column + c];
    return result;
}

//: Create a vector out of row[row_index].
template <class T>
vnl_vector<T> vnl_matrix<T>::get_row(unsigned row_index) const
{
#ifdef ERROR_CHECKING
    if (row_index >= this->num_rows)
        vnl_error_matrix_row_index ("get_row", row_index);
#endif
    
    vnl_vector<T> v = this->base_class::row(row_index);
    //for (unsigned int j = 0; j < this->num_cols; j++)    // For each element in row
    //    v[j] = this->data[row_index][j];
    return v;
}

//: Create a vector out of column[column_index].
template <class T>
vnl_vector<T> vnl_matrix<T>::get_column(unsigned column_index) const
{
#ifdef ERROR_CHECKING
    if (column_index >= this->num_cols)
        vnl_error_matrix_col_index ("get_column", column_index);
#endif
    
    vnl_vector<T> v = this->base_class::col(column_index);//(this->num_rows);
    //for (unsigned int j = 0; j < this->num_rows; j++)    // For each element in row
    //    v[j] = this->data[j][column_index];
    return v;
}

//: Create a vector out of row[row_index].
template <class T>
vnl_matrix<T> vnl_matrix<T>::get_rows(const vnl_vector<unsigned int>& i) const
{
    vnl_matrix<T> m(i.size(), this->cols());
    for (unsigned int j = 0; j < i.size(); ++j)
        m.set_row(j, this->get_row(i.get(j)));
    return m;
}

//: Create a vector out of column[column_index].
template <class T>
vnl_matrix<T> vnl_matrix<T>::get_columns(const vnl_vector<unsigned int>& i) const
{
    vnl_matrix<T> m(this->rows(), i.size());
    for (unsigned int j = 0; j < i.size(); ++j)
        m.set_column(j, this->get_column(i.get(j)));
    return m;
}

template<typename T>
vnl_vector<T> vnl_matrix<T>::get_diagonal() const
{
    size_t n = std::min(this->rows(), this->cols());
    vnl_vector<T> v(n);
    for (unsigned int j = 0; j < n; ++j)
        v[j] = (*this)(j, j);
    return v;
}

//: Flatten row-major (C-style)
template<typename T>
vnl_vector<T> vnl_matrix<T>::flatten_row_major() const
{
    const unsigned int len = this->rows() * this->cols();
    vnl_vector<T> v(len);
    std::copy(this->data(), this->data()+len, v.data());
    return v;
}

//: Flatten column-major (Fortran-style)
template<typename T>
vnl_vector<T> vnl_matrix<T>::flatten_column_major() const
{
    vnl_vector<T> v(this->rows() * this->cols());
    for (unsigned int c = 0; c < this->cols(); ++c)
        for (unsigned int r = 0; r < this->rows(); ++r)
            v[c*this->rows()+r] = (*this)(r, c);
    return v;
}

//: Make each row of the matrix have unit norm.
// All-zero rows are ignored.
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::normalize_rows()
{
    typedef typename vnl_numeric_traits<T>::abs_t Abs_t;
    typedef typename vnl_numeric_traits<T>::real_t Real_t;
    typedef typename vnl_numeric_traits<Real_t>::abs_t abs_real_t;
    for (unsigned int i = 0; i < this->rows(); ++i) {  // For each row in the Matrix
        Abs_t norm(0); // double will not do for all types.
        for (unsigned int j = 0; j < this->cols(); ++j)  // For each element in row
            norm += vnl_math::squared_magnitude((*this)(i, j));
        
        if (norm != 0) {
            abs_real_t scale = abs_real_t(1)/(std::sqrt((abs_real_t)norm));
            for (unsigned int j = 0; j < this->cols(); ++j)
                (*this)(i, j) = T(Real_t((*this)(i, j)) * scale);
        }
    }
    return *this;
}

//: Make each column of the matrix have unit norm.
// All-zero columns are ignored.
template <typename T>
vnl_matrix<T>& vnl_matrix<T>::normalize_columns()
{
    typedef typename vnl_numeric_traits<T>::abs_t Abs_t;
    typedef typename vnl_numeric_traits<T>::real_t Real_t;
    typedef typename vnl_numeric_traits<Real_t>::abs_t abs_real_t;
    for (unsigned int j = 0; j < this->cols(); j++) {  // For each column in the Matrix
        Abs_t norm(0); // double will not do for all types.
        for (unsigned int i = 0; i < this->rows(); i++)
            norm += vnl_math::squared_magnitude((*this)(i, j));
        
        if (norm != 0) {
            abs_real_t scale = abs_real_t(1)/(std::sqrt((abs_real_t)norm));
            for (unsigned int i = 0; i < this->rows(); i++)
                (*this)(i, j) = T(Real_t((*this)(i, j)) * scale);
        }
    }
    return *this;
}


//: Return location of minimum value of elements
template<typename T>
unsigned int vnl_matrix<T>::arg_min() const {
    unsigned int r, c;
    this->minCoeff(&r, &c);
    return r * this->cols() + c;
}

//: Return location of maximum value of elements
template<typename T>
unsigned int vnl_matrix<T>::arg_max() const {
    unsigned int r, c;
    this->maxCoeff(&r, &c);
    return r * this->cols() + c;
}

//: Return minimum value of elements
template<typename T>
T vnl_matrix<T>::min_value() const {
    unsigned int r, c;
    T min_v = this->minCoeff(&r, &c);
    return min_v;
}

//: Return maximum value of elements
template<typename T>
T vnl_matrix<T>::max_value() const {
    unsigned int r, c;
    T max_v = this->maxCoeff(&r, &c);
    return max_v;
}

//: Make the matrix as if it had been default-constructed.
template <class T>
void vnl_matrix<T>::clear()
{
    *this = vnl_matrix<T>();
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
