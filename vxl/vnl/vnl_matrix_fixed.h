//
//  vnl_matrix_fixed.h
//  vxl
//
//  Created by Jimmy Chen on 2020-06-14.
//  Copyright © 2020 Nowhere Planet. All rights reserved.
//

#ifndef vnl_matrix_fixed_h
#define vnl_matrix_fixed_h

#include <Eigen/Dense>

template <typename T, unsigned int num_rows, unsigned int num_cols> class vnl_matrix_fixed;
template <typename T> class vnl_vector;

template <typename T, unsigned int num_rows, unsigned int num_cols>
class vnl_matrix_fixed : public Eigen::Matrix<T, num_rows, num_cols, Eigen::RowMajor>
{
    static constexpr size_t num_elements = num_rows*num_cols;
    static constexpr size_t num_bytes = num_elements*sizeof(T);
    using base_class = Eigen::Matrix<T, num_rows, num_cols, Eigen::RowMajor>;
public:
    //: Construct an empty num_rows*num_cols matrix
    vnl_matrix_fixed() = default;
    
    //: Construct an m*n Matrix and copy rhs into it.
    //  Abort if rhs is not the same size.
    //vnl_matrix_fixed(const vnl_matrix_fixed& rhs) = default;
    
    // This constructor allows us to construct vnl_matrix_fixed from Eigen expressions
    // https://eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html
    template<typename OtherDerived>
    vnl_matrix_fixed(const Eigen::MatrixBase<OtherDerived>& other):
        base_class(other)
    {}
    
    // This method allows us to assign Eigen expressions to vnl_matrix_fixed
    template<typename OtherDerived>
    vnl_matrix_fixed& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        this->base_class::operator=(other);
        return *this;
    }
    
    
    
    //: Construct an m*n matrix and fill with value
    explicit vnl_matrix_fixed(T value)
    {
        this->setConstant(value);
    }
    
    //: Construct an m*n Matrix and copy data into it row-wise.
    explicit vnl_matrix_fixed(const T* datablck)
    {
        std::memcpy(this->data(), datablck, num_bytes);
    }
    
    //: Return a vector with the content of the (main) diagonal
    vnl_vector<T> get_diagonal() const
    {
        return this->diagonal();
    }
    
    
    // fill and copy
    vnl_matrix_fixed& fill(const T v)
    {
        this->setConstant(v);
        return *this;
    }
    
    // ----------------------- Arithmetic --------------------------------
    // note that these functions should not pass scalar as a const&.
    // Look what would happen to A /= A(0,0).
    
    //: Add \a s to each element of lhs matrix in situ
    vnl_matrix_fixed& operator+= (T s)
    {
        this->array() += s;
        return *this;
    }
    
    //: Subtract \a s from each element of lhs matrix in situ
    vnl_matrix_fixed& operator-= (T s)
    {
        this->array() -= s;
        return *this;
    }
    
    //:
    vnl_matrix_fixed& operator*= (T s)
    {
        this->array() *= s;
        return *this;
    }
    
    //:
    vnl_matrix_fixed& operator/= (T s)
    {
        this->array() /= s;
        return *this;
    }
    
    //:
    vnl_matrix_fixed& operator+= (vnl_matrix_fixed const& m)
    {
        add( this->data(), m.data(), this->data() ); return *this;
    }
    
    /*
    //:
    vnl_matrix_fixed& operator+= (vnl_matrix<T> const& m)
    {
        assert( m.rows() == rows() && m.cols() == cols() );
        self::add( data_block(), m.data_block(), data_block() ); return *this;
    }
     */
    
    //:
    vnl_matrix_fixed& operator-= (vnl_matrix_fixed const& m)
    {
        sub( this->data(), m.data(), this->data()); return *this;
    }
    
    /*
    //:
    vnl_matrix_fixed& operator-= (vnl_matrix<T> const& m)
    {
        assert( m.rows() == rows() && m.cols() == cols() );
        self::sub( data_block(), m.data_block(), data_block() );
        return *this;
    }
     */
    
    //: Negate all elements of matrix
    vnl_matrix_fixed operator- () const
    {
        vnl_matrix_fixed r;
        sub( T(0), this->data(), r.data() );
        return r;
    }
    
    //:
    vnl_matrix_fixed& operator*= (vnl_matrix_fixed<T,num_cols,num_cols> const& s)
    {
        vnl_matrix_fixed<T, num_rows, num_cols> out;
        for (unsigned i = 0; i < num_rows; ++i){
            for (unsigned j = 0; j < num_cols; ++j)
            {
                T accum = this->data()[i*num_cols] * s(0,j);
                for (unsigned k = 1; k < num_cols; ++k)
                    accum += this->data()[i*num_cols+k] * s(k,j);
                out(i,j) = accum;
            }
        }
        return *this = out;
    }

    
    
    // Helper routines for arithmetic. These routines know the size from
    // the template parameters. The vector-vector operations are
    // element-wise.
    
    static void add( const T* a, const T* b, T* r );
    static void add( const T* a, T b, T* r );
    static void sub( const T* a, const T* b, T* r );
    static void sub( const T* a, T b, T* r );
    static void sub( T a, const T* b, T* r );
    static void mul( const T* a, const T* b, T* r );
    static void mul( const T* a, T b, T* r );
    static void div( const T* a, const T* b, T* r );
    static void div( const T* a, T b, T* r );
    
    static bool equal( const T* a, const T* b );
    
};


// implement + - * / and equal
template<class T, unsigned nrows, unsigned ncols>
void
vnl_matrix_fixed<T,nrows,ncols>::add( const T* a, const T* b, T* r )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        *(r++) = *(a++) + *(b++);
}


template<class T, unsigned nrows, unsigned ncols>
void
vnl_matrix_fixed<T,nrows,ncols>::add( const T* a, T b, T* r )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        *(r++) = *(a++) + b;
}

template<class T, unsigned nrows, unsigned ncols>
void
vnl_matrix_fixed<T,nrows,ncols>::sub( const T* a, const T* b, T* r )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        *(r++) = *(a++) - *(b++);
}

template<class T, unsigned nrows, unsigned ncols>
void
vnl_matrix_fixed<T,nrows,ncols>::sub( const T* a, T b, T* r )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        *(r++) = *(a++) - b;
}

template<class T, unsigned nrows, unsigned ncols>
void
vnl_matrix_fixed<T,nrows,ncols>::sub( T a, const T* b, T* r )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        *(r++) = a - *(b++);
}

template<class T, unsigned nrows, unsigned ncols>
void
vnl_matrix_fixed<T,nrows,ncols>::mul( const T* a, const T* b, T* r )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        *(r++) = *(a++) * *(b++);
}

template<class T, unsigned nrows, unsigned ncols>
void
vnl_matrix_fixed<T,nrows,ncols>::mul( const T* a, T b, T* r )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        *(r++) = *(a++) * b;
}

template<class T, unsigned nrows, unsigned ncols>
void
vnl_matrix_fixed<T,nrows,ncols>::div( const T* a, const T* b, T* r )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        *(r++) = *(a++) / *(b++);
}

template<class T, unsigned nrows, unsigned ncols>
void
vnl_matrix_fixed<T,nrows,ncols>::div( const T* a, T b, T* r )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        *(r++) = *(a++) / b;
}

template<class T, unsigned nrows, unsigned ncols>
bool
vnl_matrix_fixed<T,nrows,ncols>::equal( const T* a, const T* b )
{
    unsigned int count = nrows*ncols;
    while ( count-- )
        if ( *(a++) != *(b++) )  return false;
    return true;
}



// --- Matrix-scalar -------------------------------------------------------------

template <class T, unsigned m, unsigned n>
inline
vnl_matrix_fixed<T,m,n> operator+( const vnl_matrix_fixed<T,m,n>& mat1, const vnl_matrix_fixed<T,m,n>& mat2 )
{
    vnl_matrix_fixed<T,m,n> r;
    vnl_matrix_fixed<T,m,n>::add( mat1.data(), mat2.data(), r.data() );
    return r;
}

template <class T, unsigned m, unsigned n>
inline
vnl_matrix_fixed<T,m,n> operator+( const vnl_matrix_fixed<T,m,n>& mat, T s )
{
    vnl_matrix_fixed<T,m,n> r;
    vnl_matrix_fixed<T,m,n>::add( mat.data(), s, r.data() );
    return r;
}

template <class T, unsigned m, unsigned n>
inline
vnl_matrix_fixed<T,m,n> operator+( const T& s,
                                  const vnl_matrix_fixed<T,m,n>& mat )
{
    vnl_matrix_fixed<T,m,n> r;
    vnl_matrix_fixed<T,m,n>::add( mat.data(), s, r.data() );
    return r;
}

template <class T, unsigned m, unsigned n>
inline
vnl_matrix_fixed<T,m,n> operator-( const vnl_matrix_fixed<T,m,n>& mat1, const vnl_matrix_fixed<T,m,n>& mat2 )
{
    vnl_matrix_fixed<T,m,n> r;
    vnl_matrix_fixed<T,m,n>::sub( mat1.data(), mat2.data(), r.data() );
    return r;
}

template <class T, unsigned m, unsigned n>
inline
vnl_matrix_fixed<T,m,n> operator-( const vnl_matrix_fixed<T,m,n>& mat, T s )
{
    vnl_matrix_fixed<T,m,n> r;
    vnl_matrix_fixed<T,m,n>::sub( mat.data(), s, r.data() );
    return r;
}

template <class T, unsigned m, unsigned n>
inline
vnl_matrix_fixed<T,m,n> operator-( const T& s,
                                  const vnl_matrix_fixed<T,m,n>& mat )
{
    vnl_matrix_fixed<T,m,n> r;
    vnl_matrix_fixed<T,m,n>::sub( s, mat.data(), r.data() );
    return r;
}


#endif /* vnl_matrix_fixed_h */
