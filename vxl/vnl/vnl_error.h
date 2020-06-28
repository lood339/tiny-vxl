// This is core/vnl/vnl_error.h
#ifndef vnl_error_h_
#define vnl_error_h_
//:
//  \file
//  \author fsm

#include <iostream>
#include "vnl/vnl_export.h"

//: Raise exception for invalid index.
inline void
vnl_error_vector_index(char const * fcn, int index)
{
    // RAISE Error, SYM(vnl_error_vector), SYM(Invalid_Index),
    std::cerr << "vnl_error_vector_index:" << fcn << ": Invalid value " << index << " specified for index.\n";
    throw 0;
}

//: Raise exception for invalid dimension.
inline void
vnl_error_vector_dimension(char const * fcn, int l1, int l2)
{
    // RAISE Error, SYM(vnl_error_vector), SYM(Invalid_Dim),
    std::cerr << "vnl_error_vector_dimension:" << fcn << ": Dimensions [" << l1 << "] and [" << l2 << "] do not match.\n";
    throw 0;
}


//: Raise exception for using class objects, or chars in (...).
inline void
vnl_error_vector_va_arg(int n)
{
    // RAISE Error, SYM(vnl_error_vector), SYM(Invalid_Va_Arg),
    std::cerr << "vnl_error_vector_va_arg: Invalid type in ..."
    << " or wrong alignment with " << n << " bytes.\n";
    throw 0;
}

//--------------------------------------------------------------------------------

//: Raise exception for invalid row index.
inline void
vnl_error_matrix_row_index(char const * fcn, unsigned r)
{
    // RAISE Error, SYM(vnl_error_matrix), SYM(Invalid_Row),
    std::cerr << "vnl_error_matrix_row_index:" << fcn << ": Invalid value " << r << " specified for row.\n";
    throw 0;
}


//: Raise exception for invalid col index.
inline void
vnl_error_matrix_col_index(char const * fcn, unsigned c)
{
    // RAISE Error, SYM(vnl_error_matrix), SYM(Invalid_Col),
    std::cerr << "vnl_error_matrix_col_index:" << fcn << ": Invalid value " << c << " specified for column.\n";
    throw 0;
}

//: Raise exception for invalid dimensions.
inline void
vnl_error_matrix_dimension(char const * fcn, int r1, int c1, int r2, int c2)
{
    // RAISE Error, SYM(vnl_error_matrix), SYM(Invalid_Dim),
    std::cerr << "vnl_error_matrix_dimension:" << fcn << ": Dimensions [" << r1 << ',' << c1 << "] and [" << r2 << ','
    << c2 << "] do not match.\n";
    throw 0;
}


//: Raise exception for a nonsquare matrix.
inline void
vnl_error_matrix_nonsquare(char const * fcn)
{
    // RAISE Error, SYM(vnl_error_matrix), SYM(Invalid_Dim),
    std::cerr << "vnl_error_matrix_nonsquare:" << fcn << ": Matrix must be square.\n";
    throw 0;
}

//: Raise exception for using class objects, or chars in (...).
inline void
vnl_error_matrix_va_arg(int n)
{
    // RAISE Error, SYM(vnl_error_matrix), SYM(Invalid_Va_Arg),
    std::cerr << "vnl_error_matrix_va_arg: Invalid type in ..."
    << " or wrong alignment with " << n << " bytes.\n";
    throw 0;
}

#endif // vnl_error_h_
