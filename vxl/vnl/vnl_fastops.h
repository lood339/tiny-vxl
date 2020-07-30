// This is core/vnl/vnl_fastops.h
#ifndef vnl_fastops_h_
#define vnl_fastops_h_
//:
//  \file
//  \brief Collection of C-style matrix functions
//  \author Andrew W. Fitzgibbon, Oxford RRG
//  \date   09 Dec 96
//
// \verbatim
//  Modifications
//   Feb.2002 -Peter Vanroose- brief doxygen comment placed on single line
//   Jun.2004 -Peter Vanroose- Added inc_X_by_ABt dec_X_by_AtB {inc,dec}_X_by_AB
//   Jun.2004 -Peter Vanroose- First step to migrate towards non-pointer args
//   Mar.2007 -Peter Vanroose- Commented deprecated versions of the functions
// \endverbatim

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_export.h>

//: Collection of C-style matrix functions for the most time-critical applications.
// In general, however one should consider using the vnl_transpose envelope-letter
// class to achieve the same results with about a 10% speed penalty.
class VNL_EXPORT vnl_fastops
{
 public:
  static void AtA(vnl_matrix<double>& out, const vnl_matrix<double>& A)
    {
        out = A.transpose()*A;
    }
  static void AB (vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        out = A*B;
    }
  static void AtB(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        out = A.transpose()*B;
    }
  static void AtB(vnl_vector<double>& out, const vnl_matrix<double>& A, const vnl_vector<double>& b)
    {
        out = A.transpose()*b;
    }
  static void Ab (vnl_vector<double>& out, const vnl_matrix<double>& A, const vnl_vector<double>& b)
    {
        out = out = A*b;
    }
  static void ABt(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        out = A*B.transpose();
    }

    static double btAb (const vnl_matrix<double>& A, const vnl_vector<double>& b);
    

  static void ABAt(vnl_matrix<double>& out, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        out = A*B*A.transpose();
    }

    // increase
  static void inc_X_by_AtA(vnl_matrix<double>& X, const vnl_matrix<double>& A)
    {
        X += A.transpose()*A;
    }
  static void inc_X_by_AB (vnl_matrix<double>& X, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        X += A*B;
    }
  static void inc_X_by_AtB(vnl_matrix<double>& X, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        X += A.transpose()*B;
    }
  static void inc_X_by_AtB(vnl_vector<double>& X, const vnl_matrix<double>& A, const vnl_vector<double>& b);
  static void inc_X_by_ABt(vnl_matrix<double>& X, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        X += A*B.transpose();
    }
  static void inc_X_by_ABAt(vnl_matrix<double>& X, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        X += A*B*A.transpose();
    }

    //decrease
  static void dec_X_by_AtA(vnl_matrix<double>& X, const vnl_matrix<double>& A)
    {
        X -= A.transpose()*A;
    }
  static void dec_X_by_AB (vnl_matrix<double>& X, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        X -= A*B;
    }
  static void dec_X_by_AtB(vnl_matrix<double>& X, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        X -= A.transpose()*B;
    }
  static void dec_X_by_AtB(vnl_vector<double>& X, const vnl_matrix<double>& A, const vnl_vector<double>& b);
    
  static void dec_X_by_ABt(vnl_matrix<double>& X, const vnl_matrix<double>& A, const vnl_matrix<double>& B)
    {
        X -= A*B.transpose();
    }

 private:
  // BLAS-like operations
  static double dot(const double* a, const double* b, unsigned int n);
};

// copy from .cpp

//: Compute $b^\top A b$ for vector b and matrix A
double
vnl_fastops::btAb(const vnl_matrix<double> & A, const vnl_vector<double> & b)
{
    const unsigned int m = A.rows();
    const unsigned int n = A.cols();
    const unsigned int l = b.size();
    
    // Verify matrices compatible
    if (m != l)
    {
        std::cerr << "vnl_fastops::btAb: argument sizes do not match: " << m << " != " << l << '\n';
        std::abort();
    }
    if (m != n)
    {
        std::cerr << "vnl_fastops::btAb: not a square matrix: " << m << " != " << n << '\n';
        std::abort();
    }
    
    
    //double const * const * a = A.data_array();
    double const * bb = b.data_block();
    
    double accum = 0;
    for (unsigned int i = 0; i < n; ++i)
        for (unsigned int j = 0; j < n; ++j)
        {
            accum += bb[j] * A[i][j] * bb[i];
        }
    return accum;
}

#endif // vnl_fastops_h_
