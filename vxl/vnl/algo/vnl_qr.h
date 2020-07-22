// This is core/vnl/algo/vnl_qr.h
#ifndef vnl_qr_h_
#define vnl_qr_h_
//:
// \file
// \brief Calculate inverse of a matrix using QR
// \author  Andrew W. Fitzgibbon, Oxford RRG
// \date   08 Dec 1996
//
// \verbatim
//  Modifications
//   081296 AWF Temporarily abandoned as I realized my problem was symmetric...
//   080697 AWF Recovered, implemented solve().
//   200897 AWF Added determinant().
//   071097 AWF Added Q(), R().
//   Christian Stoecklin, ETH Zurich, added QtB(v)
//   31-mar-2000 fsm: templated
//   28/03/2001 - dac (Manchester) - tidied up documentation
//   13 Jan.2003 - Peter Vanroose - added missing implementation for inverse(),
//                                tinverse(), solve(matrix), extract_q_and_r().
// \endverbatim

#include <iosfwd>
#include <iostream>
#include <complex>
#include <cassert>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_math.h>

#include <Eigen/Dense>

//#include <vnl/algo/vnl_algo_export.h>
//#include <vnl/vnl_complex.h>  // vnl_math::squared_magnitude()
//#include <vnl/vnl_matlab_print.h>
//#include <vnl/vnl_complex_traits.h>
//#include <vnl/algo/vnl_netlib.h> // dqrdc_(), dqrsl_()

//: Extract the Q*R decomposition of matrix M.
//  M: m x n, Q: mxm, R:mxn
//
//  The decomposition is stored in a compact and time-efficient
// packed form, which is most easily used via the "solve" and
// "determinant" methods.

template <class T>
class vnl_qr
{
    using RowMajorMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
 public:
    vnl_qr(vnl_matrix<T> const & M);
    ~vnl_qr();
    
    /*
    //: return the inverse matrix of M
    vnl_matrix<T> inverse () const;
    //: return the transpose of the inverse matrix of M
    vnl_matrix<T> tinverse () const;
     */
    //: return the original matrix M
    vnl_matrix<T> recompose () const;
     

    //: Solve equation M x = rhs for x using the computed decomposition.
    vnl_matrix<T> solve (const vnl_matrix<T>& rhs) const;
    //: Solve equation M x = rhs for x using the computed decomposition.
    vnl_vector<T> solve (const vnl_vector<T>& rhs) const;

    //: Return the determinant of M.  This is computed from M = Q R as follows:
    // |M| = |Q| |R|.
    // |R| is the product of the diagonal elements.
    // |Q| is (-1)^n as it is a product of Householder reflections.
    // So det = -prod(-r_ii).
    T determinant() const;

    //: Return residual vector d of M x = b -> d = Q'b
    //vnl_vector<T> QtB(const vnl_vector<T>& b) const;

    //: Unpack and return unitary part Q.
    vnl_matrix<T> const& Q() const;
    //: Unpack and return R.
    vnl_matrix<T> const& R() const;
    //void extract_q_and_r(vnl_matrix<T>* q, vnl_matrix<T>* r) const { *q = Q(); *r = R(); }
    
 private:
    
    int m_; // rows of input matrix
    int n_; // cols of input matrix
    
    mutable vnl_matrix<T>* Q_;
    mutable vnl_matrix<T>* R_;
    
    Eigen::HouseholderQR<RowMajorMatrix> qr_;

    // Disallow assignment.
    vnl_qr(const vnl_qr<T> &) { }
    vnl_qr<T>& operator=(const vnl_qr<T> &) { return *this; }
};

//: Compute determinant of matrix "M" using QR.
template <class T>
inline T vnl_qr_determinant(vnl_matrix<T> const& m)
{
  return vnl_qr<T>(m).determinant();
}

template <class T>
std::ostream& operator<<(std::ostream&, vnl_qr<T> const & qr);

// copy from .cpp
template <class T>
vnl_qr<T>::vnl_qr(vnl_matrix<T> const& M):
m_(M.rows()),
n_(M.columns()),
Q_(nullptr),
R_(nullptr)
{
    assert(! M.empty());
    
    RowMajorMatrix M_copy = M;
    qr_ = Eigen::HouseholderQR<RowMajorMatrix>(M_copy);
    
    //Eigen::ComputationInfo info = qr_.info();
    //if(info != Eigen::Success) {
    //    std::cerr <<"Error: QR decomposition failed\n";
    //}
    
    //std::cout<<"column permutation: \n"<<qr_.colsPermutation()<<std::endl;
    
}

template <class T>
vnl_qr<T>::~vnl_qr()
{
    delete Q_;
    delete R_;
}

//: Return the determinant of M.  This is computed from M = Q R as follows:
// |M| = |Q| |R|
// |R| is the product of the diagonal elements.
// |Q| is (-1)^n as it is a product of Householder reflections.
// So det = -prod(-r_ii).
template <class T>
T vnl_qr<T>::determinant() const
{
    vnl_matrix<T> R = this->R();
    
    int m = std::min(m_, n_);
    T det = R(0,0);

    for (int i = 1; i < m; ++i)
        det *= -R(i,i);

    return det;
    
}

//: Unpack and return unitary part Q.
template <class T>
vnl_matrix<T> const& vnl_qr<T>::Q() const
{
    if(!Q_) {
        int m = m_;
        RowMajorMatrix thinQ(m, m);
        thinQ.setIdentity();
        RowMajorMatrix Q = qr_.householderQ()*thinQ;
        assert(Q.rows() == m);
        assert(Q.cols() == m);
        
        Q_ = new vnl_matrix<T>(m, m);
        
        for(int i = 0; i<m; ++i) {
            for(int j = 0; j<m; ++j) {
                (*Q_)(i, j) = Q(i, j);
            }
        }
    }
    
    return *Q_;
}

//: Unpack and return R.
template <class T>
vnl_matrix<T> const& vnl_qr<T>::R() const
{
    if (!R_) {
        int m = m_;
        int n = n_;
        RowMajorMatrix R = qr_.matrixQR().template triangularView<Eigen::Upper>();
        assert(R.rows() == m);
        assert(R.cols() == n);
        
        R_ = new vnl_matrix<T>(m, n);
        for(int i = 0; i<m; ++i) {
            for(int j = 0; j<n; ++j) {
                (*R_)(i, j) = R(i, j);
            }
        }
    }
    
    return *R_;
}

template <class T>
vnl_matrix<T> vnl_qr<T>::recompose() const
{
    return Q() * R();
}

// JOB: ABCDE decimal
// A     B     C     D              E
// ---   ---   ---   ---            ---
// Qb    Q'b   x     norm(A*x - b)  A*x


//: Solve equation M x = b for x using the computed decomposition.
template <class T>
vnl_vector<T> vnl_qr<T>::solve(const vnl_vector<T>& b) const
{
    vnl_vector<T> x = qr_.solve(b);
    return x;
}

/*
//: Return residual vector d of M x = b -> d = Q'b
template <class T>
vnl_vector<T> vnl_qr<T>::QtB(const vnl_vector<T>& b) const
{
    long n = qrdc_out_.columns();
    long p = qrdc_out_.rows();
    const T* b_data = b.data_block();
    vnl_vector<T> Qt_B(n);
    
    // see comment above
    long JOB = 1000;
    
    long info = 0;
    vnl_linpack_qrsl(qrdc_out_.data_block(),
                     &n, &n, &p,
                     qraux_.data_block(),
                     b_data,
                     (T*)nullptr,               // A: Qb
                     Qt_B.data_block(),   // B: Q'b
                     (T*)nullptr,               // C: x
                     (T*)nullptr,               // D: residual
                     (T*)nullptr,               // E: Ax
                     &JOB,
                     &info);
    
    if (info > 0)
        std::cerr << __FILE__ ": vnl_qr<T>::QtB() -- matrix is rank-deficient by "
        << info << '\n';
    
    return Qt_B;
}

template <class T>
vnl_matrix<T> vnl_qr<T>::inverse() const
{
    unsigned int r = qrdc_out_.columns();
    assert(r > 0 && r == qrdc_out_.rows());
    vnl_matrix<T> inv(r,r);
    
    // Use solve() to compute the inverse matrix, using (00..010..00) as rhs
    vnl_vector<T> rhs(r,T(0));
    for (unsigned int i=0; i<r; ++i)
    {
        rhs(i) = T(1);
        vnl_vector<T> col = this->solve(rhs); // returns i-th column of inverse
        inv.set_column(i,col);
        rhs(i) = T(0);
    }
    return inv;
}

template <class T>
vnl_matrix<T> vnl_qr<T>::tinverse() const
{
    unsigned int r = qrdc_out_.columns();
    assert(r > 0 && r == qrdc_out_.rows());
    vnl_matrix<T> tinv(r,r);
    
    // Use solve() to compute the inverse matrix, using (00..010..00) as rhs
    vnl_vector<T> rhs(r,T(0));
    for (unsigned int i=0; i<r; ++i)
    {
        rhs(i) = T(1);
        vnl_vector<T> col = this->solve(rhs); // returns i-th column of inverse
        tinv.set_row(i,col);
        rhs(i) = T(0);
    }
    return tinv;
}

template <class T>
vnl_matrix<T> vnl_qr<T>::solve(vnl_matrix<T> const& rhs) const
{
    assert(rhs.rows() == qrdc_out_.columns()); // column-major storage
    int c = qrdc_out_.rows();
    int n = rhs.columns();
    vnl_matrix<T> result(c,n);
    
    for (int i=0; i<n; ++i)
    {
        vnl_vector<T> b = rhs.get_column(i);
        vnl_vector<T> col = this->solve(b); // returns i-th column of result
        result.set_column(i,col);
    }
    return result;
}
*/
#endif // vnl_qr_h_
