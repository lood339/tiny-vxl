//
//  vnl_svd.h
//  vxl
//
//  Created by Jimmy Chen on 2020-06-23.
//  Copyright © 2020 Nowhere Planet. All rights reserved.
//

#ifndef vnl_svd_h
#define vnl_svd_h

#include <Eigen/Dense>

#include <vnl/vnl_numeric_traits.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>

template <typename T> class vnl_svd;

template <typename T>
class vnl_svd
{
public:
    //: The singular values of a matrix of complex<T> are of type T, not complex<T>
    using singval_t = typename vnl_numeric_traits<T>::abs_t;
    
    //:
    // Construct a vnl_svd<T> object from $m \times n$ matrix $M$.  The
    // vnl_svd<T> object contains matrices $U$, $W$, $V$ such that
    // $U W V^\top = M$.
    //
    // Uses linpack routine DSVDC to calculate an ``economy-size'' SVD
    // where the returned $U$ is the same size as $M$, while $W$
    // and $V$ are both $n \times n$.  This is efficient for
    // large rectangular solves where $m > n$, typical in least squares.
    //
    // The optional argument zero_out_tol is used to mark the zero singular
    // values: If nonnegative, any s.v. smaller than zero_out_tol in
    // absolute value is set to zero.  If zero_out_tol is negative, the
    // zeroing is relative to |zero_out_tol| * sigma_max();
    
    vnl_svd(vnl_matrix<T> const &M, double zero_out_tol = 0.0);
   
    
    
private:
    
    Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> svd_;
    
    int m_, n_;              // Size of M, local cache.
    vnl_matrix<T> U_;        // Columns Ui are basis for range of M for Wi != 0
    vnl_diag_matrix<singval_t> W_;// Singular values, sorted in decreasing order
    vnl_diag_matrix<singval_t> Winverse_;
    vnl_matrix<T> V_;       // Columns Vi are basis for nullspace of M for Wi = 0
    unsigned rank_;
    bool have_max_;
    singval_t max_;
    bool have_min_;
    singval_t min_;
    double last_tol_;
    bool valid_;        // false if the NETLIB call failed.
    
    // Disallow assignment.
    vnl_svd(vnl_svd<T> const &) { }
    vnl_svd<T>& operator=(vnl_svd<T> const &) { return *this; }
};


// implementation
template<typename T>
vnl_svd<T>::vnl_svd(vnl_matrix<T> const &M, double zero_out_tol):
m_(M.rows()), n_(M.cols()), U_(m_, n_), W_(n_), Winverse_(n_), V_(n_, n_)
{
    
}



#endif /* vnl_svd_h */
