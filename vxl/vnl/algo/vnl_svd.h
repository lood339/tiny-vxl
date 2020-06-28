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
    using RowMajorMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
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
    
    virtual ~vnl_svd() = default;
    
    // Data Access---------------------------------------------------------------
    
    //: find weights below threshold tol, zero them out, and update W_ and Winverse_
    void zero_out_absolute(double tol = 1e-8); //sqrt(machine epsilon)
    
    //: find weights below tol*max(w) and zero them out
    void zero_out_relative(double tol = 1e-8); //sqrt(machine epsilon)
    int             singularities () const { return W_.rows() - rank(); }
    unsigned int    rank () const { return rank_; }
    singval_t       well_condition () const { return sigma_min()/sigma_max(); }
    
    //: Calculate determinant as product of diagonals in W.
    singval_t       determinant_magnitude () const;
    singval_t       norm() const;
    
    //: Return the matrix U.
    vnl_matrix<T>      & U()       { return U_; }
    
    //: Return the matrix U.
    vnl_matrix<T> const& U() const { return U_; }
    
    //: Return the matrix U's (i,j)th entry (to avoid svd.U()(i,j); ).
    T U(int i, int j) const { return U_(i,j); }
    
    //: Get at DiagMatrix (q.v.) of singular values, sorted from largest to smallest
    vnl_diag_matrix<singval_t>       & W()             { return W_; }
    
    //: Get at DiagMatrix (q.v.) of singular values, sorted from largest to smallest
    vnl_diag_matrix<singval_t> const & W() const       { return W_; }
    vnl_diag_matrix<singval_t>       & Winverse()      { return Winverse_; }
    vnl_diag_matrix<singval_t> const & Winverse() const { return Winverse_; }
    singval_t                   & W(int i, int j) { return W_(i,j); }
    singval_t                   & W(int i)        { return W_(i,i); }
    singval_t     sigma_max() const { return W_(0,0); }       // largest
    singval_t     sigma_min() const { return W_(n_-1,n_-1); } // smallest
    
    //: Return the matrix V.
    vnl_matrix<T>      & V()       { return V_; }
    
    //: Return the matrix V.
    vnl_matrix<T> const& V() const { return V_; }
    
    //: Return the matrix V's (i,j)th entry (to avoid svd.V()(i,j); ).
    T V(int i, int j) const { return V_(i,j); }
    
    
    //:
    inline vnl_matrix<T> inverse () const { return pinverse(); }
    
    //: pseudo-inverse (for non-square matrix) of desired rank.
    vnl_matrix<T> pinverse (unsigned int rank = ~0u) const; // ~0u == (unsigned int)-1
    
    
    //: Calculate inverse of transpose, using desired rank.
    vnl_matrix<T> tinverse (unsigned int rank = ~0u) const; // ~0u == (unsigned int)-1
    
    //: Recompose SVD to U*W*V', using desired rank.
    vnl_matrix<T> recompose (unsigned int rank = ~0u) const; // ~0u == (unsigned int)-1
    
    /*
    //: Solve the matrix equation M X = B, returning X
    vnl_matrix<T> solve (vnl_matrix<T> const& B) const;
    */
    //: Solve the matrix-vector system M x = y, returning x.
    vnl_vector<T> solve (vnl_vector<T> const& y) const;
    /*
    void          solve (T const *rhs, T *lhs) const; // min ||A*lhs - rhs||
    
    //: Solve the matrix-vector system M x = y.
    // Assuming that the singular values W have been preinverted by the caller.
    void solve_preinverted(vnl_vector<T> const& rhs, vnl_vector<T>* out) const;
    
    //: Return N such that M * N = 0
    vnl_matrix<T> nullspace() const;
    
    //: Return N such that M' * N = 0
    vnl_matrix<T> left_nullspace() const;
    
    //: Return N such that M * N = 0
    vnl_matrix<T> nullspace(int required_nullspace_dimension) const;
    
    //: Implementation to be done yet; currently returns left_nullspace(). - PVR.
    vnl_matrix<T> left_nullspace(int required_nullspace_dimension) const;
    
    //: Return the rightmost column of V.
    //  Does not check to see whether or not the matrix actually was rank-deficient -
    // the caller is assumed to have examined W and decided that to his or her satisfaction.
    vnl_vector<T> nullvector() const;
    
    //: Return the rightmost column of U.
    //  Does not check to see whether or not the matrix actually was rank-deficient.
    vnl_vector<T> left_nullvector() const;
    
    bool valid() const { return valid_; }
     */
    
private:
    
    Eigen::JacobiSVD<RowMajorMatrix> svd_;
    
    int m_, n_;              // Size of M, local cache.
    vnl_matrix<T> U_;        // Columns Ui are basis for range of M for Wi != 0
    vnl_diag_matrix<singval_t> W_;// Singular values, sorted in decreasing order @todo it is not a diagonal matrix
    vnl_diag_matrix<singval_t> Winverse_; // @todo it is not a diagonal matrixs
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
m_(M.rows()), n_(M.cols()),
U_(m_, n_), // why is mxn, but not mxm
W_(n_), Winverse_(n_), V_(n_, n_)
{
    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // copy data from vnl_matrix to
    //
    RowMajorMatrix M_copy = M;
    svd_ = Eigen::JacobiSVD<RowMajorMatrix>(M_copy, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    /*
    if (vnl_svd_test_heavily)
    {
        // Test that recomposed matrix == M
        typedef typename vnl_numeric_traits<T>::abs_t abs_t;
        abs_t recomposition_residual = std::abs((recompose() - M).fro_norm());
        abs_t n = std::abs(M.fro_norm());
        abs_t thresh = abs_t(m_) * abs_t(vnl_math::eps) * n;
        if (recomposition_residual > thresh)
        {
            std::cerr << "vnl_svd<T>::vnl_svd<T>() -- Warning, recomposition_residual = "
            << recomposition_residual << std::endl
            << "fro_norm(M) = " << n << std::endl
            << "eps*fro_norm(M) = " << thresh << std::endl
            << "Press return to continue\n";
            char x;
            std::cin.get(&x, 1, '\n');
        }
    }
     */
    
    svd_.setThreshold(zero_out_tol);
    
    // copy U S V values
    const auto& U = svd_.matrixU();
    assert(U.rows() <= U_.rows());
    assert(U.cols() <= U_.cols());
    for(int i = 0; i<U.rows(); ++i) {
        for(int j = 0; j<U.cols(); ++j) {
            U_(i, j) = U(i, j);
        }
    }
    
    const auto& S = svd_.singularValues();
    W_.fill(0);
    assert(S.size() <= W_.rows());
    for(int i = 0; i<S.size(); ++i) {
        W_(i, i) = S[i];
    }
    
    const auto& V = svd_.matrixV();
    assert(V.rows() <= V_.rows());
    assert(V.cols() <= V_.cols());
    
    V_.fill(0);
    for(int i = 0; i<V.rows(); ++i) {
        for(int j = 0; j<V.cols(); ++j) {
            V_(i, j) = V(i, j);
        }
    }
    
    
    //std::cout<<"S singular value size: "<<S.size()<<std::endl;
    //std::cout<<"S: "<<S<<std::endl;
    //std::cout<<"W: row "<<W_.rows()<<std::endl;
    //std::cout<<"W from svd:"<<svd_.singularValues().asDiagonal()<<std::endl;
    
    if (zero_out_tol >= 0) {
        // Zero out small sv's and update rank count.
        zero_out_absolute(double(+zero_out_tol));
    }
    else {
        // negative tolerance implies relative to max elt.
        zero_out_relative(double(-zero_out_tol));
    }
    std::cout<<"rank: {}"<<svd_.rank()<<std::endl;
    assert(svd_.rank() == rank_);
}

template <typename T>
void vnl_svd<T>::zero_out_absolute(double tol)
{
    last_tol_ = tol;
    rank_ = W_.rows();
    for(unsigned k = 0; k<W_.rows(); ++k) {
        const singval_t weight = W_(k, k);
        if(std::abs(weight) <= tol) {
            Winverse_(k, k) = 0;
            W_(k, k) = 0;
            --rank_;
        }
        else {
            Winverse_(k, k) = singval_t(1.0)/weight;
        }
    }
}

//: find weights below tol*max(w) and zero them out
template <typename T>
void vnl_svd<T>::zero_out_relative(double tol) // sqrt(machine epsilon)
{
    zero_out_absolute(tol * std::abs(sigma_max()));
}

//: Calculate pseudo-inverse.
template <typename T>
vnl_matrix<T> vnl_svd<T>::pinverse(unsigned int rnk) const
{
    if (rnk > rank_) rnk=rank_;
    vnl_matrix<T> W_inverse(Winverse_.rows(), Winverse_.cols());
    W_inverse.fill(T(0));
    for (unsigned int i=0;i<rnk;++i)
        W_inverse(i,i)=Winverse_(i,i);
    
    return V_ * W_inverse * U_.conjugate_transpose();
}

//: Calculate inverse of transpose, using desired rank.
template <typename T>
vnl_matrix<T> vnl_svd<T>::tinverse (unsigned int rnk) const // ~0u == (unsigned int)-1
{
    if (rnk > rank_) rnk=rank_;
    vnl_matrix<T> W_inverse(Winverse_.rows(),Winverse_.columns());
    W_inverse.fill(T(0));
    for (unsigned int i=0;i<rnk;++i)
        W_inverse(i,i)=Winverse_(i,i);
    
    return U_ * W_inverse * V_.conjugate_transpose();
    
}

//: Recompose SVD to U*W*V', using desired rank.
template <typename T>
vnl_matrix<T> vnl_svd<T>::recompose (unsigned int rnk) const // ~0u == (unsigned int)-1
{
    if (rnk > rank_) rnk=rank_;
    vnl_matrix<T> Wmatr(W_.rows(), W_.cols());
    Wmatr.fill(T(0));
    for (unsigned int i=0;i<rnk;++i)
        Wmatr(i,i)=W_(i,i);
    
    return U_*Wmatr*V_.conjugate_transpose();
}

//: Solve the matrix-vector system M x = y, returning x.
template <class T>
vnl_vector<T> vnl_svd<T>::solve(vnl_vector<T> const& y)  const
{
    // fsm sanity check :
    if (y.size() != U_.rows())
    {
        std::cerr << __FILE__ << ": size of rhs is incompatible with no. of rows in U_\n"
        << "y =" << y << '\n'
        << "m_=" << m_ << '\n'
        << "n_=" << n_ << '\n'
        << "U_=\n" << U_
        << "V_=\n" << V_
        << "W_=\n" << W_;
    }
    
    vnl_vector<T> x(V_.rows());                   // Solution matrix.
    if (U_.rows() < U_.columns()) {               // Augment y with extra rows of
        vnl_vector<T> yy(U_.rows(), T(0));          // zeros, so that it matches
        if (yy.size()<y.size()) { // fsm
            std::cerr << "yy=" << yy << std::endl
            << "y =" << y  << std::endl;
            // the update() call on the next line will abort...
        }
        yy.update(y);                               // cols of u.transpose.
        x = U_.conjugate_transpose() * yy;
    }
    else
        x = U_.conjugate_transpose() * y;
    
    for (unsigned i = 0; i < x.size(); i++) {        // multiply with diagonal 1/W
        T weight = W_(i, i), zero_(0);
        if (weight != zero_)
            x[i] /= weight;
        else
            x[i] = zero_;
    }
    return V_ * x;                                // premultiply with v.
}



#endif /* vnl_svd_h */
