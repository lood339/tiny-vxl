// This is core/vnl/algo/vnl_levenberg_marquardt.h
#ifndef vnl_levenberg_marquardt_h_
#define vnl_levenberg_marquardt_h_
//:
// \file
// \brief Levenberg Marquardt nonlinear least squares
// \author Andrew W. Fitzgibbon, Oxford RRG
// \date   31 Aug 96
//
// \verbatim
// Modifications
//  AGAP 160701 Some comments. Changed minimize to call the correct minimization
//              routine.
//  RWMC 001097 Added verbose flag to get rid of all that blathering.
//  AWF  151197 Added trace flag to increase blather.
//   Feb.2002 - Peter Vanroose - brief doxygen comment placed on single line
// \endverbatim
//

#include <iosfwd>
#include <iostream>
#include <cassert>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_nonlinear_minimizer.h>
#include <vnl/vnl_least_squares_function.h>
#include <vnl/algo/vnl_algo_export.h>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

class vnl_least_squares_function;

//: Levenberg Marquardt nonlinear least squares
//  vnl_levenberg_marquardt is an interface to the MINPACK routine lmdif,
//  and implements Levenberg Marquardt nonlinear fitting.  The function
//  to be minimized is passed as a vnl_least_squares_function object, which
//  may or may not wish to provide derivatives.  If derivatives are not
//  supplied, they are calculated by forward differencing, which costs
//  one function evaluation per dimension, but is perfectly accurate.
//  (See Hartley in ``Applications of Invariance in Computer Vision''
//  for example).

template <typename T>
class VNL_ALGO_EXPORT vnl_levenberg_marquardt : public vnl_nonlinear_minimizer
{
    // T is derived from vnl_least_squares_function
    using NumericalDiffFunctor = Eigen::NumericalDiff<T>;
    using JacobianType = vnl_least_squares_function::JacobianType;
    
protected:
    // from Eigen
    NumericalDiffFunctor num_dif_functor_;
    Eigen::LevenbergMarquardt<NumericalDiffFunctor, double> lmq_;
 public:

    //: Initialize with the function object that is to be minimized.
    vnl_levenberg_marquardt(T& f):
    num_dif_functor_(std::ref(f)), // @why need std::ref
    lmq_(num_dif_functor_)
    {
        //assert(std::is_base_of<vnl_least_squares_function, T>::value);//
        init(&f);
    }
   

    ~vnl_levenberg_marquardt() override;

    //: Minimize the function supplied in the constructor until convergence or failure.
    //  On return, x is such that f(x) is the lowest value achieved.
    //  Returns true for convergence, false for failure.
    //  Does not use the gradient even if the cost function provides one.
    bool minimize_without_gradient(vnl_vector<double>& x);

    //: Minimize the function supplied in the constructor until convergence or failure.
    //  On return, x is such that f(x) is the lowest value achieved.
    //  Returns true for convergence, false for failure.
    //  The cost function must provide a gradient.
    bool minimize_using_gradient  (vnl_vector<double>& x);

    //: Calls minimize_using_gradient() or minimize_without_gradient(),
    // depending on whether the cost function provides a gradient.
    bool minimize(vnl_vector<double>& x);
    bool minimize(vnl_vector_fixed<double,1>& x) { vnl_vector<double> y=x.extract(1); bool b=minimize(y); x=y; return b; }
    bool minimize(vnl_vector_fixed<double,2>& x) { vnl_vector<double> y=x.extract(2); bool b=minimize(y); x=y; return b; }
    bool minimize(vnl_vector_fixed<double,3>& x) { vnl_vector<double> y=x.extract(3); bool b=minimize(y); x=y; return b; }
    bool minimize(vnl_vector_fixed<double,4>& x) { vnl_vector<double> y=x.extract(4); bool b=minimize(y); x=y; return b; }

    // Coping with failure-------------------------------------------------------

    //: Provide an ASCII diagnosis of the last minimization on std::ostream.
    void diagnose_outcome(/*std::cerr*/) const;
    void diagnose_outcome(std::ostream&) const;

    //: Return J'*J computed at last minimum.
    //  it is an approximation of inverse of covariance
    vnl_matrix<double> const& get_JtJ();

 protected:
    // From Eigen
   

    T* f_;
    vnl_matrix<double> fdjac_; // Computed during lmdif/lmder
    vnl_vector<long>    ipvt_; // Also computed, both needed to get J'*J at end.

    vnl_matrix<double> inv_covar_;
    bool set_covariance_; // Set if covariance_ holds J'*J

    void init(T* f);

    /*
  // Communication with callback
  static void lmdif_lsqfun(long* m, long* n, double* x,
                           double* fx, long* iflag, void* userdata);
  static void lmder_lsqfun(long* m, long* n, double* x,
                           double* fx, double* fJ, long*, long* iflag,
                           void* userdata);
     */
};

//: Find minimum of "f", starting at "initial_estimate", and return.
vnl_vector<double> vnl_levenberg_marquardt_minimize(vnl_least_squares_function& f,
                                                    vnl_vector<double> const& initial_estimate);

// copy from .cpp
// see header
/*
vnl_vector<double>
vnl_levenberg_marquardt_minimize(vnl_least_squares_function & f, vnl_vector<double> const & initial_estimate)
{
    vnl_vector<double> x = initial_estimate;
    vnl_levenberg_marquardt lm(f);
    lm.minimize(x);
    return x;
}
 */

// ctor
template <typename T>
void vnl_levenberg_marquardt<T>::init(T * f)
{
   // lmq_ = Eigen::LevenbergMarquardt<Eigen::NumericalDiff<vnl_least_squares_function>, double>(*f);
    
    f_ = f;
    
    // If changing these defaults, check the help comments in vnl_levenberg_marquardt.h,
    // and MAKE SURE they're consistent.
    xtol = 1e-8;                                // Termination tolerance on X (solution vector)
    maxfev = 400 * f->get_number_of_unknowns(); // Termination maximum number of iterations.
    ftol = xtol * 0.01;                         // Termination tolerance on F (sum of squared residuals)
    gtol = 1e-5;                                // Termination tolerance on Grad(F)' * F = 0
    epsfcn = xtol * 0.001;                      // Step length for FD Jacobian
    
    unsigned int m = f_->get_number_of_residuals(); // I  Number of residuals, must be > #unknowns
    unsigned int n = f_->get_number_of_unknowns();  // I  Number of unknowns
    
    set_covariance_ = false;
    fdjac_.set_size(n, m);
    fdjac_.fill(0.0);
    ipvt_.set_size(n);
    ipvt_.fill(0);
    inv_covar_.set_size(n, n);
    inv_covar_.fill(0.0);
    // fdjac_ = new vnl_matrix<double>(n,m);
    // ipvt_ = new vnl_vector<int>(n);
    // covariance_ = new vnl_matrix<double>(n,n);
    
    
}

template <typename T>
vnl_levenberg_marquardt<T>::~vnl_levenberg_marquardt()
{
    // delete covariance_;
    // delete fdjac_;
    // delete ipvt_;
}

//--------------------------------------------------------------------------------

/*
void
vnl_levenberg_marquardt::lmdif_lsqfun(long * n,     // I   Number of residuals
                                      long * p,     // I   Number of unknowns
                                      double * x,   // I   Solution vector, size n
                                      double * fx,  // O   Residual vector f(x)
                                      long * iflag, // IO  0 ==> print, -1 ==> terminate
                                      void * userdata)
{
    auto * self = static_cast<vnl_levenberg_marquardt *>(userdata);
    vnl_least_squares_function * f = self->f_;
    assert(*p == (int)f->get_number_of_unknowns());
    assert(*n == (int)f->get_number_of_residuals());
    vnl_vector_ref<double> ref_x(*p, const_cast<double *>(x));
    vnl_vector_ref<double> ref_fx(*n, fx);
    
    if (*iflag == 0)
    {
        if (self->trace)
            std::cerr << "lmdif: iter " << self->num_iterations_ << " err [" << x[0] << ", " << x[1] << ", " << x[2] << ", "
            << x[3] << ", " << x[4] << ", ... ] = " << ref_fx.magnitude() << '\n';
        
        f->trace(self->num_iterations_, ref_x, ref_fx);
        ++(self->num_iterations_);
    }
    else
    {
        f->f(ref_x, ref_fx);
    }
    
    if (self->start_error_ == 0)
        self->start_error_ = ref_fx.rms();
    
    if (f->failure)
    {
        f->clear_failure();
        *iflag = -1; // fsm
    }
}
*/

// This function shouldn't be inlined, because (1) modification of the
// body will not cause the timestamp on the header to change, and so
// others will not be forced to recompile, and (2) the cost of making
// one more function call is far, far less than the cost of actually
// performing the minimisation, so the inline doesn't gain you
// anything.
//
template <typename T>
bool
vnl_levenberg_marquardt<T>::minimize(vnl_vector<double> & x)
{
    if (f_->has_gradient())
        return minimize_using_gradient(x);
    else
        return minimize_without_gradient(x);
}

//
template <typename T>
bool
vnl_levenberg_marquardt<T>::minimize_without_gradient(vnl_vector<double> & x)
{
    
    // fsm
    if (f_->has_gradient())
    {
        std::cerr << __FILE__ " : WARNING. calling minimize_without_gradient(), but f_ has gradient.\n";
    }
    
    // e04fcf
    long m = f_->get_number_of_residuals(); // I  Number of residuals, must be > #unknowns
    long n = f_->get_number_of_unknowns();  // I  Number of unknowns
    
    if (m < n)
    {
        std::cerr << "vnl_levenberg_marquardt: Number of unknowns(" << n << ") greater than number of data (" << m << ")\n";
        failure_code_ = ERROR_DODGY_INPUT;
        return false;
    }
     
    
    if (int(x.size()) != n)
    {
        std::cerr << "vnl_levenberg_marquardt: Input vector length (" << x.size() << ") not equal to num unknowns (" << n
        << ")\n";
        failure_code_ = ERROR_DODGY_INPUT;
        return false;
    }
    
    vnl_least_squares_function::InputType v_x = x;
    
    Eigen::LevenbergMarquardtSpace::Status status = lmq_.minimize(v_x);
    std::cout<<"Debug ...... LMQ status: "<<status<<std::endl;
    x = v_x;
    return true;

}

//--------------------------------------------------------------------------------
template <typename T>
bool
vnl_levenberg_marquardt<T>::minimize_using_gradient(vnl_vector<double> & x)
{
    // fsm
    if (!f_->has_gradient())
    {
        std::cerr << __FILE__ ": called method minimize_using_gradient(), but f_ has no gradient.\n";
        return false;
    }
    
    long m = f_->get_number_of_residuals(); // I  Number of residuals, must be > #unknowns
    long n = f_->get_number_of_unknowns();  // I  Number of unknowns
    
    if (m < n)
    {
        std::cerr << __FILE__ ": Number of unknowns(" << n << ") greater than number of data (" << m << ")\n";
        failure_code_ = ERROR_DODGY_INPUT;
        return false;
    }
    vnl_least_squares_function::InputType v_x = x;
    Eigen::LevenbergMarquardtSpace::Status status = lmq_.minimize(v_x);
    std::cout<<"Debug ...... LMQ status: "<<status<<std::endl;
    x = v_x;
    return true;
    
    /*
    vnl_vector<double> fx(m, 0.0); // W m   Explicitly set target to 0.0
    
    num_iterations_ = 0;
    set_covariance_ = false;
    long info;
    start_error_ = 0; // Set to 0 so first call to lmder_lsqfun will know to set it.
    
    
    double factor = 100;
    long nprint = 1;
    long mode = 1, nfev, njev;
    
    vnl_vector<double> diag(n, 0);
    vnl_vector<double> qtf(n, 0);
    vnl_vector<double> wa1(n, 0);
    vnl_vector<double> wa2(n, 0);
    vnl_vector<double> wa3(n, 0);
    vnl_vector<double> wa4(m, 0);
    */
    
    
    /*
    v3p_netlib_lmder_(lmder_lsqfun,
                      &m,
                      &n,
                      x.data_block(),
                      fx.data_block(),
                      fdjac_.data_block(),
                      &m,
                      &ftol,
                      &xtol,
                      &gtol,
                      &maxfev,
                      diag.data_block(),
                      &mode,
                      &factor,
                      &nprint,
                      &info,
                      &nfev,
                      &njev,
                      ipvt_.data_block(),
                      qtf.data_block(),
                      wa1.data_block(),
                      wa2.data_block(),
                      wa3.data_block(),
                      wa4.data_block(),
                      this);
     */
    
    
    num_evaluations_ = num_iterations_; // for lmder, these are the same.
    /*
    if (info < 0)
        info = ERROR_FAILURE;
    failure_code_ = (ReturnCodes)info;
    end_error_ = fx.rms();
    
    // Translate status code
    switch (failure_code_)
    {
        case 1: // ftol
        case 2: // xtol
        case 3: // both
        case 4: // gtol
            return true;
        default:
            return false;
    }
     */
}

//--------------------------------------------------------------------------------

template <typename T>
void
vnl_levenberg_marquardt<T>::diagnose_outcome() const
{
    diagnose_outcome(std::cerr);
}

// fsm: should this function be a method on vnl_nonlinear_minimizer?
// if not, the return codes should be moved into LM.
template <typename T>
void
vnl_levenberg_marquardt<T>::diagnose_outcome(std::ostream & s) const
{
#define whoami "vnl_levenberg_marquardt"
    // if (!verbose_) return;
    switch (failure_code_)
    {
            //  case -1:
            // have already warned.
            //    return;
        case ERROR_FAILURE:
            s << (whoami ": OIOIOI -- failure in leastsquares function\n");
            break;
        case ERROR_DODGY_INPUT:
            s << (whoami ": OIOIOI -- lmdif dodgy input\n");
            break;
        case CONVERGED_FTOL: // ftol
            s << (whoami ": converged to ftol\n");
            break;
        case CONVERGED_XTOL: // xtol
            s << (whoami ": converged to xtol\n");
            break;
        case CONVERGED_XFTOL: // both
            s << (whoami ": converged nicely\n");
            break;
        case CONVERGED_GTOL:
            s << (whoami ": converged via gtol\n");
            break;
        case TOO_MANY_ITERATIONS:
            s << (whoami ": too many iterations\n");
            break;
        case FAILED_FTOL_TOO_SMALL:
            s << (whoami ": ftol is too small. no further reduction in the sum of squares is possible.\n");
            break;
        case FAILED_XTOL_TOO_SMALL:
            s << (whoami ": xtol is too small. no further improvement in the approximate solution x is possible.\n");
            break;
        case FAILED_GTOL_TOO_SMALL:
            s << (whoami ": gtol is too small. Fx is orthogonal to the columns of the jacobian to machine precision.\n");
            break;
        default:
            s << (whoami ": OIOIOI: unkown info code from lmder.\n");
            break;
    }
    unsigned int m = f_->get_number_of_residuals();
    s << whoami ": " << num_iterations_ << " iterations, " << num_evaluations_ << " evaluations, " << m
    << " residuals.  RMS error start/end " << get_start_error() << '/' << get_end_error() << std::endl;
#undef whoami
}

// fjac is an output m by n array. the upper n by n submatrix
//         of fjac contains an upper triangular matrix r with
//         diagonal elements of nonincreasing magnitude such that
//
//                t     t           t
//               p *(jac *jac)*p = r *r,
//
//         where p is a permutation matrix and jac is the final
//         calculated jacobian. column j of p is column ipvt(j)
//         (see below) of the identity matrix. the lower trapezoidal
//         part of fjac contains information generated during
//         the computation of r.

// fdjac is target m*n

//: Get INVERSE of covariance at last minimum.
// Code thanks to Joss Knight (joss@robots.ox.ac.uk)
template <typename T>
vnl_matrix<double> const &
vnl_levenberg_marquardt<T>::get_JtJ()
{
    /*
    const int m = m_;
    const int n = n_;
    // check norm
    double fnorm = lm.fvec.blueNorm();
    double covfac = fnorm*fnorm/(m-n);
    
    // check covariance
    
    internal::covar(lm.fjac, lm.permutation.indices()); // TODO : move this as a function of lm
    
   
    
    MatrixXd cov;
    cov =  covfac*lm.fjac.topLeftCorner<n,n>();
    VERIFY_IS_APPROX( cov, cov_ref);
    */
    std::cerr<<"Error  vnl_levenberg_marquardt<T>::get_JtJ() is not supported \n";
    JacobianType J = lmq_.fjac;
    JacobianType JtJ = J.transpose()*J;
    inv_covar_ = JtJ;
    return inv_covar_;
    
    /*
    if (!set_covariance_)
    {
        std::cerr << __FILE__ ": get_covariance() not confirmed tested  yet\n";
        unsigned int n = fdjac_.rows();
        
        // matrix in FORTRAN is column-wise.
        // transpose it to get C style order
        vnl_matrix<double> r = fdjac_.extract(n, n).transpose();
        // r is upper triangular matrix according to documentation.
        // But the lower part has non-zeros somehow.
        // clear the lower part
        for (unsigned int i = 0; i < n; ++i)
            for (unsigned int j = 0; j < i; ++j)
                r(i, j) = 0.0;
        
        // compute r^T * r
        vnl_matrix<double> rtr;
        vnl_fastops::AtA(rtr, r);
        vnl_matrix<double> rtrpt(n, n);
        
        // Permute. First order columns.
        // Note, *ipvt_ contains 1 to n, not 0 to n-1
        vnl_vector<int> jpvt(n);
        for (unsigned int j = 0; j < n; ++j)
        {
            unsigned int i = 0;
            for (; i < n; i++)
            {
                if (ipvt_[i] == (int)j + 1)
                {
                    jpvt(j) = i;
                    break;
                }
            }
            rtrpt.set_column(j, rtr.get_column(i));
        }
        
        // Now order rows
        for (unsigned int j = 0; j < n; ++j)
        {
            inv_covar_.set_row(j, rtrpt.get_row(jpvt(j)));
        }
        
        set_covariance_ = true;
    }
     return inv_covar_;
     */
}

#endif // vnl_levenberg_marquardt_h_
