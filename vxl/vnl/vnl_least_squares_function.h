// This is core/vnl/vnl_least_squares_function.h
#ifndef vnl_least_squares_function_h_
#define vnl_least_squares_function_h_
//:
// \file
// \brief Abstract base for minimising functions
// \author Andrew W. Fitzgibbon, Oxford RRG
// \date   31 Aug 96
//
// \verbatim
//  Modifications
//   280697 AWF Changed return type of f from double to void, as it wasn't used, and
//              people were going to extra trouble to compute it.
//   20 Apr 1999 FSM Added failure flag so that f() and grad() may signal failure to the caller.
//   23/3/01 LSB (Manchester) Tidied documentation
//   Feb.2002 - Peter Vanroose - brief doxygen comment placed on single line
// \endverbatim
//
// not used? #include <vcl_compiler.h>
#include <string>
#include <iostream>
#include <cassert>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_export.h>


//:  Abstract base for minimising functions.
//    vnl_least_squares_function is an abstract base for functions to be minimized
//    by an optimizer.  To define your own function to be minimized, subclass
//    from vnl_least_squares_function, and implement the pure virtual f (and
//    optionally grad_f).
//
//    Whether or not f ought to be const is a problem.  Clients might well
//    want to cache some information during the call, and if they're compute
//    objects, will almost certainly be writing to members during the
//    computation.  For the moment it's non-const, but we'll see...
class VNL_EXPORT vnl_least_squares_function
{
    using RowVectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor>;
 public:
  enum  UseGradient {
    no_gradient,
    use_gradient
  };
  bool failure;

  //: Construct vnl_least_squares_function.
  // Passing number of parameters (unknowns, domain dimension) and number of
  // residuals (range dimension).
  // The optional argument should be no_gradient if the gradf function has not
  // been implemented.  Default is use_gradient.
  vnl_least_squares_function(unsigned int number_of_unknowns,
                             unsigned int number_of_residuals,
                             UseGradient g = use_gradient)
  : failure(false), p_(number_of_unknowns), n_(number_of_residuals),
    use_gradient_(g == use_gradient)
  { dim_warning(p_,n_); }

  virtual ~vnl_least_squares_function() = default;

  // the virtuals may call this to signal a failure.
  void throw_failure() { failure = true; }
  void clear_failure() { failure = false; }
    
    // Interface for LevenbergMarquardt Functor
    int inputs() const {return p_;}
    int values() const {return n_;}
    
    
    int operator()(const RowVectorXd &x, RowVectorXd &fvec) const
    {
        
        //Eigen::Map<vnl_vector<double>> v_x(x.data(), x.size());
        //Eigen::Map<vnl_vector<double>&> v_fvec(fvec.data(), fvec.size());
        vnl_vector<double> v_x = x;
        vnl_vector<double> v_fvec = fvec;
        this->f(v_x, v_fvec);
        for(int i = 0; i<fvec.size(); ++i) {
            fvec[i] = v_fvec[i];
        }
     
        return 0;
    }
     /**/
    //int df(const RowVectorXd &x, RowVectorXd &fjac) const;
    

  //: The main function.
  //  Given the parameter vector x, compute the vector of residuals fx.
  //  Fx has been sized appropriately before the call.
  virtual void f(vnl_vector<double> const& x, vnl_vector<double>& fx) const = 0;

  //: Calculate the Jacobian, given the parameter vector x.
  virtual void gradf(vnl_vector<double> const& x, vnl_matrix<double>& jacobian);

  //: Use this to compute a finite-difference gradient other than lmdif
  void fdgradf(vnl_vector<double> const& x, vnl_matrix<double>& jacobian,
               double stepsize);

  //: Use this to compute a finite-forward-difference gradient other than lmdif
  // This takes about half as many estimates as fdgradf
  void ffdgradf(vnl_vector<double> const& x, vnl_matrix<double>& jacobian,
                double stepsize);

  //: Called after each LM iteration to print debugging etc.
  virtual void trace(int iteration,
                     vnl_vector<double> const& x,
                     vnl_vector<double> const& fx);

  //: Compute the rms error at x by calling f and returning the norm of the residual vector.
  double rms(vnl_vector<double> const& x);

  //: Return the number of unknowns
  unsigned int get_number_of_unknowns() const { return p_; }

  //: Return the number of residuals.
  unsigned int get_number_of_residuals() const { return n_; }

  //: Return true if the derived class has indicated that gradf has been implemented
  bool has_gradient() const { return use_gradient_; }

 protected:
  unsigned int p_;
  unsigned int n_;
  bool use_gradient_;

  void init(unsigned int number_of_unknowns, unsigned int number_of_residuals)
  { p_ = number_of_unknowns; n_ = number_of_residuals; dim_warning(p_,n_); }
 private:
  void dim_warning(unsigned int n_unknowns, unsigned int n_residuals);
};

// copy from .cpp
void
vnl_least_squares_function::dim_warning(unsigned int number_of_unknowns, unsigned int number_of_residuals)
{
    if (number_of_unknowns > number_of_residuals)
        std::cerr << "vnl_least_squares_function: WARNING: "
        << "unknowns(" << number_of_unknowns << ") > "
        << "residuals(" << number_of_residuals << ")\n";
}

void
vnl_least_squares_function::gradf(vnl_vector<double> const & /*x*/, vnl_matrix<double> & /*jacobian*/)
{
    std::cerr << "Warning: gradf() called but not implemented in derived class\n";
}

//: Compute finite differences gradient using central differences.
void
vnl_least_squares_function::fdgradf(vnl_vector<double> const & x, vnl_matrix<double> & jacobian, double stepsize)
{
    unsigned int dim = x.size();
    unsigned int n = jacobian.rows();
    assert(dim == get_number_of_unknowns());
    assert(n == get_number_of_residuals());
    assert(dim == jacobian.columns());
    
    vnl_vector<double> tx = x;
    vnl_vector<double> fplus(n);
    vnl_vector<double> fminus(n);
    for (unsigned int i = 0; i < dim; ++i)
    {
        // calculate f just to the right of x[i]
        double tplus = tx[i] = x[i] + stepsize;
        this->f(tx, fplus);
        
        // calculate f just to the left of x[i]
        double tminus = tx[i] = x[i] - stepsize;
        this->f(tx, fminus);
        
        double h = 1.0 / (tplus - tminus);
        for (unsigned int j = 0; j < n; ++j)
            jacobian(j, i) = (fplus[j] - fminus[j]) * h;
        
        // restore tx
        tx[i] = x[i];
    }
}

//: Compute finite differences gradient using forward differences.
void
vnl_least_squares_function::ffdgradf(vnl_vector<double> const & x, vnl_matrix<double> & jacobian, double stepsize)
{
    unsigned int dim = x.size();
    unsigned int n = jacobian.rows();
    assert(dim == get_number_of_unknowns());
    assert(n == get_number_of_residuals());
    assert(dim == jacobian.columns());
    
    vnl_vector<double> tx = x;
    vnl_vector<double> fplus(n);
    vnl_vector<double> fcentre(n);
    this->f(x, fcentre);
    for (unsigned int i = 0; i < dim; ++i)
    {
        // calculate f just to the right of x[i]
        double tplus = tx[i] = x[i] + stepsize;
        this->f(tx, fplus);
        
        double h = 1.0 / (tplus - x[i]);
        for (unsigned int j = 0; j < n; ++j)
            jacobian(j, i) = (fplus[j] - fcentre[j]) * h;
        
        // restore tx
        tx[i] = x[i];
    }
}

void
vnl_least_squares_function::trace(int /* iteration */,
                                  vnl_vector<double> const & /*x*/,
                                  vnl_vector<double> const & /*fx*/)
{
    // This default implementation is empty; overloaded in derived class.
}

double
vnl_least_squares_function::rms(vnl_vector<double> const & x)
{
    vnl_vector<double> fx(n_);
    f(x, fx);
    return fx.rms();
}

#endif // vnl_least_squares_function_h_
