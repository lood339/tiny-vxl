// This is core/vgl/vgl_homg_plane_3d.h
#ifndef vgl_homg_plane_3d_h
#define vgl_homg_plane_3d_h
//:
// \file
// \brief homogeneous plane in 3D projective space
// \author Don HAMILTON Peter TU
//
// \verbatim
//  Modifications
//   Peter Vanroose  6 July 2001: Now using vgl_vector_3d for normal direction
//   Peter Vanroose  6 July 2001: Added normal(); replaced data_[4] by a_ b_ c_ d_
//   Peter Vanroose  6 July 2001: Added constructor from 3 points
//   CJB (Manchester) 16/03/2001: Tidied up the documentation
//   Peter Vanroose 15 July 2002: Added constructor from two concurrent lines
//   cbw (imorphics) 31 Mar 2008: Corrected constructor with normal and point (negated d) and added test
//   Peter Vanroose 7 March 2009: Added normalize(), similar to the one in vgl_homg_line_2d<T>
// \endverbatim

#include <iosfwd>
#include <cassert>
#include <cmath>
#include <iostream>

#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif
#include "vgl_fwd.h" // forward declare vgl_plane_3d and vgl_homg_point_3d
#include "vgl_vector_3d.h"
#include "vgl_plane_3d.h"
#include "vgl_homg_point_3d.h"
//#include "vgl_homg_line_3d_2_points.h"


//: Represents a homogeneous 3D plane
template <class Type>
class vgl_homg_plane_3d
{
  // the four homogeneous coordinates of the plane (dual of a point).
  Type a_;
  Type b_;
  Type c_;
  Type d_;

 public:
  inline vgl_homg_plane_3d () = default;

  //: Construct from four Types.
  inline vgl_homg_plane_3d(Type ta, Type tb, Type tc, Type td) : a_(ta), b_(tb), c_(tc), d_(td) {}

  //: Construct from 4-vector.
  inline vgl_homg_plane_3d(const Type v[4]) : a_(v[0]), b_(v[1]), c_(v[2]), d_(v[3]) {}

  //: Construct from non-homogeneous plane.
  vgl_homg_plane_3d (vgl_plane_3d<Type> const& pl);

  //: Construct from normal and a point
  //  The given point must not be at infinity.
  vgl_homg_plane_3d (vgl_vector_3d<Type> const& n, vgl_homg_point_3d<Type> const& p);

  //: Construct from three non-collinear points
  vgl_homg_plane_3d (vgl_homg_point_3d<Type> const& p1,
                     vgl_homg_point_3d<Type> const& p2,
                     vgl_homg_point_3d<Type> const& p3);
    
    /*
  //: Construct from two concurrent lines
  vgl_homg_plane_3d (vgl_homg_line_3d_2_points<Type> const& l1,
                     vgl_homg_line_3d_2_points<Type> const& l2);
     */

  // Data Access-------------------------------------------------------------

  //: Return \a x coefficient
  inline Type a() const {return a_;}
  inline Type nx() const {return a_;}
  //: Return \a y coefficient
  inline Type b() const {return b_;}
  inline Type ny() const {return b_;}
  //: Return \a z coefficient
  inline Type c() const {return c_;}
  inline Type nz() const {return c_;}
  //: Return homogeneous scaling coefficient
  inline Type d() const {return d_;}

  //: Set equation \a a*x+b*y+c*z+d*w=0
  inline void set(Type ta,Type tb,Type tc,Type td) { assert(ta||tb||tc||td);a_=ta;b_=tb;c_=tc;d_=td; }

  //: the comparison operator
  bool operator==( vgl_homg_plane_3d<Type> const& pl) const;
  inline bool operator!=( vgl_homg_plane_3d<Type>const& pl) const { return !operator==(pl); }

  //: Return true iff the plane is the plane at infinity.
  // The method checks that max(|a|,|b|,|c|) <= tol * |d|
  // If called without an argument, tol=0, i.e., a, b and c must be 0.
  inline bool ideal(Type tol = (Type)0) const
  {
#define vgl_Abs(x) ((x)<0?-(x):(x)) // avoid #include of vcl_cmath.h AND vcl_cstdlib.h
    return vgl_Abs(a()) <= tol * vgl_Abs(d()) &&
           vgl_Abs(b()) <= tol * vgl_Abs(d()) &&
           vgl_Abs(c()) <= tol * vgl_Abs(d());
#undef vgl_Abs
  }

  inline vgl_vector_3d<double> normal() const { return normalized(vgl_vector_3d<double>(a(),b(),c())); }

  //: divide all coefficients by sqrt(a^2 + b^2 + c^2)
  void normalize();
};

//: Return true iff p is the plane at infinity
// The method checks that max(|a|,|b|,|c|) <= tol * |d|
// \relatesalso vgl_homg_plane_3d
template <class Type>
inline bool is_ideal(vgl_homg_plane_3d<Type> const& p, Type tol=(Type)0) { return p.ideal(tol); }

// stream operators

template <class Type>
std::ostream&  operator<<(std::ostream& s, const vgl_homg_plane_3d<Type>& p);

// copy from .cpp
//: Construct from non-homogeneous plane
template <class Type>
vgl_homg_plane_3d<Type>::vgl_homg_plane_3d(vgl_plane_3d<Type> const& pl)
: a_(pl.a()), b_(pl.b()), c_(pl.c()), d_(pl.d()) {}

//: Construct from three points
template <class Type>
vgl_homg_plane_3d<Type>::vgl_homg_plane_3d (vgl_homg_point_3d<Type> const& p1,
                                            vgl_homg_point_3d<Type> const& p2,
                                            vgl_homg_point_3d<Type> const& p3)
: a_(p1.w()*(p2.y()*p3.z()-p2.z()*p3.y())
     +p2.w()*(p3.y()*p1.z()-p3.z()*p1.y())
     +p3.w()*(p1.y()*p2.z()-p1.z()*p2.y()))
, b_(p1.w()*(p2.z()*p3.x()-p2.x()*p3.z())
     +p2.w()*(p3.z()*p1.x()-p3.x()*p1.z())
     +p3.w()*(p1.z()*p2.x()-p1.x()*p2.z()))
, c_(p1.w()*(p2.x()*p3.y()-p2.y()*p3.x())
     +p2.w()*(p3.x()*p1.y()-p3.y()*p1.x())
     +p3.w()*(p1.x()*p2.y()-p1.y()*p2.x()))
, d_(p1.x()*(p2.z()*p3.y()-p2.y()*p3.z())
     +p2.x()*(p1.y()*p3.z()-p1.z()*p3.y())
     +p3.x()*(p1.z()*p2.y()-p1.y()*p2.z()))
{
    assert(a_||b_||c_||d_); // points should not be collinear or coinciding
}

/*
//: Construct from two concurrent lines
template <class Type>
vgl_homg_plane_3d<Type>::vgl_homg_plane_3d(vgl_homg_line_3d_2_points<Type> const& l1,
                                           vgl_homg_line_3d_2_points<Type> const& l2)
{
    assert(concurrent(l1,l2));
    vgl_homg_point_3d<Type> p1 = l1.point_finite();
    vgl_homg_point_3d<Type> p2 = l1.point_infinite();
    vgl_homg_point_3d<Type> p3 = l2.point_finite();
    if (collinear(p1,p2,p3)) p3 = l2.point_infinite();
    *this = vgl_homg_plane_3d<Type>(p1,p2,p3);
}
 */

//: Construct from normal and a point
// This will fail when the given point is at infinity!
template <class Type>
vgl_homg_plane_3d<Type>::vgl_homg_plane_3d(vgl_vector_3d<Type> const& n,
                                           vgl_homg_point_3d<Type> const&p)
: a_(n.x()*p.w()), b_(n.y()*p.w()), c_(n.z()*p.w()),
d_(-(n.x()*p.x()+n.y()*p.y()+n.z()*p.z())) {}

//: divide all coefficients by sqrt(a^2 + b^2 + c^2)
template <class Type>
void vgl_homg_plane_3d<Type>::normalize()
{
    double sum = a_*a_ + b_*b_ + c_*c_;
    if (sum<1e-12) // don't normalize plane at infinity
        return;
    double den = std::sqrt(sum);
    double an= (double)a()/den; a_ = (Type)an;
    double bn= (double)b()/den; b_ = (Type)bn;
    double cn= (double)c()/den; c_ = (Type)cn;
    double dn= (double)d()/den; d_ = (Type)dn;
    //standardize so that largest of (|a|,|b|,|c|) is positive
    if ((std::fabs(an)>=std::fabs(bn) && std::fabs(an)>=std::fabs(cn) && an < 0) ||
        (std::fabs(bn)>std::fabs(an) && std::fabs(bn)>=std::fabs(cn) && bn < 0) ||
        (std::fabs(cn)>std::fabs(an) && std::fabs(cn)>std::fabs(bn) && cn < 0))
        a_ = -a_, b_ = -b_, c_ = -c_, d_ = -d_;
    return;
}

template <class Type>
bool vgl_homg_plane_3d<Type>::operator==(vgl_homg_plane_3d<Type> const& p) const
{
    return (this==&p) ||
    (   (a()*p.b()==p.a()*b())
     && (a()*p.c()==p.a()*c())
     && (a()*p.d()==p.a()*d())
     && (b()*p.c()==p.b()*c())
     && (b()*p.d()==p.b()*d())
     && (c()*p.d()==p.c()*d()) );
}

template <class Type>
std::ostream& operator<<(std::ostream& s, const vgl_homg_plane_3d<Type>& p)
{
    return s << " <vgl_homg_plane_3d "
    << p.a() << " x + "
    << p.b() << " y + "
    << p.c() << " z + "
    << p.d() << " w = 0 >";
}

#endif // vgl_homg_plane_3d_h
