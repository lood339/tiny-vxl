// This is core/vgl/vgl_line_2d.h
#ifndef vgl_line_2d_h_
#define vgl_line_2d_h_
//:
// \file
// \author Don Hamilton, Peter Tu, Peter Vanroose, Francois BERTEL, Franck Bettinger
// \date   2000-02-16 Don HAMILTON, Peter TU - Creation
//
// \verbatim
//  Modifications
//   2000-02-29 Peter Vanroose    Several minor fixes
//   2000-05-05 Francois BERTEL   Several minor bugs fixed
//   2000-05-09 Peter Vanroose    dist_origin() re-implemented
//   2000-12-01 Peter Vanroose    moved dist_origin() to vgl_distance.h
//   2001-03-19 Franck Bettinger  added Manchester binary IO code
//   2001-06-27 Peter Vanroose    Added operator==
//   2001-07-05 Peter Vanroose    direction, normal in terms of vgl_vector_2d
//   2001-07-06 Peter Vanroose    Added concurrent(), added assertions
//   2009-05-21 Peter Vanroose    istream operator>> re-implemented
// \endverbatim

#include <iosfwd>
#include <cmath>
#include <iostream>
#include <cassert>

#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif
#include <cassert>
#include "vgl_fwd.h" // forward declare vgl_point_2d and vgl_homg_line_2d
#include "vgl_vector_2d.h"
#include "vgl_point_2d.h"
#include "vgl_homg_line_2d.h"


//: Represents a Euclidean 2D line
// An interface for the line coefficients, [a,b,c], is provided in terms of the
// standard implicit line equation: a*x + b*y + c = 0
template <class Type>
class vgl_line_2d
{
  // the data associated with this point
  Type a_;
  Type b_;
  Type c_;

 public:
  //: Default constructor (Line 1.y==0, the X axis)
  inline vgl_line_2d() : a_(0), b_(1), c_(0) {}

  //: Construct a vgl_line_2d from its equation, three Types.
  //  The values of a and b should not be both zero.
  inline vgl_line_2d(Type ta, Type tb, Type tc) :a_(ta),b_(tb),c_(tc) { assert(ta||tb); }

  //: Construct from its equation, a 3-vector.
  //  The values v[0] and v[1] should not be both zero.
  inline vgl_line_2d(const Type v[3]):a_(v[0]),b_(v[1]),c_(v[2]) { assert(a_||b_); }

  //: Construct from homogeneous description of line
  //  The line l should not be the line at infinity.
  vgl_line_2d (vgl_homg_line_2d<Type> const& l);

  //: Construct from two distinct points (join)
  //  The two points must be distinct!
  vgl_line_2d (vgl_point_2d<Type> const& p1, vgl_point_2d<Type> const& p2);

  //: Construct from one point and one vector
  vgl_line_2d (vgl_point_2d<Type> const& p, vgl_vector_2d<Type> const& v);

  //: the comparison operator
  inline bool operator==(vgl_line_2d<Type> const& l) const
  {
    return (this==&l) ||
           (a()*l.c()==c()*l.a() && b()*l.c()==c()*l.b() && b()*l.a()==a()*l.b());
  }

  inline bool operator!=(vgl_line_2d<Type>const& other) const { return !operator==(other); }

  //: angle with the horizontal line y=0, measured in radians.
  //  Returns values between -pi and pi, i.e., the lines x-y=0 and y-x=0
  //  return different values (pi/4 and -3pi/4 respectively) although these
  //  lines are identical.
  double slope_radians() const;

  //: angle with the horizontal line y=0, measured in 360-degrees.
  //  Returns values between -180 and 180, i.e., the lines x-y=0 and y-x=0
  //  return different values (45 and -135 respectively) although these
  //  lines are identical.
  double slope_degrees() const;

  // Data Access-------------------------------------------------------------

  //: Parameter a of line a*x + b*y + c = 0
  inline Type a() const {return a_;}
  //: Parameter b of line a*x + b*y + c = 0
  inline Type b() const {return b_;}
  //: Parameter c of line a*x + b*y + c = 0
  inline Type c() const {return c_;}

  //: unit vector describing line direction
  inline vgl_vector_2d<Type> direction() const
  { return normalized(vgl_vector_2d<Type>(b_,-a_)); }

  //: unit vector orthogonal to line
  inline vgl_vector_2d<Type> normal() const
  { return normalized(vgl_vector_2d<Type>(a_,b_)); }

  //: normalize the line coefficients s.t. a^2 + b^2 = 1
  bool normalize();

  //: Set a b c.
  //  The values of a and b should not be both zero.
  //  Note that it does not make sense to set a, b or c separately
  inline void set(Type ta, Type tb, Type tc) { assert(ta||tb); a_=ta; b_=tb; c_=tc; }

  //: Return true iff this line is the line at infinity
  //  This always returns "false"
  inline bool ideal(Type = (Type)0) const { return false; }

  //: Get two points on the line; normally the intersection with X and Y axes.
  // When the line is parallel to one of these,
  // the point with \a y=1 or \a x=1, resp. are taken.  When the line goes
  // through the origin, the second point is (b, -a).
  void get_two_points(vgl_point_2d<Type> &p1, vgl_point_2d<Type> &p2) const;
};

#define l vgl_line_2d<Type>

//: Return true iff line is the line at infinity
// \relatesalso vgl_line_2d
template <class Type> inline
bool is_ideal(l const&, Type = (Type)0) { return false; }

//: Are three lines concurrent, i.e., do they pass through a common point?
// \relatesalso vgl_line_2d
template <class Type> inline
bool concurrent(l const& l1, l const& l2, l const& l3)
{
  return l1.a()*(l2.b()*l3.c()-l3.b()*l2.c())
        +l2.a()*(l3.b()*l1.c()-l1.b()*l3.c())
        +l3.a()*(l1.b()*l2.c()-l2.b()*l1.c())==0;
}

//: Write line description to stream: "<vgl_line_2d ax+by+c>"
// \relatesalso vgl_line_2d
template <class Type>
std::ostream&  operator<<(std::ostream& s, l const& line);

#undef l

// copy from .cpp
//: line through two given points
template <class Type>
vgl_line_2d<Type>::vgl_line_2d (vgl_point_2d<Type> const& p1, vgl_point_2d<Type> const& p2)
: a_ ( p1.y() - p2.y() )
, b_ ( p2.x() - p1.x() )
, c_ ( p1.x() * p2.y() - p1.y() * p2.x() )
{
    assert(a_||b_); // two points were distinct
}

//: line defined by one point and one vector
template <class Type>
vgl_line_2d<Type>::vgl_line_2d (vgl_point_2d<Type> const& p, vgl_vector_2d<Type> const& v)
: a_ ( -v.y() )
, b_ ( v.x() )
, c_ ( -a_*p.x() - b_*p.y() )
{
}

template <class Type>
vgl_line_2d<Type>::vgl_line_2d (vgl_homg_line_2d<Type> const& l)
: a_(l.a()) , b_(l.b()) , c_(l.c())
{
    //JLM I see no reason to prohibit lines through the origin
    //  assert(c_);
}

//: Get two points on the line.
// These two points are normally the intersections
// with the Y axis and X axis, respectively.  When the line is parallel to one
// of these, the point with \a y=1 or \a x=1, resp. are taken.
// When the line goes through the origin, the second point is \a (b,-a).
template <class Type>
void vgl_line_2d<Type>::get_two_points(vgl_point_2d<Type> &p1, vgl_point_2d<Type> &p2) const
{
    if (b() == 0)       p1.set(-c()/a(), 1);
    else                p1.set(0, -c()/b());
    if (a() == 0)       p2.set(1, -c()/b());
    else if ( c() == 0) p2.set(b(), -a());
    else                p2.set(-c()/a(), 0);
}

template <class Type>
double vgl_line_2d<Type>::slope_degrees() const
{
    static const double deg_per_rad = 45.0/std::atan2(1.0,1.0);
    // do special cases separately, to avoid rounding errors:
    if (a() == 0) return b()<0 ? 0.0 : 180.0;
    if (b() == 0) return a()<0 ? -90.0 : 90.0;
    if (a() == b()) return a()<0 ? -45.0 : 135.0;
    if (a()+b() == 0) return a()<0 ? -135.0 : 45.0;
    // general case:
    return deg_per_rad * std::atan2(double(a()),-double(b()));
}

template <class Type>
double vgl_line_2d<Type>::slope_radians() const
{
    return std::atan2(double(a()),-double(b()));
}

template <class Type>
bool vgl_line_2d<Type>::normalize()
{
    double mag = a_*a_ + b_*b_;
    if (mag==1.0) return true;
    if (mag==0.0) return false;
    mag = 1.0/std::sqrt(mag);
    a_ = Type(a_*mag);
    b_ = Type(b_*mag);
    c_ = Type(c_*mag);
    mag = a_*a_ + b_*b_;
    // return false when normalisation did not succeed, e.g. when Type == int:
    return mag>0.99 && mag<1.01;
}

#define vp(os,v,s) { (os)<<' '; if ((v)>0) (os)<<'+'; if ((v)&&!(s)[0]) (os)<<(v); else { \
if ((v)==-1) (os)<<'-';\
else if ((v)!=0&&(v)!=1) (os)<<(v);\
if ((v)!=0) (os)<<' '<<(s); } }

//: Write line description to stream: "<vgl_line_2d ax+by+c=0>"
template <class Type>
std::ostream&  operator<<(std::ostream& os, vgl_line_2d<Type> const& l)
{
    os << "<vgl_line_2d"; vp(os,l.a(),"x"); vp(os,l.b(),"y"); vp(os,l.c(),"");
    return os << " = 0 >";
}

#undef vp

#endif // vgl_line_2d_h_
