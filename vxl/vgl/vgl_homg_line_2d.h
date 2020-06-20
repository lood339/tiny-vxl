// This is core/vgl/vgl_homg_line_2d.h
#ifndef vgl_homg_line_2d_h
#define  vgl_homg_line_2d_h
//:
// \file
// \brief line in projective 2D space
// \author Don Hamilton, Peter Tu
//
// \verbatim
//  Modifications
//   Peter Vanroose -  6 July 2001 - Added normal(), direction() and concurrent()
//   Peter Vanroose -  4 July 2001 - Added assertions and cstr from non-homg line
//   Peter Vanroose - 27 June 2001 - Added operator==
// \endverbatim

#include <iosfwd>
#include <cassert>

#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif
#include "vgl_fwd.h" // forward declare vgl_homg_point_2d and vgl_line_2d

#include "vgl_vector_2d.h"
#include "vgl_homg_point_2d.h"
#include "vgl_line_2d.h"

//: Represents a homogeneous 2D line.
template <class T>
class vgl_homg_line_2d
{
  //: the data associated with this line
  T a_;
  T b_;
  T c_;

 public:

  // Constructors/Initializers/Destructor------------------------------------

  //: Default constructor (Line 1.y==0, the X axis)
  inline vgl_homg_line_2d() : a_(0), b_(1), c_(0) {}

  //: Construct from three Types.
  //  The three given numbers should not be all 0
  inline vgl_homg_line_2d(T va, T vb, T vc) : a_(va), b_(vb), c_(vc) {assert(va||vb||vc);}

  //: Construct from 3-vector.
  //  The three given numbers should not be all 0
  inline vgl_homg_line_2d(const T v[3]) : a_(v[0]), b_(v[1]), c_(v[2]) {assert(a_||b_||c_);}

  //: Construct from non-homogeneous line
  vgl_homg_line_2d<T> (vgl_line_2d<T> const& p);

  //: Construct from two distinct points (join)
  //  The two points must be distinct!
  vgl_homg_line_2d(vgl_homg_point_2d<T> const& p1, vgl_homg_point_2d<T> const& p2);

  //: the comparison operator
  inline bool operator==(vgl_homg_line_2d<T> const& l) const
  {
    return (this==&l) ||
           (a()*l.c()==c()*l.a() && b()*l.c()==c()*l.b() && b()*l.a()==a()*l.b());
  }

  inline bool operator!=(vgl_homg_line_2d<T> const& other)const{return !operator==(other);}

  // Data Access-------------------------------------------------------------

  //: Parameter a of line a*x + b*y + c*w = 0
  inline T a() const {return a_;}
  //: Parameter b of line a*x + b*y + c*w = 0
  inline T b() const {return b_;}
  //: Parameter c of line a*x + b*y + c*w = 0
  inline T c() const {return c_;}

  //: unit vector describing line direction, or (0,0) if line at infinity
  inline vgl_vector_2d<double> direction() const { return normalized(vgl_vector_2d<double>(b_,-a_)); }

  //: unit vector orthogonal to line, or (0,0) if line at infinity
  inline vgl_vector_2d<double> normal() const { return normalized(vgl_vector_2d<double>(a_,b_)); }

  //: divide all coefficients by sqrt(a^2 + b^2)
  void normalize();

  //: Set a b c.
  //  The three given numbers should not be all 0
  //  Note that it does not make sense to set a, b or c separately
  inline void set(T va, T vb, T vc) {assert(va||vb||vc); a_=va; b_=vb; c_=vc;}

  //: Return true iff this line is the line at infinity
  //  This version checks (max(|a|,|b|) <= tol * |c|
  inline bool ideal(T tol = (T)0) const
  {
#define vgl_Abs(x) ((x)<0?-(x):(x)) // avoid #include of vcl_cmath.h AND vcl_cstdlib.h
    return vgl_Abs(a()) <= tol*vgl_Abs(c()) && vgl_Abs(b()) <= tol*vgl_Abs(c());
#undef vgl_Abs
  }

  //:get two points on the line
  // These two points are normally the intersections
  // with the Y axis and X axis, respectively.  When the line is parallel to one
  // of these, the point with y=1 or x=1, resp. are taken.  When the line goes
  // through the origin, the second point is (b, -a, 1).  Finally, when the line
  // is the line at infinity, the returned points are (1,0,0) and (0,1,0).
  // Thus, whenever possible, the returned points are not at infinity.
  void get_two_points(vgl_homg_point_2d<T> &p1, vgl_homg_point_2d<T> &p2) const;
};

#define l vgl_homg_line_2d<T>

//: Return true iff line is the line at infinity
//  This version checks (max(|a|,|b|) <= tol * |c|
// \relatesalso vgl_homg_line_2d
template <class T>
inline bool is_ideal(l const& line, T tol = (T)0) { return line.ideal(tol); }

//: Are three lines concurrent, i.e., do they pass through a common point?
// \relatesalso vgl_homg_line_2d
template <class T>
inline bool concurrent(l const& l1, l const& l2, l const& l3)
{
  return l1.a()*(l2.b()*l3.c()-l3.b()*l2.c())
        +l2.a()*(l3.b()*l1.c()-l1.b()*l3.c())
        +l3.a()*(l1.b()*l2.c()-l2.b()*l1.c())==0;
}

//: Print line equation to stream
// \relatesalso vgl_homg_line_2d
template <class T>
std::ostream& operator<<(std::ostream& s, l const& line);
#undef l

// copy from .cpp file
template <class Type>
vgl_homg_line_2d<Type>::vgl_homg_line_2d (vgl_line_2d<Type> const& l)
: a_(l.a()) , b_(l.b()) , c_(l.c())
{
}

//: get two points on the line.
//  These two points are normally the intersections with the Y axis and X axis,
//  respectively.  When the line is parallel to one of these, the point with
//  \a y/w=1 or \a x/w=1, resp. are taken.  When the line goes through the origin,
//  the second point is (b, -a, 1).  Finally, when the line is the line at
//  infinity, the returned points are (1,0,0) and (0,1,0).
//
//  Thus, whenever possible, the returned points are not at infinity.
//
template <class Type>
void vgl_homg_line_2d<Type>::get_two_points(vgl_homg_point_2d<Type> &p1, vgl_homg_point_2d<Type> &p2) const
{
    if (      b() == 0) p1.set(-c(), a(), a());
    else                p1.set(0, -c(), b());
    if (      a() == 0) p2.set(b(), -c(), b());
    else if ( c() == 0) p2.set(b(), -a(), 1);
    else                p2.set(-c(), 0, a());
}

template <class Type>
vgl_homg_line_2d<Type>::vgl_homg_line_2d (vgl_homg_point_2d<Type> const& p1,
                                          vgl_homg_point_2d<Type> const& p2)
{
    set(p1.y()*p2.w()-p1.w()*p2.y(),
        p1.w()*p2.x()-p1.x()*p2.w(),
        p1.x()*p2.y()-p1.y()*p2.x());
    assert(a_||b_||c_); // given points should be different
}

#define vp(os,v,s) { (os)<<' '; if ((v)>0) (os)<<'+';\
if ((v)&&!(s)[0]) (os)<<(v); else { \
if ((v)==-1) (os)<<'-';\
else if ((v)!=0&&(v)!=1) (os)<<(v);\
if ((v)!=0) (os)<<' '<<(s); } }

//: Print line equation to stream
template <class Type>
std::ostream&  operator<<(std::ostream& os, vgl_homg_line_2d<Type>const& l)
{
    os << "<vgl_homg_line_2d"; vp(os,l.a(),"x"); vp(os,l.b(),"y"); vp(os,l.c(),"w");
    return os << " = 0 >";
}

#undef vp

//: Load in line parameters from stream
template <class Type>
std::istream&  operator>>(std::istream& is, vgl_homg_line_2d<Type>& p)
{
    Type a,b,c;
    is >> a >> b >> c;
    p.set(a,b,c);
    return is;
}

template <class Type>
void vgl_homg_line_2d<Type>::normalize()
{
    double sum = a_*a_ + b_*b_;
    double den = std::sqrt(sum);
    if (den<1.0e-8)//don't normalize ideal line
        return;
    double an= (double)a()/den;
    double bn= (double)b()/den;
    double cn= (double)c()/den;
    //standardize so that a is positive unless a is smaller than b, then
    //standardize the sign of b
    if (std::fabs(an)>std::fabs(bn))
        if (an>0)
        {
            a_ = (Type)an;
            b_ = (Type)bn;
            c_ = (Type)cn;
        }
        else
        {
            a_ = -(Type)an;
            b_ = -(Type)bn;
            c_ = -(Type)cn;
        }
        else
            if (bn>0)
            {
                a_ = (Type)an;
                b_ = (Type)bn;
                c_ = (Type)cn;
            }
            else
            {
                a_ = -(Type)an;
                b_ = -(Type)bn;
                c_ = -(Type)cn;
            }
    return;
}





#endif //  vgl_homg_line_2d_h
