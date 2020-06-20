// This is core/vgl/vgl_infinite_line_3d.h
#ifndef vgl_infinite_line_3d_h_
#define vgl_infinite_line_3d_h_
//:
// \file
// \brief A 3-d infinite line with position parameterized by orthogonal plane coordinates
// \author  J.L. Mundy
//
// \verbatim
// Modifications
// Initial version July 25, 2009
// \endverbatim

#include <iosfwd>
#include <iostream>
#include <cmath>
#include <string>
#include <cassert>

#include "vgl/vgl_vector_2d.h"
#include "vgl/vgl_vector_3d.h"
#include "vgl/vgl_point_3d.h"
#include "vgl/vgl_line_segment_3d.h"
#include "vgl/vgl_line_3d_2_points.h"

//: Represents a 3-d line with position defined in the orthogonal plane passing through the origin.
//  The line direction is t_.
//  The 2-d plane coordinate system (u, v) is aligned with the 3-d coordinate
//  system (X, Y, X), where v = t x X and u = v x t.
template <class Type>
class vgl_infinite_line_3d
{
  vgl_vector_2d<Type> x0_; //!< line position vector
  vgl_vector_3d<Type> t_;  //!< line direction vector (tangent)
 public:
  //: Default constructor - does not initialise!
  inline vgl_infinite_line_3d() = default;

  //: Copy constructor
  inline vgl_infinite_line_3d(vgl_infinite_line_3d<Type> const& l)
    : x0_(l.x0()), t_(l.direction()) {}

  //: Construct from x0 and direction
  inline vgl_infinite_line_3d(vgl_vector_2d<Type> const& x_0,
                              vgl_vector_3d<Type> const& direction)
    : x0_(x_0), t_(direction) {}

  //: Construct from two points
  vgl_infinite_line_3d(vgl_point_3d<Type> const& p1,
                       vgl_point_3d<Type> const& p2);

  //: Construct from a point and direction
  vgl_infinite_line_3d(vgl_point_3d<Type> const& p,
                       vgl_vector_3d<Type> const& direction);

  //: Construct from a line segment
  inline vgl_infinite_line_3d(vgl_line_segment_3d<Type> const& ls)
  {
    vgl_infinite_line_3d<Type> inf_l(ls.point1(), ls.point2());
    x0_ = inf_l.x0(); t_ = inf_l.direction();
  }

  //: Construct from a line 2 points
  inline vgl_infinite_line_3d(vgl_line_3d_2_points<Type> const& ls)
  {
    vgl_infinite_line_3d<Type> inf_l(ls.point1(), ls.point2());
    x0_ = inf_l.x0(); t_ = inf_l.direction();
  }

  //: Destructor
  inline ~vgl_infinite_line_3d() = default;

  //: Accessors
  inline vgl_vector_2d<Type> x0() const { return x0_; } // return a copy
  inline vgl_vector_3d<Type> direction() const
  { return t_/static_cast<Type>(t_.length()); } // return a copy

  //: The comparison operator
  inline bool operator==(vgl_infinite_line_3d<Type> const& l) const
  { return (this==&l) || (direction() == l.direction() && x0() == l.x0()); }

  inline bool operator!=(vgl_infinite_line_3d<Type>const& other) const
  { return !operator==(other); }

  //: Assignment
  inline void set(vgl_vector_2d<Type> const& x_0, vgl_vector_3d<Type> const& direction)
  { x0_ = x_0; t_ = direction; }

  //: Return the point on the line closest to the origin
  vgl_point_3d<Type> point() const;

  //: Return a point on the line defined by a scalar parameter \a t.
  // \a t=0.0 corresponds to the closest point on the line to the origin
  vgl_point_3d<Type> point_t(const double t) const { return point() + t*direction(); }

  //: Check if point \a p is on the line
  bool contains(const vgl_point_3d<Type>& p ) const;

  //: The unit vectors perpendicular to the line direction
  void compute_uv_vectors(vgl_vector_3d<Type>& u, vgl_vector_3d<Type>& v) const;
};

//: Write to stream
// \relatesalso vgl_infinite_line_3d
template <class Type>
std::ostream&  operator<<(std::ostream& s, const vgl_infinite_line_3d<Type>& p);

// copy from .cpp
template <class Type>
vgl_infinite_line_3d<Type>::vgl_infinite_line_3d(vgl_point_3d<Type> const& p1,
                                                 vgl_point_3d<Type> const& p2)
{
    vgl_vector_3d<Type> dir = p2-p1;
    vgl_infinite_line_3d<Type> l(p1, dir);
    x0_ = l.x0();
    t_ = dir;
}

template <class Type>
void vgl_infinite_line_3d<Type>::
compute_uv_vectors(vgl_vector_3d<Type>& u, vgl_vector_3d<Type>& v) const
{
    // define the plane coordinate system (u, v)
    // v is given by the cross product of t with x, unless t is nearly
    // parallel to x, in which case v is given by z X t.
    vgl_vector_3d<Type> x(Type(1), Type(0), Type(0));
    v = cross_product(t_,x);
    Type vmag = static_cast<Type>(v.length());
    double vmagd = static_cast<double>(vmag);
    if (vmagd < 1.0e-8) {
        vgl_vector_3d<Type> z(Type(0), Type(0), Type(1));
        v = cross_product(z, t_);
        vmag = static_cast<Type>(v.length());
        assert(vmag>Type(0));
        v/=vmag;
    }
    else v/=vmag;
    // The other plane coordinate vector is perpendicular to both t and v
    u = cross_product(v,t_);
    Type umag = static_cast<Type>(u.length());
    u/=umag;
}

template <class Type>
vgl_infinite_line_3d<Type>::
vgl_infinite_line_3d(vgl_point_3d<Type> const& p,
                     vgl_vector_3d<Type> const& dir)
{
    // reconcile direction so that tangent is in the positive hemisphere
    double ttx = std::fabs(static_cast<double>(dir.x()));
    double tty = std::fabs(static_cast<double>(dir.y()));
    double ttz = std::fabs(static_cast<double>(dir.z()));
    double max_comp = ttx;
    double sign = static_cast<double>(dir.x());
    if (max_comp < tty) {
        max_comp = tty;
        sign = static_cast<double>(dir.y());
    }
    if (max_comp < ttz) {
        max_comp = ttz;
        sign = static_cast<double>(dir.z());
    }
    // switch sense if max component is negative
    Type sense = static_cast<Type>(sign/max_comp);
    t_ = normalized(dir*sense);
    // Define the plane perpendicular to the line passing through the origin
    // the plane normal is t_ the distance of the plane from the origin is 0
    // it follows that the intersection of the line with the perpendicular plane
    // is as follows:
    Type mag = static_cast<Type>(t_.length());
    assert(mag>Type(0));
    vgl_vector_3d<Type> pv(p.x(), p.y(), p.z());
    Type dp = dot_product(pv, t_);
    Type k = -dp/(mag*mag);
    // The intersection point
    vgl_vector_3d<Type> p0 = pv + k*t_, u, v;
    this->compute_uv_vectors(u, v);
    // The location of the intersection point in plane coordinates can now be computed
    Type u0 = dot_product(u, p0), v0 = dot_product(v, p0);
    x0_.set(u0, v0);
}

// the point on the line closest to the origin
template <class Type>
vgl_point_3d<Type> vgl_infinite_line_3d<Type>::point() const
{
    // u,v plane coordinate vectors
    vgl_vector_3d<Type> u, v, pv;
    this->compute_uv_vectors(u, v);
    pv = x0_.x()*u + x0_.y()*v;
    return vgl_point_3d<Type>(pv.x(), pv.y(), pv.z());
}

template <class Type>
bool vgl_infinite_line_3d<Type>::contains(const vgl_point_3d<Type>& p ) const
{
    vgl_point_3d<Type> point1 = this->point();
    vgl_point_3d<Type> point2 = this->point_t(Type(1));
    double seg = 1.0;
    double len1 = static_cast<double>((point1 - p).length());
    double len2 = static_cast<double>((point2 - p).length());
    // two cases: point inside (point1, point2) segment;
    //            point outside (point1, point2) segment
    double r = seg -(len1 + len2);
    if (len1>seg||len2>seg)
        r = seg - std::fabs(len1-len2);
    return r < 1e-8 && r > -1e-8;
}

// stream operators
template <class Type>
std::ostream& operator<<(std::ostream& s, vgl_infinite_line_3d<Type> const & p)
{
    return s << "<vgl_infinite_line_3d: origin " << p.x0() << " dir " << p.direction() << " >";
}

#endif // vgl_infinite_line_3d_h_
