// This is core/vgl/vgl_plane_3d.h
#ifndef vgl_plane_3d_h
#define vgl_plane_3d_h
//:
// \file
// \brief a plane in 3D nonhomogeneous space
// \author Don Hamilton, Peter Tu
// \date   Feb 15 2000
//
// \verbatim
//  Modifications
//   Peter Vanroose  6 July 2001: Added assertion in constructors
//   Peter Vanroose  6 July 2001: Now using vgl_vector_3d for normal direction
//   Peter Vanroose  6 July 2001: Implemented constructor from 3 points
//   Peter Vanroose  6 July 2001: Added normal(); replaced data_[4] by a_ b_ c_ d_
//   Peter Vanroose  6 July 2001: Added operator== and operator!=
//   Peter Vanroose 19 Aug. 2004: implementation of both constructors corrected
//   Peter Vanroose 21 May  2009: istream operator>> re-implemented
// \endverbatim

#include <iosfwd>
#include <cassert>
#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif
#include <cassert>
#include "vgl_fwd.h" // forward declare vgl_homg_plane_3d, vgl_point_3d
#include "vgl_vector_3d.h"
#include "vgl_homg_plane_3d.h"
#include "vgl_point_3d.h"
#include "vgl_point_2d.h"
#include "vgl_ray_3d.h"
#include "vgl_closest_point.h"
#include "vgl_distance.h"
#include "vgl_tolerance.h"



//: Represents a Euclidean 3D plane
//  The equation of the plane is $ a x + b y + c z + d = 0 $
template <class T>
class vgl_plane_3d
{
  // the data associated with this plane
  T a_;
  T b_;
  T c_;
  T d_;

 public:

  // Constructors/Initializers/Destructor------------------------------------

  // Default constructor: horizontal XY-plane (equation 1.z = 0)
  inline vgl_plane_3d () : a_(0), b_(0), c_(1), d_(0) {}
    
  //: Construct a vgl_plane_3d from its equation $ax+by+cz+d=0$
  //  At least one of a, b or c should be nonzero.
  inline vgl_plane_3d (T ta,T tb,T tc,T td)
    : a_(ta), b_(tb), c_(tc), d_(td) { assert(ta||tb||tc); }

  //: Construct a vgl_plane_3d from its equation $v[0]x+v[1]y+v[2]z+v[3]=0$
  //  At least one of v[0], v[1] or v[2] should be nonzero.
  inline vgl_plane_3d (const T v[4])
    : a_(v[0]), b_(v[1]), c_(v[2]), d_(v[3]) { assert(a_||b_||c_); }

  //: Construct from a homogeneous plane
  vgl_plane_3d (vgl_homg_plane_3d<T> const& p);

  //: Construct from Normal and a point
  //  The plane goes through the point \a p and will be orthogonal to \a normal.
  vgl_plane_3d (vgl_vector_3d<T> const& normal,
                vgl_point_3d<T> const& p);

  //: Construct from three non-collinear points
  //  The plane will contain all three points \a p1, \a p2 and \a p3.
  vgl_plane_3d (vgl_point_3d<T> const& p1,
                vgl_point_3d<T> const& p2,
                vgl_point_3d<T> const& p3);

  //: Construct from two non-skew rays. The rays intersect at their origins
  // or are parallel. The plane will contain the two rays
  vgl_plane_3d (vgl_ray_3d<T> const& r0, vgl_ray_3d<T> const& r1);

  // Data Access-------------------------------------------------------------

  //: Return \a x coefficient
  inline T a()  const {return a_;}
  inline T nx() const {return a_;}
  //: Return \a y coefficient
  inline T b()  const {return b_;}
  inline T ny() const {return b_;}
  //: Return \a z coefficient
  inline T c()  const {return c_;}
  inline T nz() const {return c_;}
  //: Return constant coefficient
  inline T d()  const {return d_;}

  //: Set this vgl_plane_3d to have the equation $ax+by+cz+d=0$
  inline void set(T ta,T tb,T tc,T td) { assert(ta||tb||tc); a_=ta; b_=tb; c_=tc; d_=td; }

  //: the comparison operator
  //  The equations need not be identical, but just equivalent.
  bool operator==( vgl_plane_3d<T> const& p) const;
  inline bool operator!=( vgl_plane_3d<T>const& p) const { return !operator==(p); }

  //: Return true iff the plane is the plane at infinity.
  //  Always returns false
  inline bool ideal(T = (T)0) const { return false; }

  // divide all plane coefs by sqrt(a^2 +b^2 +c^2)
  bool normalize();

  //: Return the normal direction, i.e., a unit vector orthogonal to this plane
  inline vgl_vector_3d<T> normal() const
  { return normalized(vgl_vector_3d<T>(a(),b(),c())); }

  //: Return true if p is on the plane
  bool contains(vgl_point_3d<T> const& p, T tol = (T)0) const;

  // the coordinate system on the plane, (u, v), is defined as,
  // cases:
  // 1) n not parallel to Y
  //   u = Y x n and v = n x u
  //
  // 2) n parallel to Y
  //   u = n x Z and v = u x n
  //
  // the plane origin is the point in the plane closest to the world origin

  //: Given a 3-d point, return a 2-d point in the coord. system of the plane
  // If the point is not on the plane then false is returned
  bool plane_coords(vgl_point_3d<T> const& p3d,
                    vgl_point_2d<T>& p2d, T tol=(T)0 ) const;

  //: inverse map from plane coordinates to world coordinates
  vgl_point_3d<T> world_coords(vgl_point_2d<T> const& p2d) const;

  //: plane coordinate unit vectors
  void plane_coord_vectors(vgl_vector_3d<T>& uvec,
                           vgl_vector_3d<T>& vvec) const;
};

//: Return true iff p is the plane at infinity
//  Always returns false
template <class T> inline
bool is_ideal(vgl_plane_3d<T> const&, T tol=(T)0) { return false; }


//: Write to stream
// \relatesalso vgl_plane_3d
template <class T>
std::ostream& operator<<(std::ostream& s, const vgl_plane_3d<T>& p);

// copy from .cpp
//: Construct from homogeneous plane
template <class T>
vgl_plane_3d<T>::vgl_plane_3d(vgl_homg_plane_3d<T> const& p)
: a_(p.a()), b_(p.b()), c_(p.c()), d_(p.d()) {assert (a_||b_||c_);}

//: Construct from three points
template <class T>
vgl_plane_3d<T>::vgl_plane_3d (vgl_point_3d<T> const& p1,
                               vgl_point_3d<T> const& p2,
                               vgl_point_3d<T> const& p3)
: a_(p2.y()*p3.z()-p2.z()*p3.y()
     +p3.y()*p1.z()-p3.z()*p1.y()
     +p1.y()*p2.z()-p1.z()*p2.y())
, b_(p2.z()*p3.x()-p2.x()*p3.z()
     +p3.z()*p1.x()-p3.x()*p1.z()
     +p1.z()*p2.x()-p1.x()*p2.z())
, c_(p2.x()*p3.y()-p2.y()*p3.x()
     +p3.x()*p1.y()-p3.y()*p1.x()
     +p1.x()*p2.y()-p1.y()*p2.x())
, d_(p1.x()*(p2.z()*p3.y()-p2.y()*p3.z())
     +p2.x()*(p3.z()*p1.y()-p3.y()*p1.z())
     +p3.x()*(p1.z()*p2.y()-p1.y()*p2.z()))
{
    assert(a_||b_||c_); // points should not be collinear or coinciding
}

//: Construct from normal and a point
template <class T>
vgl_plane_3d<T>::vgl_plane_3d(vgl_vector_3d<T> const& n,
                              vgl_point_3d<T> const& p)
: a_(n.x()), b_(n.y()), c_(n.z()), d_(-(n.x()*p.x()+n.y()*p.y()+n.z()*p.z()))
{
    assert(a_||b_||c_); // normal vector should not be the null vector
}


template <class T>
vgl_plane_3d<T>::vgl_plane_3d (vgl_ray_3d<T> const& r0,
                               vgl_ray_3d<T> const& r1){
    // check if the rays are parallel
    const vgl_vector_3d<T>& v0 = r0.direction();
    const vgl_vector_3d<T>& v1 = r1.direction();
    double  para = std::fabs(1.0-std::fabs(cos_angle(v0, v1)));
    bool parallel = para < vgl_tolerance<double>::position;
    // check if the ray origins are coincident
    const vgl_point_3d<T>& p0 = r0.origin();
    const vgl_point_3d<T>& p1 = r1.origin();
    double d01 = length(p1-p0);
    bool coincident = d01 < vgl_tolerance<double>::position;
    
    // assert the rays are distinct
    bool distinct = !parallel || (parallel&&!coincident);
    assert(distinct);
    // assert the rays are not skew
    bool not_skew = (parallel&&distinct) || (!parallel&&coincident);
    assert(not_skew);
    
    // Case I: coincident
    if(coincident){
        vgl_vector_3d<T> norm = cross_product(v0, v1);
        vgl_plane_3d<T> pln(norm, p0);
        a_=pln.a();   b_=pln.b();   c_=pln.c();   d_=pln.d();
        return;
    }
    // Case II: parallel
    vgl_vector_3d<T> v01 = p1-p0;
    vgl_vector_3d<T> norm = cross_product(v0, v01);
    vgl_plane_3d<T> pln(norm, p0);
    a_=pln.a();   b_=pln.b();   c_=pln.c();   d_=pln.d();
}

template <class T>
bool vgl_plane_3d<T>::normalize(){
    double sum = a_*a_ + b_*b_ + c_*c_;
    if (sum<1e-12) // don't normalize plane at infinity
        return false;
    double den = std::sqrt(sum);
    double an= (double)a()/den; a_ = (T)an;
    double bn= (double)b()/den; b_ = (T)bn;
    double cn= (double)c()/den; c_ = (T)cn;
    double dn= (double)d()/den; d_ = (T)dn;
    //standardize so that largest of (|a|,|b|,|c|) is positive
    if ((std::fabs(an)>=std::fabs(bn) && std::fabs(an)>=std::fabs(cn) && an < 0) ||
        (std::fabs(bn)>std::fabs(an) && std::fabs(bn)>=std::fabs(cn) && bn < 0) ||
        (std::fabs(cn)>std::fabs(an) && std::fabs(cn)>std::fabs(bn) && cn < 0))
        a_ = -a_, b_ = -b_, c_ = -c_, d_ = -d_;
    return true;
}
//: Return true if p is on the plane
template <class T>
bool vgl_plane_3d<T>::contains(vgl_point_3d<T> const& p, T tol) const
{
    //to maintain a consistent distance metric the plane should be normalized
    vgl_vector_3d<T> n(a_, b_, c_), pv(p.x(), p.y(), p.z());
    T dist = (dot_product(n,pv) + d_) / static_cast<T>(length(n));
    return dist >= -tol && dist <= tol;
}

template <class T>
bool vgl_plane_3d<T>::operator==(vgl_plane_3d<T> const& p) const
{
    return (this==&p) ||
    (   (a()*p.b()==p.a()*b())
     && (a()*p.c()==p.a()*c())
     && (a()*p.d()==p.a()*d())
     && (b()*p.c()==p.b()*c())
     && (b()*p.d()==p.b()*d())
     && (c()*p.d()==p.c()*d()) );
}

#define vp(os,v,s)  os<<' '<< v << ' ' <<s;

template <class T>
std::ostream& operator<<(std::ostream& os, const vgl_plane_3d<T>& p)
{
    os << "<vgl_plane_3d"; vp(os,p.a(),"x"); vp(os,p.b(),"y"); vp(os,p.c(),"z");
    vp(os,p.d(),""); return os << " = 0 >";
}

#undef vp

template <class T>
std::istream& operator>>(std::istream& is, vgl_plane_3d<T>& p)
{
    if (! is.good()) return is; // (TODO: should throw an exception)
    bool paren = false;
    bool formatted = false;
    is >> std::ws; // jump over any leading whitespace
    char ch;
    ch = is.peek();
    if(ch == '<'){
        std::string temp;
        is >> temp;
        formatted = true;
    }
    is >> std::ws;
    T a, b, c, d;
    if (is.eof()) return is; // nothing to be set because of EOF (TODO: should throw an exception)
    if (is.peek() == '(') { is.ignore(); paren=true; }
    is >> std::ws >> a >> std::ws;
    if (is.eof()) return is;
    if (is.peek() == ',') is.ignore();
    else if (is.peek() == 'x') { is.ignore(); formatted=true; }
    is >> std::ws >> b >> std::ws;
    if (is.eof()) return is;
    if (formatted) {
        if (is.eof()) return is;
        if (is.peek() == 'y') is.ignore();
        else                  return is; // formatted input incorrect (TODO: throw an exception)
    }
    else if (is.peek() == ',') is.ignore();
    is >> std::ws >> c >> std::ws;
    if (is.eof()) return is;
    if (formatted) {
        if (is.eof()) return is;
        if (is.peek() == 'z') is.ignore();
        else                  return is; // formatted input incorrect (TODO: throw an exception)
    }
    else if (is.peek() == ',') is.ignore();
    is >> std::ws >> d >> std::ws;
    if (paren) {
        if (is.eof()) return is;
        if (is.peek() == ')') is.ignore();
        else                  return is; // closing parenthesis is missing (TODO: throw an exception)
    }
    if (formatted) {
        if (is.eof()) return is;
        if (is.peek() == '=') is.ignore();
        else                  return is; //  = 0  is missing (TODO: throw an exception)
        is >> std::ws;
        if (is.peek() == '0') is.ignore();
        else                  return is; // = 0 is missing (TODO: throw an exception)
        is >> std::ws;
        if (!paren && is.peek() == '>') is.ignore();
        else                  return is; // closing > is missing (TODO: throw an exception)
    }
    p.set(a,b,c,d);
    return is;
}

template <class T>
void vgl_plane_3d<T>::
plane_coord_vectors(vgl_vector_3d<T>& uvec, vgl_vector_3d<T>& vvec) const
{
    vgl_vector_3d<T> Y((T)0, (T)1, (T)0);
    vgl_vector_3d<T> n = this->normal();
    
    // Since we have an int Template definition, we need to static cast input so VS is happy.
    // Note* currently there are only float and double Template defs. If long double is ever created,
    // this cast will need to get expanded to prevent loss of precision issues.
    T dp = (T)1 - (T)std::fabs(static_cast<double>(dot_product(n, Y)));
    T tol = ((T)1)/((T)10);
    if (dp>tol)//ok to use the Y axis to form the coordinate system
    {
        uvec = normalized(cross_product(Y, n));
        vvec = normalized(cross_product(n, uvec));
    }
    else { // the normal is parallel to the Y axis
        vgl_vector_3d<T> Z((T)0, (T)0, (T)1);
        uvec = normalized(cross_product(n, Z));
        vvec = normalized(cross_product(uvec, n));
    }
}


template <class T>
bool vgl_plane_3d<T>::plane_coords(vgl_point_3d<T> const& p3d,
                                   vgl_point_2d<T>& p2d, T tol) const
{
    // check if point is on the plane
    vgl_point_3d<T> pt_on_plane = vgl_closest_point(p3d, *this);
    double dist = vgl_distance(p3d, pt_on_plane), dtol = static_cast<double>(tol);
    if (dist>dtol)
        return false;
    // use the plane point to compute coordinates
    // construct the axis vectors
    vgl_point_3d<T> origin_pt = vgl_closest_point_origin(*this);
    vgl_vector_3d<T> p = pt_on_plane - origin_pt;
    vgl_vector_3d<T> uvec, vvec;
    this->plane_coord_vectors(uvec, vvec);
    T u = dot_product(uvec, p), v = dot_product(vvec, p);
    p2d.set(u, v);
    return true;
}

template <class T>
vgl_point_3d<T>
vgl_plane_3d<T>::world_coords(vgl_point_2d<T> const& p2d) const
{
    vgl_point_3d<T> origin_pt = vgl_closest_point_origin(*this);
    vgl_vector_3d<T> uvec, vvec;
    this->plane_coord_vectors(uvec, vvec);
    vgl_point_3d<T> p3d = origin_pt + uvec*p2d.x() + vvec*p2d.y();
    return p3d;
}

#endif // vgl_plane_3d_h
