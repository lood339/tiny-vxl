// This is core/vgl/vgl_sphere_3d.h
#ifndef vgl_sphere_3d_h
#define vgl_sphere_3d_h
//:
// \file
// \brief a sphere in 3D nonhomogeneous space
// \author Ian Scott

#include <iosfwd>

#include "vgl_fwd.h"// forward declare vgl_line_3d_2_points
#include "vgl_point_3d.h"
#include "vgl_point_3d.h"
#include "vgl_closest_point.h"
#include "vgl_line_3d_2_points.h"

//: Represents a cartesian 3D point
template <class Type>
class vgl_sphere_3d
{
  vgl_point_3d<Type> c_; //!< centre
  Type r_;               //!< radius
 public:

  // Constructors/Initializers/Destructor------------------------------------

  //: Default constructor
  inline vgl_sphere_3d (): c_(Type(0), Type(0), Type(0)), r_(Type(-1)) {}

  //: Construct from four scalars: centre and radius.
  inline vgl_sphere_3d(Type px, Type py, Type pz, Type rad) : c_(px, py, pz), r_(rad) {}

  //: Construct from a 4-array, representing centre and radius.
  inline vgl_sphere_3d (const Type v[4]): c_(v[0], v[1], v[2]), r_(v[3]) {}

  //: Construct from centre point and radius.
  vgl_sphere_3d (vgl_point_3d<Type> const& cntr, Type rad): c_(cntr), r_(rad) {}

  //: Test for equality
  inline bool operator==(const vgl_sphere_3d<Type> &s) const { return this==&s || (c_==s.c_ && r_==s.r_); }
  //: Test for inequality
  inline bool operator!=(vgl_sphere_3d<Type>const& s) const { return !operator==(s); }

  // Data Access-------------------------------------------------------------

  inline const vgl_point_3d<Type> & centre() const {return c_;}
  inline Type radius() const {return r_;}

  //: Return true if this sphere is empty
  inline bool is_empty() const {
    return r_ < 0.0;
  }

  //: Return true iff the point p is inside (or on) this sphere
  bool contains(vgl_point_3d<Type> const& p) const;

  //: Make the sphere empty.
  void set_empty() {c_.set(0,0,0); r_=-1;}

  //: Set radius \a r of this sphere (while centre unchanged)
  inline void set_radius(Type r) { r_=r; }
  //: Set centre of this sphere to \a c (while radius unchanged)
  inline void set_centre(const vgl_point_3d<Type> & c) { c_=c; }

  //: Calculate the end points of a line clipped by this sphere.
  bool clip(const vgl_line_3d_2_points<Type> & line,
            vgl_point_3d<Type> &p1, vgl_point_3d<Type> &p2) const;


  //: convert point on sphere to Cartesian coordinates, angles in radians
  void spherical_to_cartesian(Type elevation_rad, Type azimuth_rad,
                              Type& x, Type& y, Type& z) const;

  void spherical_to_cartesian(Type elevation_rad, Type azimuth_rad,
                              vgl_point_3d<Type>&  pt) const;

  //:find elevation and azimuth of closest point on the sphere to x,y,z
  void cartesian_to_spherical(Type x, Type y, Type z, Type& elevation_rad, Type& azimuth_rad) const;
  void cartesian_to_spherical(vgl_point_3d<Type> const& pt, Type& elevation_rad, Type& azimuth_rad) const;

  //: Writes "<vgl_sphere_3d centre=vgl_point_3d<x,y,z> radius=r)>" to stream
  std::ostream& print(std::ostream& os) const;

};


//: Writes "<vgl_sphere_3d centre=vgl_point_3d<x,y,z> radius=r)>" to stream
template <class Type>
std::ostream& operator<<(std::ostream& os, const vgl_sphere_3d<Type>& sph);

// copy from .cpp
//: Return true iff the point p is inside (or on) this sphere
template <class T>
bool vgl_sphere_3d<T>::contains(vgl_point_3d<T> const& p) const
{
    return r_ >= 0 && (p-c_).sqr_length() <= r_*r_;
}


//: Calculate the end points of a line clipped by this sphere.
// \return true if any of the line touches the sphere.
template <class T>
bool vgl_sphere_3d<T>::clip(const vgl_line_3d_2_points<T> & line,
                            vgl_point_3d<T> &p1, vgl_point_3d<T> &p2) const
{
    // The empty sphere does not intersect anything:
    if (r_ < 0) return false;
    
    vgl_point_3d<T> cp = vgl_closest_point(line, c_);
    
    T cp_sqr_len = (cp - c_).sqr_length();
    if (cp_sqr_len > r_*r_) return false;
    double arg = static_cast<double>(r_*r_ - cp_sqr_len);//for VC10
    T half_chord_len = static_cast<T>(std::sqrt(arg));
    
    vgl_vector_3d<T> linevec = line.direction();
    linevec *= half_chord_len / linevec.length();
    
    p1 = cp - linevec;
    p2 = cp + linevec;
    
    return true;
}


//: Writes "<vgl_sphere_3d centre=vgl_point_3d<x,y,z> radius=r)>" to stream
template <class T>
std::ostream& vgl_sphere_3d<T>::print(std::ostream& os) const
{
    return os << "<vgl_sphere_3d centre=" << c_
    << "radius=" << r_ << '>';
}


//: Read from stream, possibly with formatting.
//  Either just reads 4 blank-separated numbers,
//  or reads 4 comma-separated numbers,


template <class Type>
void vgl_sphere_3d<Type>::spherical_to_cartesian(Type elevation_rad, Type azimuth_rad,
                                                 Type& x, Type& y, Type& z) const{
    
    double el = static_cast<double>(elevation_rad), az = static_cast<double>(azimuth_rad);
    double cx = static_cast<double>(c_.x()),cy =static_cast<double>(c_.y()), cz = static_cast<double>(c_.z());
    double r = static_cast<double>(r_);
    double se = std::sin(el), ce = std::cos(el);
    double sa = std::sin(az), ca = std::cos(az);
    
    x = static_cast<Type>((r*se*ca)+cx);
    y = static_cast<Type>((r*se*sa)+cy);
    z = static_cast<Type>((r*ce)+cz);
    
}

template <class Type>
void vgl_sphere_3d<Type>::spherical_to_cartesian(Type elevation_rad, Type azimuth_rad, vgl_point_3d<Type>&  pt) const
{
    
    Type x, y, z;
    spherical_to_cartesian(elevation_rad, azimuth_rad, x, y, z);
    pt.set(x, y, z);
    
}
template <class Type>
void vgl_sphere_3d<Type>::cartesian_to_spherical(Type x, Type y, Type z, Type& elevation_rad, Type& azimuth_rad) const{
    double xd = static_cast<double>(x-c_.x()), yd = static_cast<double>(y-c_.y()), zd = static_cast<double>(z-c_.z());
    double r  = std::sqrt(xd*xd + yd*yd +zd*zd);
    elevation_rad = static_cast<Type>(std::acos(zd/r));
    azimuth_rad = static_cast<Type>(std::atan2(yd,xd));
}
template <class Type>
void vgl_sphere_3d<Type>::cartesian_to_spherical(vgl_point_3d<Type> const& pt, Type& elevation_rad, Type& azimuth_rad) const{
    return cartesian_to_spherical(pt.x(), pt.y(), pt.z(),elevation_rad,azimuth_rad);
}

//: Writes "<vgl_sphere_3d centre=vgl_point_3d<x,y,z> radius=r)>" to stream
template <class T>
std::ostream& operator<<(std::ostream& os, const vgl_sphere_3d<T>& sph)
{
    return sph.print(os);
}




#endif // vgl_sphere_3d_h
