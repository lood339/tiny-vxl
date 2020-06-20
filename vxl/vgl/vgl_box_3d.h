// This is core/vgl/vgl_box_3d.h
#ifndef vgl_box_3d_h
#define vgl_box_3d_h
//:
// \file
// \brief Contains class to represent a cartesian 3D bounding box.
// \author Don Hamilton, Peter Tu
// \date   15 Feb 2000
//
// \verbatim
//  Modifications
//   Peter Vanroose, 28 Feb 2000: lots of minor corrections
//   NPC (Manchester)14 Mar 2001: Tidied up the documentation + added binary_io
//   Peter Vanroose, 10 Jul 2001: Deprecated get_*() in favour of *(), and explicit casts
//   Peter Vanroose,  5 Oct 2001: Added operator==() and methods is_empty() and contains()
//   Peter Vanroose,  6 Oct 2001: Added method add(vgl_point_3d<T>) to enlarge a box
//   Peter Vanroose,  7 Oct 2001: Removed deprecated get_*() functions
//   Peter Vanroose,    Feb 2002: brief doxygen comment placed on single line
//   Peter Vanroose, 12 Sep 2002: Added method add(vgl_box_3d<T>) to enlarge a box
//   Peter Vanroose, 13 May 2003: Constructor interface change (compat with vgl_box_2d)
//   Peter Vanroose  15 Oct 2003: Removed deprecated constructors without 5th arg
//   Peter Vanroose  16 Oct 2003: Corner pts given to constructor may now be in any order
//   Peter Vanroose  16 Oct 2003: Added intersect(box1,box2)
//   Gamze Tunali    25 Jan 2007: Moved intersect(box1,box2) to vgl_intersection
//   Peter Vanroose  30 Mar 2007: Commented out deprecated intersect() function
//   Peter Vanroose  22 Jul 2009: Moved vgl_intersection() to vgl_intersection.h
// \endverbatim

#include <iosfwd>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif
#include <vgl/vgl_fwd.h> // forward declare vgl_point_3d

//: Represents a cartesian 3D box
//  A 3d box with sides aligned with \a x, \a y and \a z axes. Supports operations
//  required of a bounding box for geometric volume tests.
//
//  A box can be empty; this is what the default constructor creates, or what
//  is left after applying the empty() method.  Use the add() methods to enlarge
//  a box, and use the contains() methods to check for inclusion of a point or
//  an other box.
//
//  To make the convex union of two boxes, use box1.add(box2).
//  \verbatim
//                                 MaxPosition
//                       |<--width-->|
//                       O-----------O  ---
//                      /           /|   ^
//                     /           / |   |
//                    O-----------O  | height
//                    |       o   |  |   |
//                    |  centroid |  |   v
//                    |           |  O  ---
//     Y              |           | /   /_____depth
//     |   Z          |           |/   /
//     |  /           O-----------O  ---
//     | /         MinPosition
//     O-----X
// \endverbatim
// \sa vgl_box_2d

template <class Type>
class vgl_box_3d
{
 public:

  //: Default constructor (creates empty box)
  vgl_box_3d();

  //: Construct using two corner points
  vgl_box_3d(Type const corner1[3],
             Type const corner2[3]);

  //: Construct using two corner points
  vgl_box_3d(vgl_point_3d<Type> const& corner1,
             vgl_point_3d<Type> const& corner2);

  //: Construct from ranges in \a x,y,z (take care with order of inputs).
  //  The \a x range is given by the 1st and 4th coordinates,
  //  the \a y range is given by the 2nd and 5th coordinates,
  //  the \a z range is given by the 3rd and 6th coordinates.
  vgl_box_3d(Type xmin, Type ymin, Type zmin,
             Type xmax, Type ymax, Type zmax);

  enum point_type { centre=0, min_pos, max_pos };

  //: Construct a box sized width x height x depth at a given reference point.
  //  The box will either be centered at ref_point or will have ref_point
  //  as its min-position or max-position, as specified by the 5th argument.
  vgl_box_3d(Type const ref_point[3],
             Type width, Type height, Type depth,
             point_type);

  //: Construct a box sized width x height x depth at a given reference point.
  //  The box will either be centered at ref_point or will have ref_point
  //  as its min-position or max-position, as specified by the 5th argument.
  vgl_box_3d(vgl_point_3d<Type> const& ref_point,
             Type width, Type height, Type depth,
             point_type);

  //: Equality test
  inline bool operator==(vgl_box_3d<Type> const& b) const {
    // All empty boxes are equal:
    if (b.is_empty()) return is_empty();
    return min_x() == b.min_x() && min_y() == b.min_y() && min_z() == b.min_z()
        && max_x() == b.max_x() && max_y() == b.max_y() && max_z() == b.max_z();
  }

  // Data Access---------------------------------------------------------------

  //: Get width of this box (= \a x dimension)
  Type width() const;
  //: Get height of this box (= \a y dimension)
  Type height() const;
  //: Get depth of this box (= \a z dimension)
  Type depth() const;

  //: Get volume of this box
  inline Type volume() const { return width()*height()*depth(); }

  //: Get min \a x
  inline Type min_x() const { return min_pos_[0]; }
  //: Get min \a y
  inline Type min_y() const { return min_pos_[1]; }
  //: Get min \a z
  inline Type min_z() const { return min_pos_[2]; }

  //: Get max \a x
  inline Type max_x() const { return max_pos_[0]; }
  //: Get max \a y
  inline Type max_y() const { return max_pos_[1]; }
  //: Get max \a z
  inline Type max_z() const { return max_pos_[2]; }

  //: Get the centroid point
  vgl_point_3d<Type> centroid() const;
  //: Get \a x component of centroid
  Type centroid_x() const;
  //: Get \a y component of centroid
  Type centroid_y() const;
  //: Get \a z component of centroid
  Type centroid_z() const;

  //: Return lower left corner of box
  vgl_point_3d<Type> min_point() const;

  //: Return upper right corner of box
  vgl_point_3d<Type> max_point() const;

  //: Return the 8 vertices of the box
  std::vector<vgl_point_3d<Type> > vertices() const;

  // Data Control--------------------------------------------------------------

  //: Return true if this box is empty
  inline bool is_empty() const {
    return min_x() > max_x() || min_y() > max_y() || min_z() > max_z();
  }

  //: Add a point to this box.
  // Do this by possibly enlarging the box so that the point just falls within the box.
  // Adding a point to an empty box makes it a size zero box only containing p.
  void add(vgl_point_3d<Type> const& p);

  //: Make the convex union of two boxes.
  // Do this by possibly enlarging this box so that the corner points of the
  // given box just fall within the box.
  // Adding an empty box does not change the current box.
  void add(vgl_box_3d<Type> const& b);

  //: Return true iff the point p is inside this box
  bool contains(vgl_point_3d<Type> const& p) const;

  //: Return true iff the corner points of b are inside this box
  bool contains(vgl_box_3d<Type> const& b) const;

  //: Return true if \a (x,y,z) is inside this box, ie \a x_min <= \a x <= \a x_max etc
  inline bool contains(Type const& x, Type const& y, Type const& z) const {
    return x >= min_x() && x <= max_x() &&
           y >= min_y() && y <= max_y() &&
           z >= min_z() && z <= max_z();
  }

  //: Make the box empty
  void empty();

  //: Set min \a x ordinate of box (other sides unchanged)
  inline void set_min_x(Type m) { min_pos_[0]=m; }
  //: Set min \a y ordinate of box (other sides unchanged)
  inline void set_min_y(Type m) { min_pos_[1]=m; }
  //: Set min \a z ordinate of box (other sides unchanged)
  inline void set_min_z(Type m) { min_pos_[2]=m; }

  //: Set max \a x ordinate of box (other sides unchanged)
  inline void set_max_x(Type m) { max_pos_[0]=m; }
  //: Set max \a y ordinate of box (other sides unchanged)
  inline void set_max_y(Type m) { max_pos_[1]=m; }
  //: Set max \a z ordinate of box (other sides unchanged)
  inline void set_max_z(Type m) { max_pos_[2]=m; }

  //: Move box so centroid lies at cx (size unchanged)
  void set_centroid_x(Type cx);
  //: Move box so centroid lies at cy (size unchanged)
  void set_centroid_y(Type cy);
  //: Move box so centroid lies at cz (size unchanged)
  void set_centroid_z(Type cz);

  //: Set width (x), centroid unchanged
  void set_width(Type width);
  //: Set height (y), centroid unchanged
  void set_height(Type height);
  //: Set depth (z), centroid unchanged
  void set_depth(Type depth);


  //: Add to width and height, centroid unchanged.
  // Will move each side by \p expand / 2.
  void expand_about_centroid(Type expand);
  //: Scale width, height and depth, centroid unchanged.
  void scale_about_centroid(double s);
  //: Scale width, height and depth, keeping scaled position of origin unchanged.
  void scale_about_origin(double s);

  //: Modify min corner point. Max corner point only changed if necessary to avoid empty box
  void set_min_position(Type const m[3]);
  //: Modify max corner point. Min corner point only changed if necessary to avoid empty box
  void set_max_position(Type const m[3]);
  //: Modify min corner point. Max corner point only changed if necessary to avoid empty box
  void set_min_point(vgl_point_3d<Type> const& min_pt);
  //: Modify max corner point. Min corner point only changed if necessary to avoid empty box
  void set_max_point(vgl_point_3d<Type> const& max_pt);
  //: Move box so centroid lies at c (size unchanged)
  inline void set_centroid(Type const c[3]) { set_centroid_x(c[0]); set_centroid_y(c[1]); set_centroid_z(c[2]); }
  //: Move box so centroid lies at c (size unchanged)
  inline void set_centroid(vgl_point_3d<Type> const& c) { set_centroid_x(c.x()); set_centroid_y(c.y()); set_centroid_z(c.z()); }

  // INTERNALS-----------------------------------------------------------------
 protected:
  // Data Members--------------------------------------------------------------
  Type min_pos_[3];
  Type max_pos_[3];
};

//: Write box to stream
// \relatesalso vgl_box_3d
template <class Type>
std::ostream&  operator<<(std::ostream& s, vgl_box_3d<Type> const& p);

//: Calculate the bounding box of a sequence of points or boxes.
template <class T, class ITER>
void vgl_box_3d_bounds(ITER begin, ITER end, vgl_box_3d<T>& bounding_box)
{
  for (; begin != end; ++begin)
    bounding_box.add(*begin);
}

// copy from .cpp
// Constructors/Destructor---------------------------------------------------

template <class Type>
vgl_box_3d<Type>::vgl_box_3d()
{
    min_pos_[0]=min_pos_[1]=min_pos_[2]=(Type)1;
    max_pos_[0]=max_pos_[1]=max_pos_[2]=(Type)0; // empty box
}

template <class Type>
vgl_box_3d<Type>::vgl_box_3d(Type const corner1[3],
                             Type const corner2[3])
{
    min_pos_[0]=max_pos_[0]=corner1[0];
    min_pos_[1]=max_pos_[1]=corner1[1];
    min_pos_[2]=max_pos_[2]=corner1[2];
    this->add(corner2);
}

template <class Type>
vgl_box_3d<Type>::vgl_box_3d(vgl_point_3d<Type> const& corner1,
                             vgl_point_3d<Type> const& corner2)
{
    min_pos_[0]=max_pos_[0]=corner1.x();
    min_pos_[1]=max_pos_[1]=corner1.y();
    min_pos_[2]=max_pos_[2]=corner1.z();
    this->add(corner2);
}

template <class Type>
vgl_box_3d<Type>::vgl_box_3d(Type xmin, Type ymin, Type zmin,
                             Type xmax, Type ymax, Type zmax)
{
    min_pos_[0]=max_pos_[0]=xmin;
    min_pos_[1]=max_pos_[1]=ymin;
    min_pos_[2]=max_pos_[2]=zmin;
    this->add(vgl_point_3d<Type>(xmax,ymax,zmax));
    if (xmin > xmax || ymin > ymax || zmin > zmax) this->empty();
}

template <class Type>
vgl_box_3d<Type>::vgl_box_3d(Type const ref_point[3],
                             Type w, Type h, Type d,
                             typename vgl_box_3d<Type>::point_type t)
{
    if (t == vgl_box_3d<Type>::centre)
    {
        min_pos_[0]=Type(ref_point[0]-0.5*w);
        min_pos_[1]=Type(ref_point[1]-0.5*h);
        min_pos_[2]=Type(ref_point[2]-0.5*d);
        max_pos_[0]=Type(ref_point[0]+0.5*w);
        max_pos_[1]=Type(ref_point[1]+0.5*h);
        max_pos_[2]=Type(ref_point[2]+0.5*d);
    }
    else if (t == vgl_box_3d<Type>::min_pos)
    {
        min_pos_[0]=ref_point[0];
        min_pos_[1]=ref_point[1];
        min_pos_[2]=ref_point[2];
        max_pos_[0]=ref_point[0]+w;
        max_pos_[1]=ref_point[1]+h;
        max_pos_[2]=ref_point[2]+d;
    }
    else if (t == vgl_box_3d<Type>::max_pos)
    {
        min_pos_[0]=ref_point[0]-w;
        min_pos_[1]=ref_point[1]-h;
        min_pos_[2]=ref_point[2]-d;
        max_pos_[0]=ref_point[0];
        max_pos_[1]=ref_point[1];
        max_pos_[2]=ref_point[2];
    }
    else
        assert(!"point_type should be one of: centre, min_pos, max_pos");
}

template <class Type>
vgl_box_3d<Type>::vgl_box_3d(vgl_point_3d<Type> const& ref_point,
                             Type w, Type h, Type d,
                             typename vgl_box_3d<Type>::point_type t)
{
    if (t == vgl_box_3d<Type>::centre)
    {
        min_pos_[0]=Type(ref_point.x()-0.5*w);
        min_pos_[1]=Type(ref_point.y()-0.5*h);
        min_pos_[2]=Type(ref_point.z()-0.5*d);
        max_pos_[0]=Type(ref_point.x()+0.5*w);
        max_pos_[1]=Type(ref_point.y()+0.5*h);
        max_pos_[2]=Type(ref_point.z()+0.5*d);
    }
    else if (t == vgl_box_3d<Type>::min_pos)
    {
        min_pos_[0]=ref_point.x();
        min_pos_[1]=ref_point.y();
        min_pos_[2]=ref_point.z();
        max_pos_[0]=ref_point.x()+w;
        max_pos_[1]=ref_point.y()+h;
        max_pos_[2]=ref_point.z()+d;
    }
    else if (t == vgl_box_3d<Type>::max_pos)
    {
        min_pos_[0]=ref_point.x()-w;
        min_pos_[1]=ref_point.y()-h;
        min_pos_[2]=ref_point.z()-d;
        max_pos_[0]=ref_point.x();
        max_pos_[1]=ref_point.y();
        max_pos_[2]=ref_point.z();
    }
    else
        assert(!"point_type should be one of: centre, min_pos, max_pos");
}

template <class Type>
Type vgl_box_3d<Type>::width() const
{
    return (max_pos_[0] > min_pos_[0]) ? max_pos_[0] - min_pos_[0] : 0;
}

template <class Type>
Type vgl_box_3d<Type>::height() const
{
    return (max_pos_[1] > min_pos_[1]) ? max_pos_[1] - min_pos_[1] : 0;
}

template <class Type>
Type vgl_box_3d<Type>::depth() const
{
    return (max_pos_[2] > min_pos_[2]) ? max_pos_[2] - min_pos_[2] : 0;
}

template <class Type>
vgl_point_3d<Type> vgl_box_3d<Type>::centroid() const
{
    return vgl_point_3d<Type>(centroid_x(),centroid_y(),centroid_z());
}

template <class Type>
Type vgl_box_3d<Type>::centroid_x() const
{
    assert(!is_empty());
    return Type(0.5*(min_pos_[0]+max_pos_[0]));
}

template <class Type>
Type vgl_box_3d<Type>::centroid_y() const
{
    assert(!is_empty());
    return Type(0.5*(min_pos_[1]+max_pos_[1]));
}

template <class Type>
Type vgl_box_3d<Type>::centroid_z() const
{
    assert(!is_empty());
    return Type(0.5*(min_pos_[2]+max_pos_[2]));
}

template <class Type>
void vgl_box_3d<Type>::set_centroid_x(Type cx)
{
    assert(!is_empty());
    Type delta = cx - centroid_x();
    min_pos_[0]= min_pos_[0] + delta;
    max_pos_[0]= max_pos_[0] + delta;
}

template <class Type>
void vgl_box_3d<Type>::set_centroid_y(Type cy)
{
    assert(!is_empty());
    Type delta = cy - centroid_y();
    min_pos_[1]= min_pos_[1] + delta;
    max_pos_[1]= max_pos_[1] + delta;
}

template <class Type>
void vgl_box_3d<Type>::set_centroid_z(Type cz)
{
    assert(!is_empty());
    Type delta = cz - centroid_z();
    min_pos_[2]= min_pos_[2] + delta;
    max_pos_[2]= max_pos_[2] + delta;
}

template <class T>
inline void set_dim_3d(T& minv, T& maxv, T spread);

// All this code is to avoid drift in the centroid.
template <>
inline void set_dim_3d(int& minv, int& maxv, int spread)
{
    int sum = minv + maxv;
    sum = sum | (spread&1); // if width is odd, then make sum odd
    minv = int(std::floor((sum-spread)/2.0));
    maxv = minv+spread;
}

template <class T>
inline void set_dim_3d(T& minv, T& maxv, T spread)
{
    T x = minv + maxv;
    minv = T( (x-spread)*0.5 );
    maxv = minv + spread;
}

template <class Type>
void vgl_box_3d<Type>::set_width(Type w)
{
    assert(!is_empty());
    set_dim_3d(min_pos_[0], max_pos_[0], w);
}

template <class Type>
void vgl_box_3d<Type>::set_height(Type h)
{
    assert(!is_empty());
    set_dim_3d(min_pos_[1], max_pos_[1], h);
}

template <class Type>
void vgl_box_3d<Type>::set_depth(Type d)
{
    assert(!is_empty());
    set_dim_3d(min_pos_[2], max_pos_[2], d);
}

//: Add to width and height, centroid unchanged.
// Will move each side by \p expand / 2.
template <class Type>
void vgl_box_3d<Type>::expand_about_centroid(Type expand)
{
    assert(!is_empty());
    set_dim_3d(min_pos_[0], max_pos_[0], width() + expand );
    set_dim_3d(min_pos_[1], max_pos_[1], height() + expand );
    set_dim_3d(min_pos_[2], max_pos_[2], depth() + expand );
}

//: Scale width, height and depth, centroid unchanged.
template <class Type>
void vgl_box_3d<Type>::scale_about_centroid(double s)
{
    assert(!is_empty());
    set_dim_3d(min_pos_[0], max_pos_[0], static_cast<Type>(width()*s));
    set_dim_3d(min_pos_[1], max_pos_[1], static_cast<Type>(height()*s));
    set_dim_3d(min_pos_[2], max_pos_[2], static_cast<Type>(depth()*s));
}

//: Scale width, height and depth, keeping scaled position of origin unchanged.
template <class Type>
void vgl_box_3d<Type>::scale_about_origin(double s)
{
    min_pos_[0] = static_cast<Type>(min_pos_[0] * s);
    min_pos_[1] = static_cast<Type>(min_pos_[1] * s);
    min_pos_[2] = static_cast<Type>(min_pos_[2] * s);
    max_pos_[0] = static_cast<Type>(max_pos_[0] * s);
    max_pos_[1] = static_cast<Type>(max_pos_[1] * s);
    max_pos_[2] = static_cast<Type>(max_pos_[2] * s);
}

template <class Type>
void vgl_box_3d<Type>::set_min_position(Type const min_position[3])
{
    min_pos_[0]=min_position[0];
    min_pos_[1]=min_position[1];
    min_pos_[2]=min_position[2];
    if (max_pos_[0] < min_pos_[0]) max_pos_[0]=min_pos_[0];
    if (max_pos_[1] < min_pos_[1]) max_pos_[1]=min_pos_[1];
    if (max_pos_[2] < min_pos_[2]) max_pos_[2]=min_pos_[2];
}

template <class Type>
void vgl_box_3d<Type>::set_max_position(Type const max_position[3])
{
    max_pos_[0]=max_position[0];
    max_pos_[1]=max_position[1];
    max_pos_[2]=max_position[2];
    if (max_pos_[0] < min_pos_[0]) min_pos_[0]=max_pos_[0];
    if (max_pos_[1] < min_pos_[1]) min_pos_[1]=max_pos_[1];
    if (max_pos_[2] < min_pos_[2]) min_pos_[2]=max_pos_[2];
}

template <class Type>
void vgl_box_3d<Type>::set_min_point(vgl_point_3d<Type> const& min_pt)
{
    min_pos_[0]=min_pt.x(); if (max_pos_[0]<min_pos_[0]) max_pos_[0]=min_pos_[0];
    min_pos_[1]=min_pt.y(); if (max_pos_[1]<min_pos_[1]) max_pos_[1]=min_pos_[1];
    min_pos_[2]=min_pt.z(); if (max_pos_[2]<min_pos_[2]) max_pos_[2]=min_pos_[2];
}

template <class Type>
void vgl_box_3d<Type>::set_max_point(vgl_point_3d<Type> const& max_pt)
{
    max_pos_[0]=max_pt.x(); if (max_pos_[0]<min_pos_[0]) min_pos_[0]=max_pos_[0];
    max_pos_[1]=max_pt.y(); if (max_pos_[1]<min_pos_[1]) min_pos_[1]=max_pos_[1];
    max_pos_[2]=max_pt.z(); if (max_pos_[2]<min_pos_[2]) min_pos_[2]=max_pos_[2];
}



template <class Type>
vgl_point_3d<Type> vgl_box_3d<Type>::min_point() const
{
    assert(!is_empty());
    return vgl_point_3d<Type>(min_pos_[0],min_pos_[1],min_pos_[2]);
}

template <class Type>
vgl_point_3d<Type> vgl_box_3d<Type>::max_point() const
{
    assert(!is_empty());
    return vgl_point_3d<Type>(max_pos_[0],max_pos_[1],max_pos_[2]);
}

template <class Type>
std::vector<vgl_point_3d<Type> > vgl_box_3d<Type>::vertices() const
{
    assert(!is_empty());
    std::vector<vgl_point_3d<Type> > vertices;
    vertices.push_back(vgl_point_3d<Type>(min_pos_[0], min_pos_[1], min_pos_[2]));
    vertices.push_back(vgl_point_3d<Type>(max_pos_[0], min_pos_[1], min_pos_[2]));
    vertices.push_back(vgl_point_3d<Type>(max_pos_[0], max_pos_[1], min_pos_[2]));
    vertices.push_back(vgl_point_3d<Type>(min_pos_[0], max_pos_[1], min_pos_[2]));
    vertices.push_back(vgl_point_3d<Type>(min_pos_[0], min_pos_[1], max_pos_[2]));
    vertices.push_back(vgl_point_3d<Type>(max_pos_[0], min_pos_[1], max_pos_[2]));
    vertices.push_back(vgl_point_3d<Type>(max_pos_[0], max_pos_[1], max_pos_[2]));
    vertices.push_back(vgl_point_3d<Type>(min_pos_[0], max_pos_[1], max_pos_[2]));
    return vertices;
}

//: Add a point to this box.
// Do this by possibly enlarging the box so that the point just falls within the box.
// Adding a point to an empty box makes it a size zero box only containing p.
template <class Type>
void vgl_box_3d<Type>::add(vgl_point_3d<Type> const& p)
{
    if (is_empty())
    {
        min_pos_[0] = max_pos_[0] = p.x();
        min_pos_[1] = max_pos_[1] = p.y();
        min_pos_[2] = max_pos_[2] = p.z();
    }
    else
    {
        if (p.x() > max_pos_[0]) max_pos_[0] = p.x();
        if (p.x() < min_pos_[0]) min_pos_[0] = p.x();
        if (p.y() > max_pos_[1]) max_pos_[1] = p.y();
        if (p.y() < min_pos_[1]) min_pos_[1] = p.y();
        if (p.z() > max_pos_[2]) max_pos_[2] = p.z();
        if (p.z() < min_pos_[2]) min_pos_[2] = p.z();
    }
}

//: Make the convex union of two boxes
// Do this by possibly enlarging this box so that the corner points of the
// given box just fall within the box.
// Adding an empty box does not change the current box.
template <class Type>
void vgl_box_3d<Type>::add(vgl_box_3d<Type> const& b)
{
    if (b.is_empty()) return;
    add(b.min_point());
    add(b.max_point());
}

//: Return true iff the point p is inside this box
template <class Type>
bool vgl_box_3d<Type>::contains(vgl_point_3d<Type> const& p) const
{
    return contains(p.x(), p.y(), p.z());
}

//: Return true iff the corner points of b are inside this box
template <class Type>
bool vgl_box_3d<Type>::contains(vgl_box_3d<Type> const& b) const
{
    return
    contains(b.min_x(), b.min_y(), b.min_z()) &&
    contains(b.max_x(), b.max_y(), b.max_z());
}

//: Make the box empty
template <class Type>
void vgl_box_3d<Type>::empty()
{
    min_pos_[0]=min_pos_[1]=min_pos_[2]=(Type)1;
    max_pos_[0]=max_pos_[1]=max_pos_[2]=(Type)0;
}

//: Write box to stream
template <class Type>
std::ostream&  operator<<(std::ostream& s, vgl_box_3d<Type> const& p)
{
    return p.print(s);
}

//: Read box from stream
template <class Type>
std::istream&  operator>>(std::istream& is,  vgl_box_3d<Type>& p)
{
    return p.read(is);
}

#endif // vgl_box_3d_h
