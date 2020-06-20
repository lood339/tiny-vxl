// This is core/vgl/vgl_box_2d.h
#ifndef vgl_box_2d_h
#define vgl_box_2d_h
//:
// \file
// \brief Contains class to represent a cartesian 2D bounding box.
// \author Don Hamilton, Peter Tu
// \date   15 Feb 2000
//
// \verbatim
//  Modifications
//   IMS (Manchester)14 Mar 2001: Tidied up the documentation + added binary_io
//   Amitha Perera   10 Jul 2001: Deprecated get_*() in favour of *(), as agreed in Zurich.
//   Peter Vanroose   5 Oct 2001: Added operator==() and is_empty()
//   Peter Vanroose   6 Oct 2001: Added method add(vgl_point_2d<T>) to enlarge a box
//   Peter Vanroose   7 Oct 2001: Removed deprecated get_*() functions
//   Peter Vanroose     Feb 2002: brief doxygen comment placed on single line
//   Peter Vanroose  12 Sep 2002: Added method add(vgl_box_2d<T>) to enlarge a box
//   Peter Vanroose  22 Apr 2003: Interface change (centroid constructor): now in correspondence with vgl_box_3d<T>
//   Peter Vanroose  13 May 2003: Constructor interface change (backward compat)
//   Peter Vanroose  15 Oct 2003: Removed deprecated constructors without 4th arg
//   Peter Vanroose  16 Oct 2003: Corner pts given to constructor may now be in any order
//   Gamze Tunali    25 Jan 2007: Moved intersect(box1,box2) to vgl_intersection
//   Peter Vanroose  30 Mar 2007: Commented out deprecated intersect() function
//   Peter Vanroose  22 Jul 2009: Moved vgl_intersection() to vgl_intersection.h
// \endverbatim

#include <iosfwd>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>

#ifdef _MSC_VER
#  include <vcl_msvc_warnings.h>
#endif
#include "vgl_fwd.h" // forward declare vgl_point_2d
#include "vgl_point_2d.h"


//: Represents a cartesian 2D box
//  A 2d box with sides aligned with the \a x and \a y axes.
//  Also supports operations required of a bounding box for geometric region
//  tests.
//
//  A box can be empty; this is what the default constructor creates, or what
//  is left after applying the empty() method.  Use the add() methods to enlarge
//  a box, and use the contains() methods to check for inclusion of a point or
//  an other box.
//
//  To make the convex union of two boxes, use box1.add(box2).
//  \verbatim
//                                  MaxPosition
//                    O-------------O
//                    |             |
//                    |             |
//                    |  Centroid   |
//                    |      o      |
//                    |             |
//        Y           |             |
//        |           |             |
//        |           O-------------O
//        |       MinPosition
//        O------X
// \endverbatim
// If you are using a vgl_box_2d<int> to indicate a window on an image, do not forget
// that your axes will be flipped. You could think of the window as follows.
//  \verbatim
//        O------X
//        |       MinPosition
//        |             O-------------O
//        |             |             |
//        Y             |             |
//                      |  Centroid   |
//                      |      o      |
//                      |             |
//                      |             |
//                      |             |
//                      O-------------O
//                               MaxPosition
// \endverbatim
// \sa vgl_box_3d
//  Note: area is not defined on the box class to keep a clean interface
//  see vgl_area<T>
template <class Type>
class vgl_box_2d
{
 public:

  //: Default constructor (creates empty box)
  vgl_box_2d();

  //: Construct using two corner points
  vgl_box_2d(Type const corner1[2],
             Type const corner2[2]);

  //: Construct using two corner points
  vgl_box_2d(vgl_point_2d<Type> const& corner1,
             vgl_point_2d<Type> const& corner2);

  //: Construct using ranges in \a x (first two args) and \a y (last two)
  vgl_box_2d(Type xmin, Type xmax, Type ymin, Type ymax);

  enum point_type { centre=0, min_pos, max_pos };

  //: Construct a box sized width x height at a given reference point.
  //  The box will either be centered at ref_point or will have ref_point
  //  as its min-position or max-position, as specified by the 4th argument.
  vgl_box_2d(Type const ref_point[2],
             Type width, Type height,
             point_type);

  //: Construct a box sized width x height at a given reference point.
  //  The box will either be centered at ref_point or will have ref_point
  //  as its min-position or max-position, as specified by the 4th argument.
  vgl_box_2d(vgl_point_2d<Type> const& ref_point,
             Type width, Type height,
             point_type);

  //: Equality test
  inline bool operator==(vgl_box_2d<Type> const& b) const {
    // All empty boxes are equal:
    if (b.is_empty()) return is_empty();
    return  min_x() == b.min_x() && min_y() == b.min_y()
         && max_x() == b.max_x() && max_y() == b.max_y();
  }

  // Data Access---------------------------------------------------------------

  //: Get width of this box (= \a x dimension)
  Type width() const;
  //: Get height of this box (= \a y dimension)
  Type height() const;

  //: Get "volume" (=area) of this box
  Type volume() const { return width()*height(); }

  //: Get min \a x
  inline Type min_x() const { return min_pos_[0]; }
  //: Get min \a y
  inline Type min_y() const { return min_pos_[1]; }
  //: Get max \a x
  inline Type max_x() const { return max_pos_[0]; }
  //: Get max \a y
  inline Type max_y() const { return max_pos_[1]; }

  //: Get the centroid point
  vgl_point_2d<Type> centroid() const;
  //: Get \a x component of centroid
  Type centroid_x() const;
  //: Get \a y component of centroid
  Type centroid_y() const;

  //: Return lower left corner of box
  vgl_point_2d<Type> min_point() const;

  //: Return upper right corner of box
  vgl_point_2d<Type> max_point() const;

  // Data Control--------------------------------------------------------------

  //: Return true if this box is empty
  inline bool is_empty() const {
    return min_x() > max_x() || min_y() > max_y();
  }

  //: Add a point to this box.
  // Do this by possibly enlarging the box so that the point just falls within the box.
  // Adding a point to an empty box makes it a size zero box only containing p.
  void add(vgl_point_2d<Type> const& p);

  //: Make the convex union of two boxes.
  // Do this by possibly enlarging this box so that the corner points of the
  // given box just fall within the box.
  // Adding an empty box does not change the current box.
  void add(vgl_box_2d<Type> const& b);

  //: Return true iff the point p is inside this box
  bool contains(vgl_point_2d<Type> const& p) const;

  //: Return true iff the corner points of b are inside this box
  bool contains(vgl_box_2d<Type> const& b) const;

  //: Return true if \a (x,y) inside box, ie \a x_min <= \a x <= \a x_max etc
  inline bool contains(Type const& x, Type const& y) const {
    return x >= min_x() && x <= max_x() && y >= min_y() && y <= max_y();
  }

  //: Make the box empty
  void empty();

  //: Set left side of box (other side ordinates unchanged)
  inline void set_min_x(Type m) { min_pos_[0]=m; }
  //: Set bottom of box (other side ordinates unchanged)
  inline void set_min_y(Type m) { min_pos_[1]=m; }
  //: Set right side (other side ordinates unchanged)
  inline void set_max_x(Type m) { max_pos_[0]=m; }
  //: Set top (other side ordinates unchanged)
  inline void set_max_y(Type m) { max_pos_[1]=m; }

  //: Move box so centroid lies at cx (width and height unchanged)
  void set_centroid_x(Type cx);
  //: Move box so centroid lies at cy (width and height unchanged)
  void set_centroid_y(Type cy);

  //: Modify width, retaining centroid at current position
  void set_width(Type width);
  //: Modify height, retaining centroid at current position
  void set_height(Type height);

  //: Add to width and height, centroid unchanged.
  // Will move each side by \p expand / 2.
  void expand_about_centroid(Type expand);
  //: Scale width and height, centroid unchanged.
  void scale_about_centroid(double s);
  //: Scale width and height, keeping scaled position of origin unchanged.
  void scale_about_origin(double s);

  //: Modify bottom left. Top right only changed if necessary to avoid empty box
  void setmin_position(Type const min_position[2]);
  //: Modify top right. Bottom left only changed if necessary to avoid empty box
  void setmax_position(Type const max_position[2]);
  //: Modify bottom left. Top right only changed if necessary to avoid empty box
  void set_min_point(vgl_point_2d<Type> const& min_pt);
  //: Modify top right. Bottom left only changed if necessary to avoid empty box
  void set_max_point(vgl_point_2d<Type> const& max_pt);

  //: Move box so centroid lies at c (width, height unchanged)
  inline void set_centroid(Type const c[2]) { set_centroid_x(c[0]); set_centroid_y(c[1]); }
  //: Move box so centroid lies at c (width, height unchanged)
  inline void set_centroid(vgl_point_2d<Type> const& c) { set_centroid_x(c.x()); set_centroid_y(c.y()); }
  
    //: Write box to stream
    // \relatesalso vgl_box_2d
    friend std::ostream&  operator<<(std::ostream& s, vgl_box_2d<Type> const& p)
    {
        if (p.is_empty())
            return s << "<vgl_box_2d (empty)>";
        else
            return s << "<vgl_box_2d "
            << p.min_pos_[0] << ',' << p.min_pos_[1] << " to "
            << p.max_pos_[0] << ',' << p.max_pos_[1] << '>';
    }
  // INTERNALS-----------------------------------------------------------------
 protected:
  // Data Members--------------------------------------------------------------
  Type min_pos_[2];
  Type max_pos_[2];
};



//: Calculate the bounding box of a sequence of points or boxes.
template <class T, class ITER>
void vgl_box_2d_bounds(ITER begin, ITER end, vgl_box_2d<T>& bounding_box)
{
  for (; begin != end; ++begin)
    bounding_box.add(*begin);
}


// Constructors/Destructor---------------------------------------------------

template <class Type>
vgl_box_2d<Type>::vgl_box_2d()
{
    min_pos_[0]=min_pos_[1]=(Type)1;
    max_pos_[0]=max_pos_[1]=(Type)0; // empty box
}

template <class Type>
vgl_box_2d<Type>::vgl_box_2d(Type const corner1[2],
                             Type const corner2[2])
{
    min_pos_[0]=max_pos_[0]=corner1[0];
    min_pos_[1]=max_pos_[1]=corner1[1];
    this->add(corner2);
}

template <class Type>
vgl_box_2d<Type>::vgl_box_2d(vgl_point_2d<Type> const& corner1,
                             vgl_point_2d<Type> const& corner2)
{
    min_pos_[0]=max_pos_[0]=corner1.x();
    min_pos_[1]=max_pos_[1]=corner1.y();
    this->add(corner2);
}

template <class Type>
vgl_box_2d<Type>::vgl_box_2d(Type xmin, Type xmax, Type ymin, Type ymax)
{
    min_pos_[0]=max_pos_[0]=xmin;
    min_pos_[1]=max_pos_[1]=ymin;
    this->add(vgl_point_2d<Type>(xmax,ymax));
    if (xmin > xmax || ymin > ymax) this->empty();
}

template <class Type>
vgl_box_2d<Type>::vgl_box_2d(Type const ref_point[2],
                             Type w, Type h,
                             typename vgl_box_2d<Type>::point_type t)
{
    if (t == vgl_box_2d<Type>::centre)
    {
        min_pos_[0]=Type(ref_point[0]-0.5*w);
        min_pos_[1]=Type(ref_point[1]-0.5*h);
        max_pos_[0]=Type(ref_point[0]+0.5*w);
        max_pos_[1]=Type(ref_point[1]+0.5*h);
    }
    else if (t == vgl_box_2d<Type>::min_pos)
    {
        min_pos_[0]=ref_point[0];
        min_pos_[1]=ref_point[1];
        max_pos_[0]=ref_point[0]+w;
        max_pos_[1]=ref_point[1]+h;
    }
    else if (t == vgl_box_2d<Type>::max_pos)
    {
        min_pos_[0]=ref_point[0]-w;
        min_pos_[1]=ref_point[1]-h;
        max_pos_[0]=ref_point[0];
        max_pos_[1]=ref_point[1];
    }
    else
        assert(!"point_type should be one of: centre, min_pos, max_pos");
}

template <class Type>
vgl_box_2d<Type>::vgl_box_2d(vgl_point_2d<Type> const& ref_point,
                             Type w, Type h,
                             typename vgl_box_2d<Type>::point_type t)
{
    if (t == vgl_box_2d<Type>::centre)
    {
        min_pos_[0]=Type(ref_point.x()-0.5*w);
        min_pos_[1]=Type(ref_point.y()-0.5*h);
        max_pos_[0]=Type(ref_point.x()+0.5*w);
        max_pos_[1]=Type(ref_point.y()+0.5*h);
    }
    else if (t == vgl_box_2d<Type>::min_pos)
    {
        min_pos_[0]=ref_point.x();
        min_pos_[1]=ref_point.y();
        max_pos_[0]=ref_point.x()+w;
        max_pos_[1]=ref_point.y()+h;
    }
    else if (t == vgl_box_2d<Type>::max_pos)
    {
        min_pos_[0]=ref_point.x()-w;
        min_pos_[1]=ref_point.y()-h;
        max_pos_[0]=ref_point.x();
        max_pos_[1]=ref_point.y();
    }
    else
        assert(!"point_type should be one of: centre, min_pos, max_pos");
}

template <class Type>
Type vgl_box_2d<Type>::centroid_x() const
{
    assert(!is_empty());
    return Type(0.5*(min_pos_[0] + max_pos_[0]));
}

template <class Type>
Type vgl_box_2d<Type>::centroid_y() const
{
    assert(!is_empty());
    return Type(0.5*(min_pos_[1] + max_pos_[1]));
}

template <class Type>
Type vgl_box_2d<Type>::width() const
{
    return (max_pos_[0] > min_pos_[0]) ? max_pos_[0] - min_pos_[0] : 0;
}

template <class Type>
Type vgl_box_2d<Type>::height() const
{
    return (max_pos_[1] > min_pos_[1]) ? max_pos_[1] - min_pos_[1] : 0;
}

template <class Type>
vgl_point_2d<Type> vgl_box_2d<Type>::min_point() const
{
    assert(!is_empty());
    return vgl_point_2d<Type>(min_pos_[0],min_pos_[1]);
}

template <class Type>
vgl_point_2d<Type> vgl_box_2d<Type>::max_point() const
{
    assert(!is_empty());
    return vgl_point_2d<Type>(max_pos_[0],max_pos_[1]);
}

template <class Type>
vgl_point_2d<Type> vgl_box_2d<Type>::centroid() const
{
    assert(!is_empty());
    return vgl_point_2d<Type>(centroid_x(),centroid_y());
}

template <class Type>
void vgl_box_2d<Type>::set_centroid_x(Type cent_x)
{
    assert(!is_empty());
    Type delta = cent_x - centroid_x();
    min_pos_[0]= min_pos_[0] + delta;
    max_pos_[0]= max_pos_[0] + delta;
}

template <class Type>
void vgl_box_2d<Type>::set_centroid_y(Type cent_y)
{
    assert(!is_empty());
    Type delta = cent_y - centroid_y();
    min_pos_[1]= min_pos_[1] + delta;
    max_pos_[1]= max_pos_[1] + delta;
}

template <class T>
inline void set_dim_2d(T & minv, T& maxv, T spread);

// All this code is to avoid drift in the centroid.
template <>
inline void set_dim_2d(int & minv, int& maxv, int spread)
{
    int sum = minv + maxv;
    sum = sum | (spread&1); // if width is odd, then make sum odd
    minv = int(std::floor((sum-spread)/2.0));
    maxv = minv+spread;
}

template <class T>
inline void set_dim_2d(T & minv, T& maxv, T spread)
{
    T x = minv + maxv;
    minv = T( (x-spread)*0.5 );
    maxv = minv + spread;
}

//: Modify width, retaining centroid at current position
// For integer types, centroid might change slightly, but
// repeat calls to set_height will not cause centroid drift.
template <class Type>
void vgl_box_2d<Type>::set_width(Type w)
{
    assert(!is_empty());
    set_dim_2d(min_pos_[0], max_pos_[0], w);
}

//: Modify height, retaining centroid at current position
// For integer types, centroid might change slightly, but
// repeat calls to set_height will not cause centroid drift.
template <class Type>
void vgl_box_2d<Type>::set_height(Type h)
{
    assert(!is_empty());
    set_dim_2d(min_pos_[1], max_pos_[1], h);
}


//: Add to width and height, centroid unchanged.
// Will move each side by \p expand / 2.
template <class Type>
void vgl_box_2d<Type>::expand_about_centroid(Type expand)
{
    assert(!is_empty());
    set_dim_2d(min_pos_[0], max_pos_[0], width() + expand );
    set_dim_2d(min_pos_[1], max_pos_[1], height() + expand );
}

//: Scale width and height, centroid unchanged.
template <class Type>
void vgl_box_2d<Type>::scale_about_centroid(double s)
{
    assert(!is_empty());
    set_dim_2d(min_pos_[0], max_pos_[0], static_cast<Type>(width()*s));
    set_dim_2d(min_pos_[1], max_pos_[1], static_cast<Type>(height()*s));
}


//: Scale width and height, keeping scaled position of origin unchanged.
template <class Type>
void vgl_box_2d<Type>::scale_about_origin(double s)
{
    min_pos_[0] = static_cast<Type>(min_pos_[0] * s);
    min_pos_[1] = static_cast<Type>(min_pos_[1] * s);
    max_pos_[0] = static_cast<Type>(max_pos_[0] * s);
    max_pos_[1] = static_cast<Type>(max_pos_[1] * s);
}

template <class Type>
void vgl_box_2d<Type>::setmin_position(Type const min_position[2])
{
    min_pos_[0]=min_position[0];
    min_pos_[1]=min_position[1];
    if (max_pos_[0] < min_pos_[0]) {
        max_pos_[0]=min_pos_[0];
    }
    if (max_pos_[1] < min_pos_[1]) {
        max_pos_[1]=min_pos_[1];
    }
}

template <class Type>
void vgl_box_2d<Type>::setmax_position(Type const max_position[2])
{
    max_pos_[0]=max_position[0];
    max_pos_[1]=max_position[1];
    if (max_pos_[0] < min_pos_[0])
        min_pos_[0]=max_pos_[0];
    if (max_pos_[1] < min_pos_[1])
        min_pos_[1]=max_pos_[1];
}

template <class Type>
void vgl_box_2d<Type>::set_min_point(vgl_point_2d<Type> const& min_pt)
{
    min_pos_[0]=min_pt.x(); if (max_pos_[0]<min_pos_[0]) max_pos_[0]=min_pos_[0];
    min_pos_[1]=min_pt.y(); if (max_pos_[1]<min_pos_[1]) max_pos_[1]=min_pos_[1];
}

template <class Type>
void vgl_box_2d<Type>::set_max_point(vgl_point_2d<Type> const& max_pt)
{
    max_pos_[0]=max_pt.x(); if (max_pos_[0]<min_pos_[0]) min_pos_[0]=max_pos_[0];
    max_pos_[1]=max_pt.y(); if (max_pos_[1]<min_pos_[1]) min_pos_[1]=max_pos_[1];
}


//: Add a point to this box.
// Do this by possibly enlarging the box so that the point just falls within the box.
// Adding a point to an empty box makes it a size zero box only containing p.
template <class Type>
void vgl_box_2d<Type>::add(vgl_point_2d<Type> const& p)
{
    if (is_empty())
    {
        min_pos_[0] = max_pos_[0] = p.x();
        min_pos_[1] = max_pos_[1] = p.y();
    }
    else
    {
        if (p.x() > max_pos_[0]) max_pos_[0] = p.x();
        if (p.x() < min_pos_[0]) min_pos_[0] = p.x();
        if (p.y() > max_pos_[1]) max_pos_[1] = p.y();
        if (p.y() < min_pos_[1]) min_pos_[1] = p.y();
    }
}

//: Make the convex union of two boxes
// Do this by possibly enlarging this box so that the corner points of the
// given box just fall within the box.
// Adding an empty box does not change the current box.
template <class Type>
void vgl_box_2d<Type>::add(vgl_box_2d<Type> const& b)
{
    if (b.is_empty()) return;
    add(b.min_point());
    add(b.max_point());
}

//: Return true iff the point p is inside this box
template <class Type>
bool vgl_box_2d<Type>::contains(vgl_point_2d<Type> const& p) const
{
    return contains(p.x(), p.y());
}

//: Return true iff the corner points of b are inside this box
template <class Type>
bool vgl_box_2d<Type>::contains(vgl_box_2d<Type> const& b) const
{
    return
    contains(b.min_x(), b.min_y()) &&
    contains(b.max_x(), b.max_y());
}

//: Make the box empty
template <class Type>
void vgl_box_2d<Type>::empty()
{
    min_pos_[0]=min_pos_[1]=(Type)1;
    max_pos_[0]=max_pos_[1]=(Type)0;
}

#endif // vgl_box_2d_h
