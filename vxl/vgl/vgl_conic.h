// This is core/vgl/vgl_conic.h
#ifndef vgl_conic_h_
#define vgl_conic_h_
//:
// \file
// \brief A quadratic plane curve
//
//  This example tells you the type of the given conic equation,
//  and prints the equation in readable form:
// \code
//   vgl_conic<double> c(1, 0, 2, 0, 0, -3);
//   std::cout << c.real_type() << '\n'; // prints "real ellipse"
//   std::cout << c << '\n'; // prints the equation: X^2 + 2 Y^2 - 3 = 0
// \endcode
//
// \verbatim
// Modifications
//   Peter Vanroose, 10 sep 1996 wrote description and example file.
//   Peter Vanroose, 17 jun 1998 added PolarLine() and PolarPoint().
//   Peter Vanroose, 18 jun 1998 added Hyperbola and Circle interface.
//   Peter Vanroose, 19 jun 1998 added dual space functions.
//   Peter Vanroose, 21 jun 1998 added Parabola interface.
//   Peter Vanroose, 27 jun 1998 added ComputeParabolaParameters().
//   M.Vergauwen & P.Vanroose, 4 jul 1998 added Intersect() & CommonTangents()
//   Peter Vanroose, 29 aug 2001 ported from Geometry to vgl
//   Peter Vanroose, 30 aug 2001 complete rewrite of most of the code
//   Peter Vanroose, 31 aug 2001 added extensive testing + fixed some bugs
//   Feb.2002 - Peter Vanroose - brief doxygen comment placed on single line
//   Ricardo Fabbri, 08 nov 2008 added curvature_at() method.
// \endverbatim
//
//-----------------------------------------------------------------------------

#include <list>
#include <string>
#include <vector>
#include <iosfwd>
#include <cmath>
#include <iostream>

#include "vgl_homg_point_2d.h"
#include "vgl_homg_line_2d.h"

//: A quadratic plane curve
//
//  A conic is either an ellipse (or circle), a hyperbola, or a parabola.
//  It is represented by a quadratic equation in two nonhomogeneous
//  or three homogeneous coordinates.  Conversely, every quadratic
//  equation represents a conic, be it that it can be degenerate:
//  either in two (intersecting or parallel) lines, or in two
//  coincident lines.  Also, it can have no "visible", real points,
//  when it is an imaginary ellipse, or consist of two complementary
//  imaginary lines in which case it only has one real point, which could
//  still be at infinity.
//
//  These 11 cases are the possible values of vgl_conic::real_type().
//  The default constructor sets the type to "invalid conic";
//  otherwise the correct type is automatically set when the equation
//  of the conic is given to the constructor that takes 6 numeric values
//  (a,b,c,d,e,f): the cartesian equation is then
//  $ax^2 + bxy + cy^2 + dx + ey + f = 0$; the homogeneous equation is
//  $ax^2 + bxy + cy^2 + dxw + eyw + fw^2 = 0$.  (Sometimes with $z$ for $w$.)
//  The numeric type (typically double or float) is the template argument
//  of this class.
//
//  When the conic is degenerate and consists of two lines, the method
//  components() returns a list of two (possibly identical) lines.
//  Otherwise, this method returns an empty list.

template <class T>
class vgl_conic
{
 public:
  enum vgl_conic_type {
    no_type=0,
    real_ellipse,
    real_circle,
    imaginary_ellipse,
    imaginary_circle,
    hyperbola,
    parabola,
    real_intersecting_lines,
    complex_intersecting_lines,
    real_parallel_lines,
    complex_parallel_lines,
    coincident_lines,
    num_conic_types // is here to enable iterating through this list
  };
    
  static const std::vector<std::string> vgl_conic_name;
 private:
  // DATA MEMBERS

  vgl_conic_type type_;
  T a_; //!< coefficient of \a x^2
  T b_; //!< coefficient of \a xy
  T c_; //!< coefficient of \a y^2
  T d_; //!< coefficient of \a xw
  T e_; //!< coefficient of \a yw
  T f_; //!< coefficient of \a w^2

 public:
  inline vgl_conic_type type() const { return type_; }

  //: Returns the type of the conic as a string.
  // Possible returned strings are:
  // "real ellipse", "real circle", "imaginary ellipse", "imaginary circle",
  // "hyperbola", "parabola",
  // "real intersecting lines", "complex intersecting lines",
  // "real parallel lines", "complex parallel lines", "coincident lines".
  // The default constructor sets the type to "invalid conic".
  std::string real_type() const;

  //: Returns the internal enum value corresponding to the string argument.
  // Useful for comparison purposes, or for use in "case" statements.
  static vgl_conic_type type_by_name(std::string const& name);

  //: Converts the conic type from enum (internal representation) to string.
  static std::string type_by_number(vgl_conic_type type);

  //: Returns the coefficient of \f$X^2\f$
  inline T a() const { return  a_; }

  //: Returns the coefficient of \f$XY\f$
  inline T b() const { return  b_; }

  //: Returns the coefficient of \f$Y^2\f$
  inline T c() const { return  c_; }

  //: Returns the coefficient of \f$XW\f$
  inline T d() const { return  d_; }

  //: Returns the coefficient of \f$YW\f$
  inline T e() const { return  e_; }

  //: Returns the coefficient of \f$W^2\f$
  inline T f() const { return  f_; }

  // CONSTRUCTORS AND RELATED STUFF

  // default constructor
  vgl_conic() : type_(no_type) {}
    
  //: constructor using polynomial coefficients.
  //  The order of the coefficients is: $X^2$, $XY$, $Y^2$, $XW$, $YW$, $W^2$,
  //  where $W$ is the homogeneous coordinate (sometimes denoted by $Z$).
  vgl_conic(T a, T b, T c, T d, T e, T f);

  //: constructor using polynomial coefficients, given as a C array.
  //  The order of the coefficients is: $X^2$, $XY$, $Y^2$, $XW$, $YW$, $W^2$,
  //  where $W$ is the homogeneous coordinate (sometimes denoted by $Z$).
  vgl_conic(T const coeff[]);

  //: constructor using centre, signed radii, and angle.
  //  This constructor can only be used for non-degenerate, real
  //  conics: If the centre point c is a finite point and rx and ry
  //  have the same sign, an ellipse is defined (any ellipse can
  //  uniquely be specified this way); rx is the length of one main
  //  axis, ry of the other axis.  Hyperbolas are obtained if rx and
  //  ry have opposite sign; the positive one determines the distance
  //  from bots tops to the centre, and the other one specified the
  //  'minor' axis length. The rotation is about the centre of the
  //  ellipse or hyperbola, measured counterclockwise from the X axis.
  //  A parabola is obtained when the centre has w()=0,
  //  i.e., is a point at infinity. In that case (rx,ry) is the
  //  top, and theta is an eccentricity parameter (since the centre
  //  already specifies the direction of the symmetry axis).
  vgl_conic(vgl_homg_point_2d<T> const& c, T rx, T ry, T theta);

  //: set or reset the conic using polynomial coefficients.
  //  The order of the coefficients is: $X^2$, $XY$, $Y^2$, $XW$, $YW$, $W^2$,
  //  where $W$ is the homogeneous coordinate (sometimes denoted by $Z$).
  void set(T a, T b, T c, T d, T e, T f);

  //: comparison operator.
  //  Comparison is on the conic, not the equation coefficients.  Hence two
  //  conics are identical if their coefficient vectors are multiples of
  //  each other.
  bool operator==(vgl_conic<T> const& c) const;

  // UTILITY FUNCTIONS

  //: Returns true if this conic is degenerate, i.e., if it consists of 2 lines.
  bool is_degenerate() const;

  //: Returns true if a central conic, i.e., an ellipse, circle, or hyperbola.
  //  Also the degenerate versions of these return true.
  //  Returns false if a parabola or two parallel or coinciding lines.
  bool is_central() const;

  //: Returns true if the point pt belongs to the conic.
  //  I.e., if it \e exactly satisfies the conic equation.
  //  Beware of rounding for floating point type T! An "almost" returns false!
  bool contains(vgl_homg_point_2d<T> const& pt) const;

  //: Returns the list of component lines, when degenerate and real components.
  //  Otherwise returns an empty list.
  //  If two coinciding lines, the list contains two identical elements.
  //  Hence this list always has length 0 or 2.
  std::list<vgl_homg_line_2d<T> > components() const;

  // Elementary geometric functions ----------------------------------

  //: Returns the polar line of the given point, w.r.t. this conic.
  //  For a non-degenerate conic, the polar line of a point outside of the conic
  //  is the connection line of the two points on the conic that form the conic
  //  "contour" as seen from that point, i.e., the touch points of the two
  //  tangents to the conic going through the given point.
  //
  //  For a point on the conic, it is just the tangent in that point.
  //
  //  And for a point inside the conic, it is the set of all polar points of
  //  the lines through the given point. This set happens to be a straight line.
  vgl_homg_line_2d<T> polar_line(vgl_homg_point_2d<T> const& p) const;

  //: Returns the polar point of the given line, w.r.t. this conic.
  //  For a non-degenerate conic, the polar point of a line that intersects the
  //  conic in two points is the intersection point of the two tangent lines
  //  though those two points.  Hence it is the point of which this line is
  //  the polar line.
  //
  //  For a tangent line to the conic, it is just the tangent point.
  //
  //  And for a line not intersecting the conic, it is the common intersection
  //  point (inside the conic) of the polar lines of all points of that line.
  vgl_homg_point_2d<T>  polar_point(vgl_homg_line_2d<T> const& l) const;

  //: Returns the tangent to the conic in the point p, if p is on the conic.
  //  In general, returns the polar line of the point w.r.t. the conic.
  vgl_homg_line_2d<T> tangent_at(vgl_homg_point_2d<T> const& p) const { return polar_line(p); }

  //: Returns the centre of the conic, or its point at infinity if a parabola.
  // When two intersecting or parallel lines, returns their intersection point.
  // In all cases this is the polar point of the line at infinity.
  vgl_homg_point_2d<T> centre() const { return polar_point(vgl_homg_line_2d<T>(0,0,1)); }

  //: Returns the curvature of the conic at point p, assuming p is on the conic.
  double curvature_at(vgl_point_2d<T> const& p) const;

  //: Converts the coefficients to a geometric description of an ellipse.
  //  Returns false if the conic is not an ellipse. Double is appropriate
  //  since integer coefficients can produce non-integer ellipse parameters.
  bool ellipse_geometry(double& xc, double& yc, double& major_axis_length,
                        double& minor_axis_length, double& angle_in_radians) const;

  // Functions related to dual space ---------------------------------

  //: Returns the dual or tangential representation of this conic.
  //  The homogeneous coordinates of the points belonging to the dual conic
  //  are the coefficients of the equations of all tangents to the original
  //  conic.
  vgl_conic dual_conic() const;

  //: Returns the dual or tangential representation of this conic.
  vgl_conic tangential_form() const { return dual_conic(); }

  //: Modify this conic by translating it over distance \a x in the \a X direction and distance \a y in the \a Y direction.
  void translate_by(T x, T y);

 private:
  //--------------------------------------------------------------------------
  //: set conic type from polynomial coefficients and store in member type_
  // This method must be called by all constructors (except the default
  // constructor) and all methods that change the coefficients.
  void set_type_from_equation();
};

//: Write "<vgl_conic aX^2+bXY+cY^2+dXW+eYW+fW^2>" to stream
// \relatesalso vgl_conic
template <class T>
std::ostream&  operator<<(std::ostream& s, vgl_conic<T> const& c);

// copy from .cpp file

//--------------------------------------------------------------
//: Returns the type name of the conic.

template <class T>
const std::vector<std::string> vgl_conic<T>::vgl_conic_name = std::vector<std::string>({
    "invalid conic",
    "real ellipse",
    "real circle",
    "imaginary ellipse",
    "imaginary circle",
    "hyperbola",
    "parabola",
    "real intersecting lines",
    "complex intersecting lines",
    "real parallel lines",
    "complex parallel lines",
    "coincident lines"});

template <class T>
std::string vgl_conic<T>::real_type() const { return vgl_conic_name[(int)type_]; }

template <class T>
typename vgl_conic<T>::vgl_conic_type vgl_conic<T>::type_by_name(std::string const& name)
{
    for (int i = (int)no_type; i < num_conic_types; i++)
        if (name == vgl_conic_name[i])
            return (typename vgl_conic<T>::vgl_conic_type)i;
    return no_type; // should never reach this point
}

template <class T>
std::string vgl_conic<T>::type_by_number(typename vgl_conic<T>::vgl_conic_type type)
{
    if (type <= 0 || type >= num_conic_types) return vgl_conic_name[no_type];
    return vgl_conic_name[type];
}

//-------------------------------------------------------------
//: equality test
template <class T>
bool vgl_conic<T>::operator==(vgl_conic<T> const& that) const
{
    if ( type() != that.type() ) return false;
    return   a()*that.b() == b()*that.a()
    && a()*that.c() == c()*that.a()
    && a()*that.d() == d()*that.a()
    && a()*that.e() == e()*that.a()
    && a()*that.f() == f()*that.a()
    && b()*that.c() == c()*that.b()
    && b()*that.d() == d()*that.b()
    && b()*that.e() == e()*that.b()
    && b()*that.f() == f()*that.b()
    && c()*that.d() == d()*that.c()
    && c()*that.e() == e()*that.c()
    && c()*that.f() == f()*that.c()
    && d()*that.e() == e()*that.d()
    && d()*that.f() == f()*that.d()
    && e()*that.f() == f()*that.e();
}

//-------------------------------------------------------------
//: set values
template <class T>
void vgl_conic<T>::set(T ta, T tb, T tc, T td, T te, T tf)
{
    a_ = ta; b_ = tb; c_ = tc; d_ = td; e_ = te; f_ = tf;
    set_type_from_equation();
}

//-------------------------------------------------------------
//: constructor using polynomial coefficients.
template <class T>
vgl_conic<T>::vgl_conic(T const co[])
: type_(no_type), a_(co[0]), b_(co[1]), c_(co[2]), d_(co[3]), e_(co[4]), f_(co[5])
{
    set_type_from_equation();
}

//-------------------------------------------------------------
//: constructor using polynomial coefficients.
template <class T>
vgl_conic<T>::vgl_conic(T ta, T tb, T tc, T td, T te, T tf)
: type_(no_type), a_(ta), b_(tb), c_(tc), d_(td), e_(te), f_(tf)
{
    set_type_from_equation();
}

//: ctor using centre, signed radii, and angle, or (for parabola) top + eccentricity
template <class T>
vgl_conic<T>::vgl_conic(vgl_homg_point_2d<T> const& co, T rx, T ry, T theta)
{
    if (co.w() == 0) { // This is a parabola
        a_ = co.y()*co.y();
        b_ = -2*co.x()*co.y();
        c_ = co.x()*co.x();
        // polar line of (rx,ry) must have direction (co.y(),-co.x()), hence
        // 2*a_*rx + b_*ry + d_ = 2*t*co.x() and b_*rx + 2*c_*ry +e_ = 2*t*co.y() :
        theta /= std::sqrt(co.x()*co.x()+co.y()*co.y()); // cannot be 0
        d_ = -2*a_*rx - b_*ry + 2*theta*co.x();
        e_ = -2*c_*ry - b_*rx + 2*theta*co.y();
        // conic must go through (rx,ry):
        f_ = -a_*rx*rx-b_*rx*ry-c_*ry*ry-d_*rx-e_*ry;
    }
    else { // hyperbola or ellipse
        rx = (rx < 0) ? (-rx*rx) : rx*rx; // abs of square
        ry = (ry < 0) ? (-ry*ry) : ry*ry; // idem
        
        double ct = std::cos(-theta); // rotate counterclockwise over theta
        double st = std::sin(-theta);
        T u = co.x();
        T v = co.y();
        a_ = T(rx*st*st + ry*ct*ct);
        b_ = T(2*(rx-ry)*ct*st);
        c_ = T(rx*ct*ct + ry*st*st);
        d_ = T(-2*(rx*st*st + ry*ct*ct)*u - 2*(rx-ry)*ct*st*v);
        e_ = T(-2*(rx-ry)*ct*st*u - 2*(rx*ct*ct + ry*st*st)*v);
        f_ = T((rx*st*st +ry*ct*ct)*u*u + 2*(rx-ry)*ct*st*u*v + (rx*ct*ct + ry*st*st)*v*v - rx*ry);
    }
    set_type_from_equation();
}

//--------------------------------------------------------------------------
//: set conic type from polynomial coefficients.
// This method must be called by all constructors (except the default
// constructor) and all methods that change the coefficients.

template <class T>
void vgl_conic<T>::set_type_from_equation()
{
    T A = a_, B = b_/2, C = c_, D = d_/2, E = e_/2, F = f_;
    
    /* determinant, subdeterminants and trace values */
    T det = A*(C*F - E*E) - B*(B*F - D*E) + D*(B*E - C*D); // determinant
    T J = A*C - B*B;  // upper 2x2 determinant
    T K = (C*F - E*E) + (A*F - D*D); // sum of two other 2x2 determinants
    T I = A + C; // trace of upper 2x2
    
    if (det != 0) {
        if (J > 0) {
            if (det*I < 0) {
                if (A==C && B==0)      type_ = real_circle;
                else                   type_ = real_ellipse;
            }
            else {
                if (A==C && B==0)      type_ = imaginary_circle;
                else                   type_ = imaginary_ellipse;
            }
        }
        else if (J < 0)            type_ = hyperbola;
        else /* J == 0 */          type_ = parabola;
    }
    else {    // limiting cases
        if (J < 0)                 type_ = real_intersecting_lines;
        else if (J > 0)            type_ = complex_intersecting_lines;
        else /* J == 0 */ {
            if ( A == 0 && B == 0 && C == 0 ) { // line at infinity is component
                if ( D !=0 || E != 0 ) type_ = real_intersecting_lines;
                else if (F != 0)       type_ = coincident_lines; // 2x w=0
                else                   type_ = no_type; // all coefficients are 0
            }
            else if (K < 0)          type_ = real_parallel_lines;
            else if (K > 0)          type_ = complex_parallel_lines;
            else                     type_ = coincident_lines;
        }
    }
}

template <class T>
bool vgl_conic<T>::contains(vgl_homg_point_2d<T> const& p) const
{
    return p.x()*p.x()*a_+p.x()*p.y()*b_+p.y()*p.y()*c_+p.x()*p.w()*d_+p.y()*p.w()*e_+p.w()*p.w()*f_ == 0;
}

template <class T>
bool vgl_conic<T>::is_degenerate() const
{
    T A = a_, B = b_/2, C = c_, D = d_/2, E = e_/2, F = f_;
    T det = A*(C*F - E*E) - B*(B*F - D*E) + D*(B*E - C*D);
    return det==0;
}

//: return a geometric description of the conic if an ellipse
// The centre of the ellipse is (xc, yc)
template <class T>
bool vgl_conic<T>::
ellipse_geometry(double& xc, double& yc, double& major_axis_length,
                 double& minor_axis_length, double& angle_in_radians) const
{
    if (type_!=real_ellipse && type_ != real_circle)
        return false;
    
    // Cast to double and half the non-diagonal (non-quadratic) entries B, D, E.
    double A = static_cast<double>(a_), B = static_cast<double>(b_)*0.5,
    C = static_cast<double>(c_), D = static_cast<double>(d_)*0.5,
    F = static_cast<double>(f_), E = static_cast<double>(e_)*0.5;
    if (A < 0)
        A=-A, B=-B, C=-C, D=-D, E=-E, F=-F;
    
    double det = A*(C*F - E*E) - B*(B*F- D*E) + D*(B*E-C*D);
    double D2 =  A*C - B*B;
    xc = (E*B - C*D)/D2;
    yc = (D*B - A*E)/D2;
    
    double trace = A + C;
    double disc = std::sqrt(trace*trace - 4.0*D2);
    double cmaj = (trace+disc)*D2/(2*det); if (cmaj < 0) cmaj = -cmaj;
    double cmin = (trace-disc)*D2/(2*det); if (cmin < 0) cmin = -cmin;
    minor_axis_length = 1.0/std::sqrt(cmaj>cmin?cmaj:cmin);
    major_axis_length = 1.0/std::sqrt(cmaj>cmin?cmin:cmaj);
    
    // Find the angle that diagonalizes the upper 2x2 sub-matrix
    angle_in_radians  = -0.5 * std::atan2(2*B, C-A);
    //                  ^
    // and return the negative of this angle
    return true;
}


//: Returns the list of component lines, when degenerate and real components.
//  Otherwise returns an empty list.
template <class T>
std::list<vgl_homg_line_2d<T> > vgl_conic<T>::components() const
{
    if (!is_degenerate() ||
        type() == complex_intersecting_lines ||
        type() == complex_parallel_lines)
        return std::list<vgl_homg_line_2d<T> >(); // no real components
    
    T A = a_, B = b_/2, C = c_, D = d_/2, E = e_/2, F = f_;
    
    if (type() == coincident_lines) {
        // coincident lines: rank of the matrix of this conic must be 1
        vgl_homg_line_2d<T> l;
        if (A!=0 || B!=0 || D!=0)
            l = vgl_homg_line_2d<T>(A,B,D);
        else if (C!=0 || E!=0)
            l = vgl_homg_line_2d<T>(B,C,E);
        else // only F!=0 : 2x line at infinity w=0
            l = vgl_homg_line_2d<T>(D,E,F);
        return std::list<vgl_homg_line_2d<T> >(2,l); // two identical lines
    }
    
    // Both component lines must pass through the centre of this conic
    vgl_homg_point_2d<T> cntr = centre();
    
    if (type() == real_parallel_lines)
    {
        // In this case the centre lies at infinity.
        // Either these lines both intersect the X axis, or both intersect the Y axis:
        if (A!=0 || D!=0) { // X axis: intersections satisfy y=0 && Axx+2Dxw+Fww=0:
            vgl_homg_line_2d<T> l1(cntr, vgl_homg_point_2d<T>(-D+std::sqrt(D*D-A*F),0,A)),
            l2(cntr, vgl_homg_point_2d<T>(-D-std::sqrt(D*D-A*F),0,A));
            std::list<vgl_homg_line_2d<T> > v(1,l1); v.push_back(l2);
            return v;
        }
        else { // Y axis: x=0 && Cyy+2Eyw+Fww=0:
            vgl_homg_line_2d<T> l1(cntr, vgl_homg_point_2d<T>(0,-E+std::sqrt(E*E-C*F),C)),
            l2(cntr, vgl_homg_point_2d<T>(0,-E-std::sqrt(E*E-C*F),C));
            std::list<vgl_homg_line_2d<T> > v(1,l1); v.push_back(l2);
            return v;
        }
    }
    
    // Only remaining case: type() == real_intersecting_lines.
    if (A==0 && B==0 && C==0) { // line at infinity (w=0) is a component
        std::list<vgl_homg_line_2d<T> > v(1,vgl_homg_line_2d<T>(0,0,1));
        v.push_back(vgl_homg_line_2d<T>(d_,e_,f_));
        return v;
    }
    // If not, the two component lines are determined by cntr and their direction,
    // i.e., they pass through one of the two pts satisfying w=0 && Axx+2Bxy+Cyy=0:
    if (A==0 && C==0) { // components are vertical and horizontal, resp.
        vgl_homg_line_2d<T> l1(cntr, vgl_homg_point_2d<T>(0,1,0)),
        l2(cntr, vgl_homg_point_2d<T>(1,0,0));
        std::list<vgl_homg_line_2d<T> > v(1,l1); v.push_back(l2);
        return v;
    }
    else {
        vgl_homg_line_2d<T> l1(cntr, vgl_homg_point_2d<T>(-B+std::sqrt(B*B-A*C),A,0)),
        l2(cntr, vgl_homg_point_2d<T>(-B-std::sqrt(B*B-A*C),A,0));
        std::list<vgl_homg_line_2d<T> > v(1,l1); v.push_back(l2);
        return v;
    }
}

//: Return true if a central conic
//  (This is an affine property, not a projective one.)
//  Equivalent to saying that the line at infinity does not touch the conic.
template <class T>
bool vgl_conic<T>::is_central() const
{
    return type_ == real_ellipse|| type_ == imaginary_ellipse|| type_ == hyperbola
    || type_ == real_circle || type_ == imaginary_circle
    || type_ == real_intersecting_lines|| type_ == complex_intersecting_lines;
}

//--------------------------------------------------------------------------------
template <class T>
void vgl_conic<T>::translate_by(T x, T y)
{
    d_ += 2*a_*x + b_*y;
    f_ += c_ * y*y - a_ * x*x + d_ * x + e_ * y;
    e_ += 2*c_*y + b_*x;
    // This does not change the type, so no need to run set_type_from_equation()
}

template <class T>
vgl_conic<T> vgl_conic<T>::dual_conic() const
{
    T A = a_, B = b_/2, C = c_, D = d_/2, E = e_/2, F = f_;
    return vgl_conic<T>(E*E-C*F, 2*(B*F-D*E), D*D-A*F, 2*(C*D-B*E), 2*(A*E-B*D), B*B-A*C);
}

//: return the polar line of the homogeneous 2-D point p.
template <class T>
vgl_homg_line_2d<T> vgl_conic<T>::polar_line(vgl_homg_point_2d<T> const& p) const
{
    return vgl_homg_line_2d<T> (p.x()*a_  +p.y()*b_/2+p.w()*d_/2,
                                p.x()*b_/2+p.y()*c_  +p.w()*e_/2,
                                p.x()*d_/2+p.y()*e_/2+p.w()*f_  );
}

//: return the polar point of the homogeneous 2-D line l.
template <class T>
vgl_homg_point_2d<T> vgl_conic<T>::polar_point(vgl_homg_line_2d<T> const& l) const
{
    if (!is_degenerate()) {
        vgl_conic<T> co = this->dual_conic();
        return vgl_homg_point_2d<T> (l.a()*co.a()  +l.b()*co.b()/2+l.c()*co.d()/2,
                                     l.a()*co.b()/2+l.b()*co.c()  +l.c()*co.e()/2,
                                     l.a()*co.d()/2+l.b()*co.e()/2+l.c()*co.f()  );
    }
    else // a degenerate conic has no dual; in this case, return the centre:
        if (a_==0 && b_==0 && d_==0) // two horizontal lines
            return vgl_homg_point_2d<T>(1,0,0);
        else if (a_*c_*4==b_*b_ && a_*e_*2==b_*d_)
            return vgl_homg_point_2d<T>(b_*f_*2-e_*d_, d_*d_-a_*f_*4, a_*e_*2-b_*d_);
        else
            return vgl_homg_point_2d<T>(b_*e_-c_*d_*2, b_*d_-a_*e_*2, a_*c_*4-b_*b_);
}

//: Returns the curvature of the conic at point p, assuming p is on the conic.
template <class T>
double vgl_conic<T>::curvature_at(vgl_point_2d<T> const& p) const
{
    // Shorthands
    const T &a_xx = a_;
    const T &a_xy = b_;
    const T &a_yy = c_;
    const T &a_xw = d_;
    const T &a_yw = e_;
    
    const T x = p.x();
    const T y = p.y();
    
    double f_x  = 2*a_xx*x + a_xy*y + a_xw;
    double f_y  = 2*a_yy*y + a_xy*x + a_yw;
    double f_xy = a_xy;
    double f_xx = 2*a_xx;
    double f_yy = 2*a_yy;
    
    double f_x_2 = f_x*f_x;
    double f_y_2 = f_y*f_y;
    double denom = f_x_2 + f_y_2;
    denom = std::sqrt(denom*denom*denom);
    
    // Divergent of the unit normal grad f/|grad f|
    return (f_xx*f_y_2 - 2*f_x*f_y*f_xy + f_yy*f_x_2) / denom;
}


//: Write "<vgl_conic aX^2+bXY+cY^2+dXW+eYW+fW^2=0>" to stream
template <class T>
std::ostream& operator<<(std::ostream& s, vgl_conic<T> const& co)
{
    s << "<vgl_conic ";
    if (co.a() == 1) s << "X^2";
    else if (co.a() == -1) s << "-X^2";
    else if (co.a() != 0) s << co.a() << "X^2";
    if (co.b() > 0) s << '+';
    if (co.b() == 1) s << "XY";
    else if (co.b() == -1) s << "-XY";
    else if (co.b() != 0) s << co.b() << "XY";
    if (co.c() > 0) s << '+';
    if (co.c() == 1) s << "Y^2";
    else if (co.c() == -1) s << "-Y^2";
    else if (co.c() != 0) s << co.c() << "Y^2";
    if (co.d() > 0) s << '+';
    if (co.d() == 1) s << "XW";
    else if (co.d() == -1) s << "-XW";
    else if (co.d() != 0) s << co.d() << "XW";
    if (co.e() > 0) s << '+';
    if (co.e() == 1) s << "YW";
    else if (co.e() == -1) s << "-YW";
    else if (co.e() != 0) s << co.e() << "YW";
    if (co.f() > 0) s << '+';
    if (co.f() == 1) s << "W^2";
    else if (co.f() == -1) s << "-W^2";
    else if (co.f() != 0) s << co.f() << "W^2";
    return s << "=0 " << co.real_type() << "> ";
}


#endif // vgl_conic_h_
