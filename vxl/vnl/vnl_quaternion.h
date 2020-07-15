// This is core/vnl/vnl_quaternion.h
#ifndef vnl_quaternion_h_
#define vnl_quaternion_h_
//:
// \file
// \brief Unit quaternion represents rotation in 3D.
// \author awf@robots.ox.ac.uk
// \date   16 Mar 00
//
// \verbatim
//  Modifications
//   20-05-2000 fsm. changed FLOAT to T since gcc will barf at
//              the very reasonable forward declaration
//              template <class FLOAT> class vnl_quaternion;
//   23-3-2001 LSB (Manchester) Tidied documentation
//   13-1-2003 Peter Vanroose - removed unimplemented method rotation_matrix()
//   20-2-2006 Ian Scott - Added conversion to from Euler angles
//   06-5-2006 Peter Vanroose - replaced all vnl_vector by vnl_vector_fixed
// \endverbatim

#include <iostream>
#include <cmath>
#include <limits>

#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_export.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_math.h>

//: 4-element vector that represents rotation in 3D.
// vnl_quaternion is a 4-element vector with 1 real and 3 imaginary
// components:
// \code
//    q = r + (i*x + j*y + k*z)
//    r = cos(theta/2)
//    (x, y, z) = sin(theta/2) (kx, ky, kz)
// \endcode
// where theta and k are respectively the angle and axis of rotation.
//
// 3D vectors can be thought of as pure imaginary quaternions, and so a
// quaternion is represented as a vnl_vector_fixed<T,4> with the imaginary
// part before the real part for 1-1 alignment.
//
// Unit quaternions (i.e., for which $x^2 + y^2 + z^2 + r^2 = 1$)
// provide a more efficient representation for rotation
// than the usual orthonormal matrix that has nine
// parameters and six orthonormal constraints.  The unit
// quaternion has only one unit magnitude constraint.  Composing
// rotations with quaternions results in fewer multiplications
// and less error.  To insure valid rotation results, the
// nearest unit quaternion is computed, and this is much easier
// than finding the nearest orthonormal matrix.  Transforming
// vectors with a quaternion requires more operations compared
// to multiplication with the equivalent orthonormal matrix.
//
// \sa
// vnl_vector_fixed and vnl_matrix_fixed for basic operations on vectors and matrices.
// \sa
// Envelope for envelope-letter scheme that avoids deep copy on
// return by value in arithmetic expressions like: q1 * q2 * q3 *...
//

template <class T>
class VNL_EXPORT vnl_quaternion : public vnl_vector_fixed<T, 4>
{
 private:
  using Base = vnl_vector_fixed<T,4>;
 public:

  //: Constructor for null quaternion
  vnl_quaternion() = default;

  //: Construct quaternion from components x,y,z,r
  vnl_quaternion(T x, T y, T z, T r);

  //: Construct quaternion from Euler Angles,
  // That is a rotation about the X axis, followed by Y, followed by
  // the Z axis, using a fixed reference frame.
  vnl_quaternion(T theta_X, T theta_Y, T theta_Z);

  //: Construct quaternion from axis and angle of rotation.
  // \note If you specify an angle in [0, 2pi], then methods angle() and axis() will return the same values.
  // However, if you specify an angle in [-2pi, 0], then methods angle() and axis() will return values with opposite signs.
  // \sa vnl_quaternion::angle()
  // \sa vnl_quaternion::axis()
  vnl_quaternion(vnl_vector_fixed<T,3> const& axis, double angle);

  //: Construct quaternion from 3x3 row-major matrix
  explicit vnl_quaternion(vnl_matrix_fixed<T,3,3> const& transform);

  //: Construct quaternion from a 3D vector
  vnl_quaternion(vnl_vector_fixed<T,3> const& vec);

  //: Construct quaternion from a 4D vector
  vnl_quaternion (vnl_vector_fixed<T,4> const& vec);

  //: Copy constructor -- Creates a copy of from quaternion.
  inline vnl_quaternion(vnl_quaternion<T> const& from) : Base(from) {}

  //: Free internal array
  inline ~vnl_quaternion() = default; // vnl_vector_fixed will free data array

  //:  Overloads assignment operator to copy rhs quaternion into lhs quaternion.
  inline vnl_quaternion& operator= (vnl_quaternion<T> const& rhs) { Base::operator=(rhs); return *this; }

  //: Imaginary component, parallel to axis of rotation.
  // Use this accessor to both get and set the component.
  inline T& x() { return this->operator()(0); }
  //: Imaginary component, parallel to axis of rotation.
  // Use this accessor to both get and set the component.
  inline T& y() { return this->operator()(1); }
  //: Imaginary component, parallel to axis of rotation.
  // Use this accessor to both get and set the component.
  inline T& z() { return this->operator()(2); }
  //: Real component.
  // Use this accessor to both get and set the component.
  inline T& r() { return this->operator()(3); }

  //: Imaginary component, parallel to axis of rotation.
  // Use this accessor to get the component.
  inline T x() const { return this->operator()(0); }
  //: Imaginary component, parallel to axis of rotation.
  // Use this accessor to get the component.
  inline T y() const { return this->operator()(1); }
  //: Imaginary component, parallel to axis of rotation.
  // Use this accessor to get the component.
  inline T z() const { return this->operator()(2); }
  //: Real component.
  // Use this accessor to get the component.
  inline T r() const { return this->operator()(3); }

  //: Copies and returns the real part.
  inline T real() const { return (*this)[3]; }

  //: Copies and returns the imaginary part.
  inline vnl_vector_fixed<T,3> imaginary() const { return this->extract(3,0); }

  //: Axis of rotation.
  // \note Axis not well defined for theta==0. In such a case (or if provided axis==(0,0,0)), this function returns (0,0,1).
  vnl_vector_fixed<T,3> axis() const;

  //: Angle of rotation.
  // \note Returned angle lies in [0, 2*pi]
  double angle() const;

  //: 3x3 rotation matrix
  // The orthonormal vectors are the rows of the matrix, not its columns
  vnl_matrix_fixed<T,3,3> rotation_matrix_transpose() const;

  //: 4x4 rotation matrix
  vnl_matrix_fixed<T,4,4> rotation_matrix_transpose_4() const;

  //: Same real, opposite img part
  vnl_quaternion<T> conjugate() const;

  //: Inverse for nonzero quat
  vnl_quaternion<T> inverse() const;

  vnl_quaternion<T> operator* (vnl_quaternion<T> const&) const;

  //: Rotate 3D v
  // The quaternion must be normalised first.
  vnl_vector_fixed<T,3> rotate(vnl_vector_fixed<T,3> const&) const;

  //: Rotation representation in Euler angles.
  // The angles returned will be [theta_X,theta_Y,theta_Z]
  // where the final rotation is found be first applying theta_X radians
  // about the X axis, then theta_Y about the Y-axis, etc.
  // The axes stay in a fixed reference frame.
  // The quaternion mut be normalised first.
  vnl_vector_fixed<T,3> rotation_euler_angles() const;
};

//: operator<<
// \relatesalso vnl_quaternion
template <class T>
inline std::ostream& operator<< (std::ostream& os, vnl_quaternion<T> const& q)
{
  return os << *((const vnl_vector_fixed<T,4>*) &q);
}

// copy from .cpp
//: Creates a quaternion from its ordered components.
// x, y, z denote the imaginary part, which are the  coordinates
// of the rotation axis multiplied by the sine of half the
// angle of rotation. r denotes  the  real  part,  or  the
// cosine  of  half  the  angle of rotation. Default is to
// create a null quaternion, corresponding to a null rotation
// or  an  identity  transform,  which has undefined
// rotation axis.

template <class T>
vnl_quaternion<T>::vnl_quaternion (T tx, T ty, T tz, T rea)
{
    this->operator[](0) = tx;  // 3 first elements are
    this->operator[](1) = ty;  // imaginary parts
    this->operator[](2) = tz;
    this->operator[](3) = rea;  // last element is real part
}

//: Creates a quaternion from the normalized axis direction and the angle of rotation in radians.

template <class T>
vnl_quaternion<T>::vnl_quaternion(vnl_vector_fixed<T,3> const& Axis, double Angle)
{
    double a = Angle * 0.5;  // half angle
    T s = T(std::sin(a));
    for (int i = 0; i < 3; i++)            // imaginary vector is sine of
        this->operator[](i) = T(s * Axis(i));// half angle multiplied with axis
    this->operator[](3) = T(std::cos(a));   // real part is cosine of half angle
}

//: Creates a quaternion from a vector.
// 3D vector is converted into an imaginary quaternion with same
// (x, y, z) components.

template <class T>
vnl_quaternion<T>::vnl_quaternion(vnl_vector_fixed<T,3> const& vec)
{
    for (unsigned int i = 0; i < 3; ++i)
        this->operator[](i) = vec(i);
    this->operator[](3) = T(0);
}

//: Creates a quaternion from a vector.
// 4D vector is assumed to be a 4-element quaternion, to
// provide casting between vector and quaternion

template <class T>
vnl_quaternion<T>::vnl_quaternion(vnl_vector_fixed<T,4> const& vec)
{
    for (unsigned int i = 0; i < 4; ++i) // 1-1 layout between vector & quaternion
        this->operator[](i) = vec[i];
}


//: Creates a quaternion from a rotation matrix.
// Its orthonormal basis vectors are the matrix rows.
// \note this matrix \e must have determinant +1; this is not verified!
// \warning Takes the transpose of the rotation matrix, i.e.,
// the orthonormal vectors must be the rows of the matrix, not the columns.
template <class T>
vnl_quaternion<T>::vnl_quaternion(vnl_matrix_fixed<T,3,3> const& rot)
{
    double d0 = rot(0,0), d1 = rot(1,1), d2 = rot(2,2);
    double xx = 1.0 + d0 - d1 - d2;      // from the diagonal of the rotation
    double yy = 1.0 - d0 + d1 - d2;      // matrix, find the terms in
    double zz = 1.0 - d0 - d1 + d2;      // each Quaternion component
    double rr = 1.0 + d0 + d1 + d2;      // (using the fact that rr+xx+yy+zz=4)
    
    double max = rr;                     // find the maximum of all terms;
    if (xx > max) max = xx;              // dividing by the maximum makes
    if (yy > max) max = yy;              // the computations more stable
    if (zz > max) max = zz;              // and avoid division by zero
    
    if (rr == max) {
        T r4 = T(std::sqrt(rr)*2);
        this->r() = r4 / 4;
        r4 = T(1) / r4;
        this->x() = (rot(1,2) - rot(2,1)) * r4;     // find other components from
        this->y() = (rot(2,0) - rot(0,2)) * r4;     // off diagonal terms of
        this->z() = (rot(0,1) - rot(1,0)) * r4;     // rotation matrix.
    }
    else if (xx == max) {
        T x4 = T(std::sqrt(xx)*2);
        this->x() = x4 / 4;
        x4 = T(1) / x4;
        this->y() = (rot(0,1) + rot(1,0)) * x4;
        this->z() = (rot(0,2) + rot(2,0)) * x4;
        this->r() = (rot(1,2) - rot(2,1)) * x4;
    }
    else if (yy == max) {
        T y4 = T(std::sqrt(yy)*2);
        this->y() =  y4 / 4;
        y4 = T(1) / y4;
        this->x() = (rot(0,1) + rot(1,0)) * y4;
        this->z() = (rot(1,2) + rot(2,1)) * y4;
        this->r() = (rot(2,0) - rot(0,2)) * y4;
    }
    else {
        T z4 = T(std::sqrt(zz)*2);
        this->z() =  z4 / 4;
        z4 = T(1) / z4;
        this->x() = (rot(0,2) + rot(2,0)) * z4;
        this->y() = (rot(1,2) + rot(2,1)) * z4;
        this->r() = (rot(0,1) - rot(1,0)) * z4;
    }
}


//: Construct quaternion from Euler Angles
// That is a rotation about the X axis, followed by Y, followed by
// the Z axis, using a fixed reference frame.
template <class T>
vnl_quaternion<T>::vnl_quaternion(T theta_X, T theta_Y, T theta_Z)
{
    vnl_quaternion<T> Rx(T(std::sin(double(theta_X)*0.5)), 0, 0, T(std::cos(double(theta_X)*0.5)));
    vnl_quaternion<T> Ry(0, T(std::sin(double(theta_Y)*0.5)), 0, T(std::cos(double(theta_Y)*0.5)));
    vnl_quaternion<T> Rz(0, 0, T(std::sin(double(theta_Z)*0.5)), T(std::cos(double(theta_Z)*0.5)));
    *this = Rz * Ry * Rx;
}

//: Rotation representation in Euler angles.
// The angles returned will be [theta_X,theta_Y,theta_Z]
// where the final rotation is found be first applying theta_X radians
// about the X axis, then theta_Y about the Y-axis, etc.
// The axes stay in a fixed reference frame.
template <class T>
vnl_vector_fixed<T,3> vnl_quaternion<T>::rotation_euler_angles() const
{
    vnl_vector_fixed<T,3> angles;
    
    vnl_matrix_fixed<T,4,4> rotM = rotation_matrix_transpose_4();
    T xy = T(std::sqrt(double(vnl_math::sqr(rotM(0,0)) + vnl_math::sqr(rotM(0,1)))));
    if (xy > std::numeric_limits<T>::epsilon() * T(8))
    {
        angles(0) = T(std::atan2(double(rotM(1,2)), double(rotM(2,2))));
        angles(1) = T(std::atan2(double(-rotM(0,2)), double(xy)));
        angles(2) = T(std::atan2(double(rotM(0,1)), double(rotM(0,0))));
    }
    else
    {
        angles(0) = T(std::atan2(double(-rotM(2,1)), double(rotM(1,1))));
        angles(1) = T(std::atan2(double(-rotM(0,2)), double(xy)));
        angles(2) = T(0);
    }
    return angles;
}


//: Queries the rotation angle of the quaternion.
//  Returned angle lies in [0, 2*pi]
template <class T>
double vnl_quaternion<T>::angle() const
{
    return 2 * std::atan2(double(this->imaginary().magnitude()),
                          double(this->real()));    // angle is always positive
}

//: Queries the direction of the rotation axis of the quaternion.
//  A null quaternion will return zero for angle and k direction for axis.
template <class T>
vnl_vector_fixed<T,3> vnl_quaternion<T>::axis() const
{
    vnl_vector_fixed<T,3> direc = this->imaginary(); // direc parallel to imag. part
    T mag = direc.magnitude();
    if (mag == T(0)) {
        std::cout << "Axis not well defined for zero Quaternion. Using (0,0,1) instead.\n";
        direc[2] = T(1);                    // or signal exception here.
    }
    else
        direc /= mag;                       // normalize direction vector
    return direc;
}


//: Converts a normalized quaternion into a square rotation matrix with dimension dim.
//  This is the reverse counterpart of constructing a quaternion from a transformation matrix.
// WARNING this is inconsistent with the quaternion docs and q.rotate()

template <class T>
vnl_matrix_fixed<T,3,3> vnl_quaternion<T>::rotation_matrix_transpose() const
{
    T x2 = x() * x(),  xy = x() * y(),  rx = r() * x(),
    y2 = y() * y(),  yz = y() * z(),  ry = r() * y(),
    z2 = z() * z(),  zx = z() * x(),  rz = r() * z(),
    r2 = r() * r();
    vnl_matrix_fixed<T,3,3> rot;
    rot(0,0) = r2 + x2 - y2 - z2;         // fill diagonal terms
    rot(1,1) = r2 - x2 + y2 - z2;
    rot(2,2) = r2 - x2 - y2 + z2;
    rot(0,1) = 2 * (xy + rz);             // fill off diagonal terms
    rot(0,2) = 2 * (zx - ry);
    rot(1,2) = 2 * (yz + rx);
    rot(1,0) = 2 * (xy - rz);
    rot(2,0) = 2 * (zx + ry);
    rot(2,1) = 2 * (yz - rx);
    
    return rot;
}


template <class T>
vnl_matrix_fixed<T,4,4> vnl_quaternion<T>::rotation_matrix_transpose_4() const
{
    vnl_matrix_fixed<T,4,4> rot;
    return rot.set_identity().update(this->rotation_matrix_transpose().as_matrix());
}

//: Returns the conjugate of given quaternion, having same real and opposite imaginary parts.

template <class T>
vnl_quaternion<T> vnl_quaternion<T>::conjugate() const
{
    return vnl_quaternion<T> (-x(), -y(), -z(), r());
}

//: Returns the inverse of given quaternion.
//  For unit quaternion representing rotation, the inverse is the
// same as the conjugate.

template <class T>
vnl_quaternion<T> vnl_quaternion<T>::inverse() const
{
    vnl_quaternion<T> inv = this->conjugate();
    //inv /= vnl_c_vector<T>::dot_product(this->data_, this->data_, 4);
    inv /= doc_product(Base::as_vector(), Base::as_vector());
    return inv;
}

//: Returns  the product of two quaternions.
// Multiplication of two quaternions is not symmetric and has
// fewer  operations  than  multiplication  of orthonormal
// matrices. If object is rotated by r1, then by r2,  then
// the  composed  rotation (r2 o r1) is represented by the
// quaternion (q2 * q1), or by the matrix (m1 * m2).  Note
// that  matrix  composition  is reversed because matrices
// and vectors are represented row-wise.

template <class T>
vnl_quaternion<T> vnl_quaternion<T>::operator* (vnl_quaternion<T> const& rhs) const
{
    T r1 = this->real();                  // real and img parts of args
    T r2 = rhs.real();
    vnl_vector_fixed<T,3> i1 = this->imaginary();
    vnl_vector_fixed<T,3> i2 = rhs.imaginary();
    T real_v = (r1 * r2) - ::dot_product(i1, i2); // real&img of product q1*q2
    vnl_vector_fixed<T,3> img = vnl_cross_3d(i1, i2);
    img += (i2 * r1) + (i1 * r2);
    return vnl_quaternion<T>(img[0], img[1], img[2], real_v);
}

//: Rotates 3D vector v with source quaternion and stores the rotated vector back into v.
// For speed and greater accuracy, first convert quaternion into an orthonormal
// matrix,  then  use matrix multiplication to rotate many vectors.

template <class T>
vnl_vector_fixed<T,3> vnl_quaternion<T>::rotate(vnl_vector_fixed<T,3> const& v) const
{
    T rea = this->real();
    vnl_vector_fixed<T,3> i = this->imaginary();
    vnl_vector_fixed<T,3> i_x_v(vnl_cross_3d(i, v));
    return v + i_x_v * T(2*rea) - vnl_cross_3d(i_x_v, i) * T(2);
}

#endif // vnl_quaternion_h_
