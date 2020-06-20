// This is core/vgl/vgl_frustum_3d.h
#ifndef vgl_frustum_3d_h
#define vgl_frustum_3d_h
//:
// \file
// \brief A polygonal cone truncated by parallel planes
// \author J.L. Mundy
// \date   December 1, 2013
//
// \verbatim
//  Modifications
// None
// \endverbatim

#include <iosfwd>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>

#include "vgl_point_3d.h" // forward declare vgl datatypes
#include "vgl_frustum_3d.h"
#include "vgl_ray_3d.h"
#include "vgl_plane_3d.h"
#include "vgl_tolerance.h"
#include "vgl_intersection.h"

//: A 3D frustum is the portion of a solid (normally a cone or pyramid)
// that lies between two parallel planes cutting it.
//
template <class Type>
class vgl_frustum_3d
{
 public:
  //: default constructor
 vgl_frustum_3d(): near_plane_(0), far_plane_(0), n_top_bot_face_verts_(0){}

  //: Construct the frustum from rays
  // the corner rays intersect in a common origin point, i.e., the apex
  // of the frustum cone. It is assumed the rays are in sorted order
  // around the cone surface, so that sequential rays are coplanar
  // with a cone surface plane and the cross product r[i].dir X r[i+1].dir
  // defines an outward-pointing normal. norm is the vector
  // perpendicular to the parallel frustum faces. d0 is the distance in
  // the norm vector dirction from the apex to the closest parallel face.
  // d1 the distance in the norm direction from the apex to the far face.
  //
  vgl_frustum_3d(std::vector<vgl_ray_3d<Type> >  const& corner_rays,
                 vgl_vector_3d<Type> const& norm, Type d0, Type d1);


  //: Equality test
  inline bool operator==(vgl_frustum_3d<Type> const& other) const;

  // Data Access---------------------------------------------------------------
  const vgl_point_3d<Type>& apex() const{
    return apex_;}

  const std::vector<vgl_plane_3d<Type> >& surface_planes() const
  {return surface_planes_;}

  const vgl_plane_3d<Type>& near_plane() const
  {return surface_planes_[near_plane_];}

  const vgl_plane_3d<Type>& far_plane() const
  {return surface_planes_[far_plane_];}

  const std::vector<vgl_point_3d<Type> >& verts() const
  {return verts_;}

  const std::map<int, std::vector<int> >& faces() const
  {return faces_;}

  vgl_box_3d<Type> bounding_box() const;

  //: Get the centroid point
  vgl_point_3d<Type> centroid() const;

  //: test if the frustum is convex
  bool is_convex() const;

  // Data Control--------------------------------------------------------------

  //: Return true iff the point p is inside this frustum
  // assumes that the frustum is a convex solid
  bool contains(vgl_point_3d<Type> const& p) const;

  //: Return true if \a (x,y,z) is inside this frustum
  // assumes that the frustum is a convex solid
  bool contains(Type const& x, Type const& y, Type const& z) const;

  // I/O-----------------------------------------------------------------------


  // INTERNALS-----------------------------------------------------------------
 protected:
  // Data Members--------------------------------------------------------------
  vgl_point_3d<Type> apex_;
  vgl_vector_3d<Type> norm_;
  //: planes bounding the frustum volume including near and far plane
  std::vector<vgl_plane_3d<Type> > surface_planes_;
  int near_plane_;
  int far_plane_;
  std::vector<vgl_point_3d<Type> > verts_;
  //: the number of verts in the top or bottom face
  int n_top_bot_face_verts_;
  // key corresponds to plane index, value is clockwise verts on face boundary
  // clockwise with respect the the face normal
  std::map<int, std::vector<int> > faces_;
};

//: Write frustum to stream
// \relatesalso vgl_frustum_3d
template <class Type>
std::ostream&  operator<<(std::ostream& s, vgl_frustum_3d<Type> const& p);

// copy from .cpp
template <class Type>
vgl_frustum_3d<Type>::
vgl_frustum_3d(std::vector<vgl_ray_3d<Type> >  const& corner_rays,
               vgl_vector_3d<Type> const& norm, Type d0, Type d1){
    
    //Construct cone bounding planes. Their surface normals point outward from
    //the interior of the frustum. It is assumed the rays are in order
    //around the boundaries of the parallel faces, so that adjacent rays
    //define the cone surface planes with outward pointing normals
    int nc = static_cast<int>(corner_rays.size());
    assert(nc>=3);//need to enclose a volume
    n_top_bot_face_verts_= nc;
    norm_ = norm;
    // contruct the near plane
    vgl_ray_3d<Type> r0 = corner_rays[0];
    apex_ = r0.origin();
    // find a point on the near plane
    vgl_vector_3d<Type> dir = r0.direction();
    // ray shouldn't be perpendicular to near plane normal
    double er = std::fabs(cos_angle(norm, dir));
    assert(er > vgl_tolerance<double>::position);
    double dp = dot_product(norm, dir);
    // the normal to the near parallel face must point to the apex
    // since the apex is, by defintion, outside the frustum volume
    vgl_vector_3d<Type> snorm = norm;
    if(dp>0) snorm = -norm;
    double dp_mag = std::fabs(dp);
    double tau_near = d0/dp_mag;
    vgl_point_3d<Type> p_near = apex_ + (tau_near*dir);
    vgl_plane_3d<Type> pln_near(snorm, p_near);
    near_plane_ = static_cast<int>(surface_planes_.size());
    surface_planes_.push_back(pln_near);
    
    //define the verts on the near plane
    verts_.push_back(p_near);
    for(int i = 1; i<nc; ++i){
        // find the intersection point of ri, i>0, with the top plane
        vgl_point_3d<Type> pint;
        bool good = vgl_intersection<Type>(corner_rays[i], pln_near, pint);
        assert(good);
        verts_.push_back(pint);
    }
    // check order
    vgl_vector_3d<Type> v01 = verts_[1]-verts_[0];
    vgl_vector_3d<Type> v12 = verts_[2]-verts_[1];
    vgl_vector_3d<Type> out = cross_product(v01, v12);
    bool rays_ccw_top = dot_product(norm, out) > 0;
    if(rays_ccw_top)
        for(int i = 0; i<nc; ++i)
            faces_[near_plane_].push_back(i);
    else
        for(int i = nc-1; i>=0; --i)
            faces_[near_plane_].push_back(i);
    // for now, the test for inside assumes a convex frustum FIX_ME
    assert(this->is_convex());
    // define the far plane
    double tau_far = d1/dp_mag;
    vgl_point_3d<Type> p_far = apex_ + (tau_far*dir);
    int far_indx = static_cast<int>(verts_.size());
    verts_.push_back(p_far);
    // the normal to the far parallel face must point away from the apex
    vgl_plane_3d<Type> pln_far(-snorm, p_far);
    far_plane_ = static_cast<int>(surface_planes_.size());
    surface_planes_.push_back(pln_far);
    // compute the verts on the far plane
    for(int i = 1; i<nc; ++i){
        // find the intersection point of ri, i>0, with the top plane
        vgl_point_3d<Type> pi;
        bool good = vgl_intersection<Type>(corner_rays[i], pln_far, pi);
        assert(good);
        verts_.push_back(pi);
    }
    // form the far face.  order is reversed on far face so normal is pointing out
    if(rays_ccw_top)
        for(int i = nc-1; i>=0; --i)
            faces_[far_plane_].push_back(i+far_indx);
    else
        for(int i = 0; i<nc; ++i)
            faces_[near_plane_].push_back(i+far_indx);
    
    // find the side cone planes and associated face verts
    // each side face has four vertices, two on the near face and two
    // on the far face
    std::vector<int>& top_verts = faces_[near_plane_];
    for(int i = 0; i<nc; ++i){
        int j = top_verts[i];
        int j_next = top_verts[(i+1)%nc];
        const vgl_ray_3d<Type>& ra = corner_rays[j];
        const vgl_ray_3d<Type>& rb = corner_rays[j_next];
        vgl_plane_3d<Type> pln(ra, rb);
        int pln_index = static_cast<int>(surface_planes_.size());
        surface_planes_.push_back(pln);
        faces_[pln_index].push_back(j);  faces_[pln_index].push_back(j_next);
        faces_[pln_index].push_back(j_next+far_indx);
        faces_[pln_index].push_back(j+far_indx);
    }
}

template <class Type>
bool vgl_frustum_3d<Type>::operator==(vgl_frustum_3d<Type> const& other) const{
    // check addresses
    if(this == &other)
        return true;
    // check apex
    if(!(apex_ == other.apex()))
        return false;
    int n = static_cast<int>(verts_.size());
    const std::vector<vgl_point_3d<Type> >& o_verts = other.verts();
    // could be round off error
    for(int i = 0; i<n; ++i){
        double dif = (verts_[i] - o_verts[i]).length();
        if( dif > vgl_tolerance<double>::position)
            return false;
    }
    return true;
}
template <class Type>
vgl_box_3d<Type> vgl_frustum_3d<Type>::bounding_box() const{
    vgl_box_3d<Type> box;
    int n = static_cast<int>(verts_.size());
    for(int i = 0; i<n; ++i)
        box.add(verts_[i]);
    return box;
}
template <class Type>
vgl_point_3d<Type> vgl_frustum_3d<Type>::centroid() const{
    int n = static_cast<int>(verts_.size());
    double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
    for(int i = 0; i<n; ++i){
        x_sum += static_cast<double>(verts_[i].x());
        y_sum += static_cast<double>(verts_[i].y());
        z_sum += static_cast<double>(verts_[i].z());
    }
    x_sum /= n; y_sum /= n; z_sum /= n;
    Type x = static_cast<Type>(x_sum);
    Type y = static_cast<Type>(y_sum);
    Type z = static_cast<Type>(z_sum);
    return vgl_point_3d<Type>(x, y, z);
}
// traverse the top face and compute the cross product of sequential
// face edges
template <class Type>
bool vgl_frustum_3d<Type>::is_convex() const{
    int n = n_top_bot_face_verts_;
    if(n<3) return false;
    if(n==3) return true;
    // n > 3
    // get the vertex map for the top face
    std::map<int, std::vector<int> >::const_iterator vit;
    vit = faces_.find(near_plane_);//use find since [] is non-const
    if(vit==faces_.end()) return false;
    const std::vector<int>& vindx = (*vit).second;
    vgl_vector_3d<Type> v = verts_[vindx[1]]-verts_[vindx[0]];
    vgl_vector_3d<Type> v_pre = verts_[vindx[2]]-verts_[vindx[1]];
    vgl_vector_3d<Type> cr = cross_product(v, v_pre);
    Type dp = dot_product(cr, norm_);
    bool pos = dp > vgl_tolerance<Type>::position;
    for(int i = 2; i<n; ++i){
        int j = (vindx[i]+1)%n;
        vgl_vector_3d<Type> v_nxt = verts_[j]-verts_[vindx[i]];
        cr = cross_product(v_pre, v_nxt);
        dp = dot_product(cr, norm_);
        bool pos_i = dp > vgl_tolerance<Type>::position;
        if(pos_i != pos) return false;
        v_pre = v_nxt;
    }
    return true;
}

template <class Type>
bool vgl_frustum_3d<Type>::
contains(Type const& x, Type const& y, Type const& z) const{
    
    // the point must be on the inside of all the faces,
    // assuming that the fustrum is a convex solid.
    int n = static_cast<int>(surface_planes_.size());
    bool inside = true;
    for(int i = 0; i<n&&inside; ++i){
        Type a = surface_planes_[i].a(), b = surface_planes_[i].b();
        Type c = surface_planes_[i].c(), d = surface_planes_[i].d();
        Type sign = a*x + b*y + c*z +d;
        inside = sign < vgl_tolerance<Type>::position;
    }
    return inside;
}

template <class Type>
bool vgl_frustum_3d<Type>::contains(vgl_point_3d<Type> const& p) const{
    return contains(p.x(), p.y(), p.z());
}

template <class Type>
std::ostream&  operator<<(std::ostream& s, vgl_frustum_3d<Type> const& f){
    s << "<vgl_frustum_3d [\n";
    const std::vector<vgl_point_3d<Type> >& verts = f.verts();
    int n = static_cast<int>(verts.size());
    for(int i = 0; i<n; ++i)
        s << verts[i] << '\n';
    s << "] >\n";
    return s;
}



#endif // vgl_frustum_3d_h
