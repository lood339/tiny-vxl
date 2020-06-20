#ifndef VGL_H
#define VGL_H

// General
#include "vgl/vgl_fwd.h"
#include "vgl/vgl_tolerance.h"
#include "vgl/vgl_box_2d.h"
#include "vgl/vgl_box_3d.h"
#include "vgl/vgl_vector_2d.h"
#include "vgl/vgl_vector_3d.h"
#include "vgl/vgl_1d_basis.h"
#include "vgl/vgl_homg.h"

// Points
#include "vgl/vgl_homg_point_1d.h"
#include "vgl/vgl_point_2d.h"
#include "vgl/vgl_homg_point_2d.h"
#include "vgl/vgl_point_3d.h"
#include "vgl/vgl_homg_point_3d.h"
#include "vgl/vgl_pointset_3d.h"

// Lines
#include "vgl/vgl_line_2d.h"
#include "vgl/vgl_homg_line_2d.h"
#include "vgl/vgl_homg_line_3d_2_points.h"
#include "vgl/vgl_line_3d_2_points.h"
#include "vgl/vgl_line_segment_2d.h"
#include "vgl/vgl_line_segment_3d.h"
#include "vgl/vgl_infinite_line_3d.h"
#include "vgl/vgl_ray_3d.h"

// Other curves
#include "vgl/vgl_conic.h"
#include "vgl/vgl_conic_segment_2d.h"
#include "vgl/vgl_polygon.h"
#include "vgl/vgl_sphere_3d.h"
#include "vgl/vgl_cylinder.h"

// Planes
#include "vgl/vgl_plane_3d.h"
#include "vgl/vgl_homg_plane_3d.h"

// Spline
#include "vgl/vgl_cubic_spline_2d.h"
#include "vgl/vgl_cubic_spline_3d.h"

// Functions
#include "vgl/vgl_closest_point.h"
#include "vgl/vgl_distance.h"
#include "vgl/vgl_clip.h"
#include "vgl/vgl_area.h"
#include "vgl/vgl_convex.h"
#include "vgl/vgl_intersection.h"
#include "vgl/vgl_bounding_box.h"
#include "vgl/vgl_oriented_box_2d.h"
//#include "vgl/vgl_fit_oriented_box_2d.h" // moved to algo

#include "vgl/vgl_lineseg_test.h"
#include "vgl/vgl_triangle_test.h"
#include "vgl/vgl_polygon_test.h"
#include "vgl/vgl_frustum_3d.h"
#include "vgl/vgl_affine_coordinates.h"

//#include "vgl/vgl_region_scan_iterator.h"
//#include "vgl/vgl_triangle_scan_iterator.h"



using vgl_point_2d_d = vgl_point_2d<double>;
using vgl_vgl_line_segment_2d_d = vgl_line_segment_2d<double>;
using vgl_line_segment_3d_d = vgl_line_segment_3d<double>;
using vgl_box_2d_d = vgl_box_2d<double>;


template class vgl_conic<double>;
template std::ostream& operator<<(std::ostream&, const vgl_conic<double>&);

template class vgl_homg_line_2d<double>;
template std::ostream& operator<<(std::ostream&, vgl_homg_line_2d<double>const&);


template class vgl_homg_point_2d<double>;
template double cross_ratio(vgl_homg_point_2d<double >const&, vgl_homg_point_2d<double >const&, \
                            vgl_homg_point_2d<double >const&, vgl_homg_point_2d<double >const&); \
template std::ostream& operator<<(std::ostream&, vgl_homg_point_2d<double >const&);


/*
template class vgl_box_3d<Type >;\
template std::ostream& operator<<(std::ostream&, vgl_box_3d<Type > const& p);\
template std::istream& operator>>(std::istream&, vgl_box_3d<Type >& p)
*/

/*
template class vgl_homg_point_3d<double >;

template bool collinear(vgl_homg_point_3d<double >const&,
                        vgl_homg_point_3d<double >const&,
                        vgl_homg_point_3d<double >const&);

template double cross_ratio(vgl_homg_point_3d<double >const&,
                            vgl_homg_point_3d<double >const&,
                            vgl_homg_point_3d<double >const&,
                            vgl_homg_point_3d<double >const&);

template std::ostream& operator<<(std::ostream&, vgl_homg_point_3d<double >const&);
 */

/*
template class vgl_homg_plane_3d<double>;
template std::ostream& operator<<(std::ostream&, vgl_homg_plane_3d<double >const&);
 */

/*
template class vgl_plane_3d<double >;
template std::ostream& operator<<(std::ostream&, vgl_plane_3d<double >const&);
 */
#endif
