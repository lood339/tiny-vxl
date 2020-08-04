// This is core/vgl/algo/vgl_h_matrix_2d_compute_4point.h
#ifndef vgl_h_matrix_2d_compute_4point_h_
#define vgl_h_matrix_2d_compute_4point_h_
//--------------------------------------------------------------
//:
// \file
//
// vgl_h_matrix_2d_compute_linear contains a linear method to calculate
// the plane projectivity which relates four 2D point correspondences.
// The returned $H$ is such that
// $H ~ [p_1 ~ p_2 ~ p_3 ~ p_4 ] \sim [p'_1 ~ p'_2 ~ p'_3 ~ p'_4 ]$
// where the $p_i$ are the homogeneous points in the first view, and the
// $p'_i$ their images.
//
// \verbatim
//  Modifications
//   08-02-98 FSM obsoleted bool compute(vgl_h_matrix_2d<double>  *)
//   Mar 26, 2003 JLM Preparing to move up to vgl
//   Jun 23, 2003 Peter Vanroose - made compute_pl() etc. pure virtual
//   Jun 23, 2003 Peter Vanroose - added rough first impl. for compute_l()
// \endverbatim

#include <cassert>
#include <vgl/algo/vgl_h_matrix_2d_compute.h>
#include <vnl/vnl_inverse.h>

class vgl_h_matrix_2d_compute_4point : public vgl_h_matrix_2d_compute
{
 public:
  int minimum_number_of_correspondences() const override { return 4; }

 protected:
  //: compute from matched points

  inline bool compute_p(std::vector<vgl_homg_point_2d<double> > const& points1,
                 std::vector<vgl_homg_point_2d<double> > const& points2,
                 vgl_h_matrix_2d<double>& H) override;

  //:compute from matched lines

  inline bool compute_l(std::vector<vgl_homg_line_2d<double> > const& lines1,
                 std::vector<vgl_homg_line_2d<double> > const& lines2,
                 vgl_h_matrix_2d<double>& H) override;

  //:compute from matched lines with weight vector

  inline bool compute_l(std::vector<vgl_homg_line_2d<double> > const& lines1,
                 std::vector<vgl_homg_line_2d<double> > const& lines2,
                 std::vector<double> const& weights,
                 vgl_h_matrix_2d<double>& H) override;

  //:compute from matched points and lines

  inline bool compute_pl(std::vector<vgl_homg_point_2d<double> > const& points1,
                  std::vector<vgl_homg_point_2d<double> > const& points2,
                  std::vector<vgl_homg_line_2d<double> > const& lines1,
                  std::vector<vgl_homg_line_2d<double> > const& lines2,
                  vgl_h_matrix_2d<double>& H) override;
};

// copy from .cpp
//-----------------------------------------------------------------------------
//
//: Compute a plane-plane projectivity using 4 point correspondences.
// Returns false if the calculation fails or there are fewer than four point
// matches in the list.
//
// The algorithm determines the transformation $H_i$ from each pointset to the
// canonical projective basis (see h_matrix_2d::projective_basis), and
// returns the combined transform $H = H_2^{-1} H_1$.
bool
vgl_h_matrix_2d_compute_4point::compute_p(std::vector<vgl_homg_point_2d<double>> const & points1,
                                          std::vector<vgl_homg_point_2d<double>> const & points2,
                                          vgl_h_matrix_2d<double> & H)
{
    vgl_h_matrix_2d<double> H1, H2;
    if (!H1.projective_basis(points1))
        return false;
    if (!H2.projective_basis(points2))
        return false;
    H.set(vnl_inverse(H2.get_matrix()) * H1.get_matrix());
    return true;
}

//-----------------------------------------------------------------------------
//
//: Compute a plane-plane projectivity using 4 line correspondences.
// Returns false if the calculation fails or there are fewer than four line
// matches in the list.
//
// Implementation is the dual of the implementation of compute_p()
bool
vgl_h_matrix_2d_compute_4point::compute_l(std::vector<vgl_homg_line_2d<double>> const & lines1,
                                          std::vector<vgl_homg_line_2d<double>> const & lines2,
                                          vgl_h_matrix_2d<double> & H)
{
    vgl_h_matrix_2d<double> H1, H2;
    if (!H1.projective_basis(lines1))
        return false;
    if (!H2.projective_basis(lines2))
        return false;
    H.set(vnl_inverse(H2.get_matrix()) * H1.get_matrix());
    return true;
}

bool
vgl_h_matrix_2d_compute_4point::compute_l(std::vector<vgl_homg_line_2d<double>> const & lines1,
                                          std::vector<vgl_homg_line_2d<double>> const & lines2,
                                          std::vector<double> const &,
                                          vgl_h_matrix_2d<double> & H)
{
    vgl_h_matrix_2d<double> H1, H2;
    if (!H1.projective_basis(lines1))
        return false;
    if (!H2.projective_basis(lines2))
        return false;
    H.set(vnl_inverse(H2.get_matrix()) * H1.get_matrix());
    return true;
}

bool
vgl_h_matrix_2d_compute_4point::compute_pl(std::vector<vgl_homg_point_2d<double>> const & /*points1*/,
                                           std::vector<vgl_homg_point_2d<double>> const & /*points2*/,
                                           std::vector<vgl_homg_line_2d<double>> const & /*lines1*/,
                                           std::vector<vgl_homg_line_2d<double>> const & /*lines2*/,
                                           vgl_h_matrix_2d<double> &)
{
    assert(!"vgl_h_matrix_2d_compute_4point::compute_pl() NYI");
    return false;
}

#endif // vgl_h_matrix_2d_compute_4point_h_
