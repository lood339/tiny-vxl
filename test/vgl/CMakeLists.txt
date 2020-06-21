# This is core/vgl/tests/CMakeLists.txt

project(vgl_test_all)

enable_testing()

include_directories(${VXL_SRC_DIR})

add_executable(vgl_test_all
  test_affine_coordinates.cpp
  test_area.cpp
  test_clip.cpp
  test_closest_point.cpp 
  test_convex.cpp
  test_frustum_3d.cpp
  test_infinite_line_3d.cpp
  test_line_segment_3d.cpp
  test_line_3d_2_points.cpp
  test_oriented_box_2d.cpp
  test_plane_3d.cpp
  test_polygon.cpp
  test_quadric.cpp
  test_ray_3d.cpp  
  test_sphere_3d.cpp
)

target_link_libraries(vgl_test_all gtest gmock_main)

add_test(NAME vgl_test_all COMMAND vgl_test_all)


