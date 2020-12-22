# This is core/vgl/tests/CMakeLists.txt

project(vpgl_test_all)

enable_testing()

include_directories(${VXL_SRC_DIR})
include_directories(${Eigen_SRC_DIR})

add_executable(vpgl_test_all
    test_affine_camera.cpp    
    test_calibration_matrix.cpp
    test_generic_camera.cpp
    test_perspective_camera.cpp
    test_proj_camera.cpp
)

target_link_libraries(vpgl_test_all gtest gmock_main)

add_test(NAME vpgl_test_all COMMAND vpgl_test_all)


