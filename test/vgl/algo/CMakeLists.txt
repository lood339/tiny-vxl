# This is core/vgl/tests/CMakeLists.txt

project(vgl_algo_test_all)

enable_testing()

include_directories(${VXL_SRC_DIR})
include_directories(${Eigen_SRC_DIR})

add_executable(vgl_algo_test_all
    test_h_matrix_2d.cpp
    test_h_matrix_3d.cpp
    test_rotation_3d.cpp
)

target_link_libraries(vgl_algo_test_all gtest gmock_main)

add_test(NAME vgl_algo_test_all COMMAND vgl_algo_test_all)



