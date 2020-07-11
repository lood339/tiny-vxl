# This is core/vnl/tests/CMakeLists.txt

project(vnl_test_all)

enable_testing()

include_directories(${VXL_SRC_DIR})
include_directories(${Eigen_SRC_DIR})

add_executable(vnl_test_all
  test_diag_matrix.cpp 
  test_matrix.cpp
  test_matrix_fixed.cpp
  test_vector.cpp
  test_vector_fixed.cpp
)

target_link_libraries(vnl_test_all gtest gmock_main)
add_test(NAME vnl_test_all COMMAND vnl_test_all)

add_subdirectory(algo)




