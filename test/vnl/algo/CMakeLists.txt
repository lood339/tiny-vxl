# This is core/vnl/algo/tests/CMakeLists.txt

project(vnl_algo_test_all)

enable_testing()

include_directories(${VXL_SRC_DIR})
include_directories(${Eigen_SRC_DIR})

add_executable(vnl_algo_test_all
  test_svd.cpp
)

target_link_libraries(vnl_algo_test_all gtest gmock_main)
add_test(NAME vnl_algo_test_all COMMAND vnl_algo_test_all)





