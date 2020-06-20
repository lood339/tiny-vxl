# cmake file for vxl (head file only version)

cmake_minimum_required(VERSION 3.8)


project(vxl LANGUAGES CXX) #

# c++ standard
if(NOT CMAKE_CXX_STANDARD)
   set(CMAKE_CXX_STANDARD "11")
endif()
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT CMAKE_CXX_EXTENSIONS)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
set(VALID_CXX_STANDARDS "11" "14" "17")
if(NOT CMAKE_CXX_STANDARD IN_LIST VALID_CXX_STANDARDS )
   MESSAGE(FATAL_ERROR "CMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD} not in know standards list\n ${VALID_CXX_STANDARDS}.")
endif()

# c standard
if(NOT CMAKE_C_STANDARD)
   set(CMAKE_C_STANDARD "11")
endif()
if(NOT CMAKE_C_STANDARD_REQUIRED)
  set(CMAKE_C_STANDARD_REQUIRED ON)
endif()
if(NOT CMAKE_C_EXTENSIONS)
  set(CMAKE_C_EXTENSIONS OFF)
endif()
set(VALID_C_STANDARDS "90" "99" "11")
if(NOT CMAKE_C_STANDARD IN_LIST VALID_C_STANDARDS )
   MESSAGE(FATAL_ERROR "CMAKE_C_STANDARD:STRING=${CMAKE_C_STANDARD} not in know standards list\n ${VALID_C_STANDARDS}.")
endif()


set(VXL_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/vxl)

#add_library(vxl MODULE ${VXL_SRC_DIR}/vgl.hpp)
#target_include_directories(vxl PUBLIC include)
#target_compile_features(vxl PUBLIC cxx_std_11)
#set_target_properties(vxl PROPERTIES LINKER_LANGUAGE CXX)


# test
option(BUILD_TESTING "Enable creation of Eigen tests." ON)
if(BUILD_TESTING)
  # For Google Test
  # From https://github.com/Crascit/DownloadProject
  include(CTest)

  if (CMAKE_VERSION VERSION_LESS 3.2)
      set(UPDATE_DISCONNECTED_IF_AVAILABLE "")
  else()
      set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")
  endif()

  include(cmake/DownloadProject.cmake)
  download_project(PROJ                googletest
                  GIT_REPOSITORY      https://github.com/google/googletest.git
                  GIT_TAG             master
                  ${UPDATE_DISCONNECTED_IF_AVAILABLE}
  )

  # Prevent GoogleTest from overriding our compiler/linker options
  # when building with Visual Studio
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})

  # When using CMake 2.8.11 or later, header path dependencies
  # are automatically added to the gtest and gmock targets.
  # For earlier CMake versions, we have to explicitly add the
  # required directories to the header search path ourselves.
  if (CMAKE_VERSION VERSION_LESS 2.8.11)
      include_directories("${gtest_SOURCE_DIR}/include"
                          "${gmock_SOURCE_DIR}/include")
  endif()
  
  # testing examples are vxl
  add_subdirectory(vxl)
endif()

