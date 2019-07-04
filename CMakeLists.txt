cmake_minimum_required (VERSION 3.14)
project (matrix CXX)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif(NOT DEFINED CMAKE_CXX_STANDARD)

if(NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif(NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)

if(NOT DEFINED CMAKE_CXX_COMPILER)
  if(APPLE)
    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
    message(STATUS "compiler - ${CMAKE_CXX_COMPILER}")
  endif(APPLE)
endif(NOT DEFINED CMAKE_CXX_COMPILER)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif(NOT DEFINED CMAKE_CXX_STANDARD)

option(BUILD_MATRIX_TESTS "build unittests as well" ON)

# Add sub-directories
add_subdirectory(matrix_impl)

include(cmake/googletest.cmake)
include(cmake/CodeCoverage.cmake)
fetch_googletest(
  ${PROJECT_SOURCE_DIR}/cmake
  ${PROJECT_BINARY_DIR}/googletest
  )

if(BUILD_MATRIX_TESTS)
  enable_testing()
  
  add_subdirectory(tests)
endif()

