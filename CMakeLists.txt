cmake_minimum_required (VERSION 3.14)
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()
project (matrix CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_TESTS "build unittests as well" ON)

# Add sub-directories
add_subdirectory(matrix_impl)

include(cmake/googletest.cmake)
fetch_googletest(
  ${PROJECT_SOURCE_DIR}/cmake
  ${PROJECT_BINARY_DIR}/googletest
  )

if(BUILD_TESTS)
  enable_testing()
  
  add_subdirectory(tests)
endif()

