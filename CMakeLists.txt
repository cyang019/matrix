cmake_minimum_required (VERSION 3.8)
project (matrix CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add sub-directories
add_subdirectory(matrix_impl)

include(cmake/googletest.cmake)
fetch_googletest(
  ${PROJECT_SOURCE_DIR}/cmake
  ${PROJECT_BINARY_DIR}/googletest
  )

enable_testing()

add_subdirectory(tests)

