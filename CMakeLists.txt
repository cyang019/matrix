cmake_minimum_required (VERSION 3.9)
project (matrix 
  LANGUAGES CXX
  VERSION 0.1)

message(STATUS "Project matrix...")

if(NOT DEFINED CMAKE_CXX_STANDARD)
  message(STATUS "CPP standard - ${CMAKE_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD 17)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_CXX_COMPILER)
  if(APPLE)
    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
    message(STATUS "compiler - ${CMAKE_CXX_COMPILER}")
  endif()
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

option(BUILD_TESTS "build unittests as well" ON)

# Add sub-directories
add_subdirectory(matrix_impl)

include(cmake/googletest.cmake)
if(NDEBUG)
  include(cmake/CodeCoverage.cmake)
endif()
fetch_googletest(
  ${PROJECT_SOURCE_DIR}/cmake
  ${PROJECT_BINARY_DIR}/googletest
  )

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

option(BUILD_BENCHMARK "build benchmark as well" ON)

if(BUILD_BENCHMARK)
  add_subdirectory(benchmark_matrix)
endif()

