cmake_minimum_required (VERSION 3.8)
project(matrix)

add_library(matrix INTERFACE)
#  "include/Matrix.h"
#  "include/MatrixImpl.hpp"
#)

if (APPLE)
  message("OSX detected. Link cblas...")
  target_link_libraries(matrix INTERFACE cblas)
  #target_link_libraries(matrix INTERFACE "-framework Foundation")

  target_include_directories(matrix
    INTERFACE
    "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers"
    )
  #unset(matrix_LIB_DEPENDS CACHE)
endif()

target_include_directories(matrix
  INTERFACE
  "${PROJECT_SOURCE_DIR}/matrix_impl/include"
  )

include("include/CMakeLists.txt")
include("src/CMakeLists.txt")

