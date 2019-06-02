add_library(matrix SHARED "")

set(USE_DEBUG_MODEL ON)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

target_compile_options(matrix PUBLIC -Wall -Wextra -Wpedantic -Werror)

if (APPLE)
  message("OSX detected. Link cblas...")
  #target_link_libraries(matrix INTERFACE "-framework Foundation")

  set(BLAS_INCLUDE_DIR
    "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers"
  )

  #add_compile_definitions(WITH_VECLIB_CBLAS)
  #unset(matrix_LIB_DEPENDS CACHE)
else()
  find_package(BLAS)
  if(BLAS_FOUND)
    message("Found Blas located at:")
    list(GET BLAS_LIBRARIES 0 BLAS_LIBRARY_DIR)
    message("${BLAS_LIBRARY_DIR}")
    find_path(BLAS_INCLUDE_DIRS cblas.h
      ${BLAS_LIBRARY_DIR}/include)
    list(GET BLAS_INCLUDE_DIRS 0 BLAS_INCLUDE_DIR)
    message("Blas Include Dir: ${BLAS_INCLUDE_DIR}")
  else()
    message("Use OpenBLAS-0.3.6")
    configure_file("${PROJECT_SOURCE_DIR}/dependencies/OpenBLAS-0.3.6/cmake/OpenBLASConfig.cmake.in" "${PROJECT_BINARY_DIR}/OpenBLASConfig.cmake")
    include(${PROJECT_BINARY_DIR}/OpenBLASConfig.cmake)
    find_package(BLAS REQUIRED)
    if(BLAS_FOUND)
      message("Found Blas located at:")
      list(GET BLAS_LIBRARIES 0 BLAS_LIBRARY_DIR)
      message("${BLAS_LIBRARY_DIR}")
      find_path(BLAS_INCLUDE_DIRS cblas.h
        ${BLAS_LIBRARY_DIR}/include)
      list(GET BLAS_INCLUDE_DIRS 0 BLAS_INCLUDE_DIR)
      message("Blas Include Dir: ${BLAS_INCLUDE_DIR}")
    endif()
  endif()
endif()

target_link_libraries(matrix PRIVATE cblas clapack)
target_include_directories(matrix
  PUBLIC
  "${BLAS_INCLUDE_DIR}"
  "${PROJECT_SOURCE_DIR}/matrix_impl/include"
)

include("include/CMakeLists.txt")
include("src/CMakeLists.txt")

