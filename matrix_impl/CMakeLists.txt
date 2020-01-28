add_library(matrix SHARED "")

if(NOT DEFINED USE_DEBUG_MODEL) 
  set(USE_DEBUG_MODEL ON)
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_DEBUG)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_RELEASE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
endif()

target_compile_options(matrix PUBLIC -Wall -Wextra -Wpedantic -Werror)

set(_possible_paths ${CMAKE_PREFIX_PATH}/include /usr/include /usr/local/include)

if (APPLE)
  message( STATUS "OSX detected. Link cblas..." )
  #target_link_libraries(matrix INTERFACE "-framework Foundation")
  find_package(BLAS)
  if(BLAS_FOUND)
    message( STATUS "BLAS found: ${BLAS_LIBRARIES}" )
    message( STATUS BLAS_HOME: ${BLAS_HOME} )
    message( STATUS BLAS_DIR: ${BLAS_DIR} )
    message( STATUS BLAS_LIBDIR: ${BLAS_LIBDIR} )
    list(GET BLAS_LIBRARIES 0 BLAS_LIBRARY_DIR)
    find_path(BLAS_INCLUDE_DIRS 
      NAMES 
        cblas.h
      HINTS
        /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
        ${_possible_paths}
        $ENV{BLAS_HOME}/include)
    list(GET BLAS_INCLUDE_DIRS 0 BLAS_INCLUDE_DIR)
    INCLUDE_DIRECTORIES(${BLAS_INCLUDE_DIR})
    message( STATUS "cblas.h => ${BLAS_INCLUDE_DIR}" )
  else()
    message( FATAL_ERROR "Missing BLAS.")
  endif()

  # ===== Look For Lapack =====
  find_package(LAPACK)
  if(LAPACK_FOUND)
    message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
    message( STATUS LAPACK_HOME: ${LAPACK_HOME} )
    message( STATUS LAPACK_DIR: ${LAPACK_DIR} )
    message( STATUS LAPACK_LIBDIR: ${LAPACK_LIBDIR} )
    find_path(CLAPACK_INCLUDE_DIRS 
      NAMES 
        clapack.h
      PATHS
        /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
        ${_possible_paths}
        $ENV{LAPACK_HOME}/include)
  else()
    message( FATAL_ERROR "LAPACK c interface headers missing.")
  endif()
  #set(BLAS_INCLUDE_DIR
  #  "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers"
  #)
  option(HAVE_APPLE_LAPACK "APPLE_LAPACK" ON)

  #add_compile_definitions(WITH_VECLIB_CBLAS)
  #unset(matrix_LIB_DEPENDS CACHE)
elseif(UNIX)
  # ===== Look For BLAS =====
  find_package(BLAS)
  if(BLAS_FOUND)
    message( STATUS "BLAS found: ${BLAS_LIBRARIES}" )
    message( STATUS BLAS_HOME: ${BLAS_HOME} )
    message( STATUS BLAS_DIR: ${BLAS_DIR} )
    message( STATUS BLAS_LIBDIR: ${BLAS_LIBDIR} )
    list(GET BLAS_LIBRARIES 0 BLAS_LIBRARY_DIR)
    find_path(BLAS_INCLUDE_DIRS 
      NAMES 
        cblas.h
      HINTS
        ${_possible_paths}
        $ENV{BLAS_HOME}/include)
    list(GET BLAS_INCLUDE_DIRS 0 BLAS_INCLUDE_DIR)
    INCLUDE_DIRECTORIES(${BLAS_INCLUDE_DIR})
    message( STATUS "cblas.h => ${BLAS_INCLUDE_DIR}" )
  else()
    message( FATAL_ERROR "Missing BLAS.")
  endif()

  # ===== Look For Lapack =====
  find_package(LAPACK)
  if(LAPACK_FOUND)
    message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
    message( STATUS LAPACK_HOME: ${LAPACK_HOME} )
    message( STATUS LAPACK_DIR: ${LAPACK_DIR} )
    message( STATUS LAPACK_LIBDIR: ${LAPACK_LIBDIR} )
    find_path(CLAPACK_INCLUDE_DIRS 
      NAMES 
        clapack.h
      PATHS
        ${_possible_paths}
        $ENV{LAPACK_HOME}/include)
      if(NOT CLAPACK_INCLUDE_DIRS)
        find_path(LAPACKE_INCLUDE_DIRS
        NAMES
          lapacke.h
        PATHS
          /usr/include
          /usr/local/include
          $ENV{LAPACK_HOME}/include)
        if(NOT LAPACKE_INCLUDE_DIRS)
          message( FATAL_ERROR "LAPACK c interface headers missing.")
        else()
          message( STATUS "lapacke.h => ${LAPACKE_INCLUDE_DIRS}" )
          list(GET LAPACKE_INCLUDE_DIRS 0 LAPACK_INCLUDE_DIR)
          option(HAVE_LAPACKE "LAPACKE" ON)
          option(HAVE_CLAPACK OFF)
        endif()
      else()
        message( STATUS "clapack.h => ${CLAPACK_INCLUDE_DIRS}" )
        list(GET CLAPACK_INCLUDE_DIRS 0 LAPACK_INCLUDE_DIR)
        option(HAVE_CLAPACK "CLAPACK" ON)
          option(HAVE_LAPACKE OFF)
      endif()
  else()
    message( FATAL_ERROR "Missing LAPACK.")
  endif()

  INCLUDE_DIRECTORIES(${LAPACK_INCLUDE_DIR})
  # else()
  #   message( STATUS "Use OpenBLAS-0.3.6" )
  #   configure_file("${PROJECT_SOURCE_DIR}/dependencies/OpenBLAS-0.3.6/cmake/OpenBLASConfig.cmake.in" "${PROJECT_BINARY_DIR}/OpenBLASConfig.cmake")
  #   include(${PROJECT_BINARY_DIR}/OpenBLASConfig.cmake)
  #   find_package(BLAS REQUIRED)
  #   if(BLAS_FOUND)
  #     message("Found Blas located at:")
  #     list(GET BLAS_LIBRARIES 0 BLAS_LIBRARY_DIR)
  #     message("${BLAS_LIBRARY_DIR}")
  #     find_path(BLAS_INCLUDE_DIRS cblas.h
  #       ${BLAS_LIBRARY_DIR}/include)
  #     list(GET BLAS_INCLUDE_DIRS 0 BLAS_INCLUDE_DIR)
  #     message("Blas Include Dir: ${BLAS_INCLUDE_DIR}")
  #   endif()
  # endif()
else()
endif()

if(HAVE_APPLE_LAPACK)
  message(STATUS "Link veclib lapack...")
  target_link_libraries(matrix PRIVATE cblas clapack)
else()
  message(STATUS "Link lapacke...")
  target_link_libraries(matrix PRIVATE cblas lapacke)
endif()

target_include_directories(matrix
  PRIVATE
  "${BLAS_INCLUDE_DIR}"
  "${LAPACK_INCLUDE_DIR}"
  PUBLIC
  "${CMAKE_CURRENT_LIST_DIR}/include"
)

include("include/CMakeLists.txt")
include("src/CMakeLists.txt")

configure_file(${CMAKE_CURRENT_LIST_DIR}/configure_matrix.h.in ${PROJECT_BINARY_DIR}/configure_matrix.h)
message( STATUS "configure_matrix.h => ${PROJECT_BINARY_DIR}" )

include_directories(${PROJECT_BINARY_DIR})

