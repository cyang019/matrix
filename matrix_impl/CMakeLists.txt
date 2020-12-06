if(NOT USE_EMCC)
  add_library(matrix SHARED "")
else()
  add_library(matrix STATIC "")
endif()

if(NOT DEFINED USE_DEBUG_MODEL) 
  set(USE_DEBUG_MODEL ON)
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_DEBUG)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
endif()

if(NOT DEFINED CMAKE_CXX_FLAGS_RELEASE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
endif()

message(STATUS "CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}")

target_compile_options(matrix PUBLIC -Wall -Wextra -Wpedantic)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

#set(_possible_paths ${CMAKE_PREFIX_PATH}/include /usr/include /usr/local/include /usr/include/x86_64-linux-gnu)
#
#if(USE_EMCC)
#  #find_library(CBLAS ${CMAKE_CURRENT_LIST_DIR}/libopenblasp-r0.3.7.a)
#  set(CBLAS_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/include")
#  set(LAPACK_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/include")
#else()
#  find_package(OPENBLAS)
#  if(OPENBLAS_FOUND)
#    message( STATUS "OPENBLAS: ${OPENBLAS}" )
#    message( STATUS "OPENBLAS_LIBRARIES: ${OPENBLAS_LIBRARIES}" )
#    message( STATUS "OPENBLAS_INCLUDE_DIRS: ${OPENBLAS_INCLUDE_DIRS}" )
#    list(GET OPENBLAS_LIBRARIES 0 OPENBLAS_LIBRARY_DIR)
#    list(GET OPENBLAS_INCLUDE_DIRS  0 OPENBLAS_INCLUDE_DIR)
#    include_directories(${OPENBLAS_INCLUDE_DIR})
#    message( STATUS "openblas/cblas.h => ${OPENBLAS_INCLUDE_DIR}" )
#  else()
#    find_package(CBLAS)
#    if(NOT CBLAS_FOUND)
#      message(FATAL_ERROR "Need at least one BLAS library")
#    else()
#      message( STATUS "CBLAS: ${CBLAS}" )
#      message( STATUS "CBLAS_LIBRARIES: ${CBLAS_LIBRARIES}" )
#      message( STATUS "CBLAS_INCLUDE_DIRS: ${CBLAS_INCLUDE_DIRS}" )
#      list(GET CBLAS_LIBRARIES 0 CBLAS_LIBRARY_DIR)
#      list(GET CBLAS_INCLUDE_DIRS  0 CBLAS_INCLUDE_DIR)
#      include_directories(${CBLAS_INCLUDE_DIR})
#      message( STATUS "cblas.h => ${CBLAS_INCLUDE_DIR}" )
#    endif()
#  endif()
#endif()
#
##  find_path(BLAS_INCLUDE_DIRS 
##    NAMES 
##      cblas.h
##    HINTS
##      /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
##      ${_possible_paths}
##      $ENV{BLAS_HOME}/include)
##  list(GET BLAS_INCLUDE_DIRS 0 BLAS_INCLUDE_DIR)
#
#if(USE_EMCC)
#  # these are NOT tested.
#  target_link_libraries(matrix PRIVATE CBLAS)
#  option(HAVE_LAPACKE "LAPACKE" ON)
#else()
#  if (APPLE)
#    message( STATUS "OSX detected..." )
#  
#    # ===== Look For Lapack =====
#    find_package(LAPACK REQUIRED)
#    message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
#    message( STATUS LAPACK_HOME: ${LAPACK_HOME} )
#    message( STATUS LAPACK_DIR: ${LAPACK_DIR} )
#    message( STATUS LAPACK_LIBDIR: ${LAPACK_LIBDIR} )
#    find_path(CLAPACK_INCLUDE_DIRS 
#      NAMES 
#        clapack.h
#      PATHS
#        /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
#        ${_possible_paths}
#        $ENV{LAPACK_HOME}/include)
#    #set(BLAS_INCLUDE_DIR
#    #  "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers"
#    #)
#    option(HAVE_APPLE_LAPACK "APPLE_LAPACK" ON)
#  
#    #add_compile_definitions(WITH_VECLIB_CBLAS)
#    #unset(matrix_LIB_DEPENDS CACHE)
#  elseif(UNIX)
#    message( STATUS "Linux detected..." )
#  
#    # ===== Look For Lapack =====
#    find_package(LAPACK REQUIRED)
#    message( STATUS "LAPACK found: ${LAPACK_LIBRARIES}" )
#    message( STATUS LAPACK_HOME: ${LAPACK_HOME} )
#    message( STATUS LAPACK_DIR: ${LAPACK_DIR} )
#    message( STATUS LAPACK_LIBDIR: ${LAPACK_LIBDIR} )
#    find_path(CLAPACK_INCLUDE_DIRS 
#      NAMES 
#        clapack.h
#      PATHS
#        ${_possible_paths}
#        $ENV{LAPACK_HOME}/include)
#      if(NOT CLAPACK_INCLUDE_DIRS)
#        find_path(LAPACKE_INCLUDE_DIRS
#        NAMES
#          lapacke.h
#        PATHS
#          /usr/include
#          /usr/local/include
#          $ENV{LAPACK_HOME}/include)
#        if(NOT LAPACKE_INCLUDE_DIRS)
#          message( FATAL_ERROR "LAPACK c interface headers missing.")
#        else()
#          message( STATUS "lapacke.h => ${LAPACKE_INCLUDE_DIRS}" )
#          list(GET LAPACKE_INCLUDE_DIRS 0 LAPACK_INCLUDE_DIR)
#          option(HAVE_LAPACKE "LAPACKE" ON)
#          option(HAVE_CLAPACK OFF)
#        endif()
#      else()
#        message( STATUS "clapack.h => ${CLAPACK_INCLUDE_DIRS}" )
#        list(GET CLAPACK_INCLUDE_DIRS 0 LAPACK_INCLUDE_DIR)
#        option(HAVE_CLAPACK "CLAPACK" ON)
#          option(HAVE_LAPACKE OFF)
#      endif()
#    include_directories(${LAPACK_INCLUDE_DIR})
#  else()
#  endif()
# 
#  if(HAVE_APPLE_LAPACK)
#    message(STATUS "Link veclib lapack...")
#    target_link_libraries(matrix PRIVATE cblas clapack)
#  elseif(HAVE_CLAPACK)
#    message(STATUS "Link clapack...")
#    if(OPENBLAS_FOUND)
#      target_link_libraries(matrix PRIVATE openblas clapack)
#    else()
#      target_link_libraries(matrix PRIVATE cblas clapack)
#    endif()
#  else()
#    message(STATUS "Link lapacke...")
#    if(OPENBLAS_FOUND)
#      target_link_libraries(matrix PRIVATE openblas lapacke)
#    else()
#      target_link_libraries(matrix PRIVATE cblas lapacke)
#    endif()
#  endif()
#endif()

if(APPLE)
  set(BLA_VENDER Apple)
elseif(MSVC)
  set(BLA_VENDER Intel10_64lp)
endif()

set(CBLAS_PATHS
  /usr/local
  /usr/local/opt
  /usr
  /opt
  /System/Library/Frameworks
  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks
  /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.15.sdk/System/Library/Frameworks
)

find_path(CBLAS_INCLUDE_DIRS
  NAMES cblas.h
  PATH_SUFFIXES
    openblas/include include/x86_64-linux-gnu
    include inc include/x86_64 include/x64
    include/blis blis/include blis/include/blis
    Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers
    Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
  PATHS ${CBLAS_PATHS}
  DOC "BLAS include header cblas.h"
)
message(STATUS "CBLAS_INCLUDE_DIRS => ${CBLAS_INCLUDE_DIRS}")

find_path(LAPACKE_INCLUDE_DIRS
  NAMES lapacke.h
  PATH_SUFFIXES
    openblas/include include/x86_64-linux-gnu
    include inc include/x86_64 include/x64 
    include/blis blis/include blis/include/blis
    Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers
    Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
  PATHS ${CBLAS_PATHS}
  DOC "Netlib LAPACK include header lapacke.h"
)

find_path(LAPACK_DIR
  NAMES lapacke.h clapack.h
  PATH_SUFFIXES
    openblas/include include/x86_64-linux-gnu
    include inc include/x86_64 include/x64 
    include/blis blis/include blis/include/blis
    Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers
    Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
  PATHS ${CBLAS_PATHS}
  DOC "Netlib LAPACK include header lapacke.h"
)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(NOT LAPACKE_INCLUDE_DIRS)
  option(HAVE_LAPACKE "LAPACKE" OFF)
  find_path(CLAPACK_INCLUDE_DIRS REQUIRED
    NAMES clapack.h
    PATH_SUFFIXES
      openblas/include include/x86_64-linux-gnu
      include inc include/x86_64 include/x64
      include/blis blis/include blis/include/blis
      Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers
      Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers
    PATHS ${CBLAS_PATHS}
    DOC "Netlib LAPACK include header clapack.h"
  )
  if(NOT CLAPACK_INCLUDE_DIRS)
    message(FATAL_ERROR "CLAPACK_INCLUDE_DIRS not found.")
  endif()
  message(STATUS "CLAPACK_INCLUDE_DIRS => ${CLAPACK_INCLUDE_DIRS}")
  option(HAVE_CLAPACK "CLAPACK" ON)
else()
  option(HAVE_CLAPACK "CLAPACK" OFF)
  option(HAVE_LAPACKE "LAPACKE" ON)
  message(STATUS "LAPACKE_INCLUDE_DIRS => ${LAPACKE_INCLUDE_DIRS}")
endif()

if(APPLE)
  option(HAVE_APPLE_LAPACK "APPLE_LAPACK" ON)
endif()


message(STATUS "BLAS_LIBRARIES => ${BLAS_LIBRARIES}")
message(STATUS "LAPACK_LIBRARIES => ${LAPACK_LIBRARIES}")

if(HAVE_LAPACKE)
  target_include_directories(matrix
    PRIVATE
    PUBLIC
    "${CMAKE_CURRENT_LIST_DIR}/include"
    "${CBLAS_INCLUDE_DIRS}"
    "${LAPACKE_INCLUDE_DIRS}"
  )

  target_link_libraries(matrix PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} lapacke)
else()
  target_include_directories(matrix
    PRIVATE
    PUBLIC
    "${CMAKE_CURRENT_LIST_DIR}/include"
    "${CBLAS_INCLUDE_DIRS}"
    "${CLAPACK_INCLUDE_DIRS}"
  )

  target_link_libraries(matrix PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} clapack)
endif()



include("include/CMakeLists.txt")
include("src/CMakeLists.txt")

configure_file(${CMAKE_CURRENT_LIST_DIR}/configure_matrix.h.in ${PROJECT_BINARY_DIR}/configure_matrix.h)
message( STATUS "configure_matrix.h => ${PROJECT_BINARY_DIR}" )

include_directories(${PROJECT_BINARY_DIR})

