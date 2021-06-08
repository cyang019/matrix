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

if(MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  target_compile_options(matrix PUBLIC -Wall)
else()
  target_compile_options(matrix PUBLIC -Wall -Wextra -Wpedantic)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

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

if(MSVC)
  set(MKL_PATH_MATRIX "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
  option(HAVE_LAPACKE "LAPACKE" ON)
  target_include_directories(matrix
    PUBLIC
    "${CMAKE_CURRENT_LIST_DIR}/include"
    "${MKL_PATH_MATRIX}/include"
  )

  target_link_libraries(matrix PUBLIC
    ${MKL_PATH_MATRIX}/lib/intel64/mkl_intel_lp64_dll.lib
    ${MKL_PATH_MATRIX}/lib/intel64/mkl_sequential_dll.lib
    ${MKL_PATH_MATRIX}/lib/intel64/mkl_core_dll.lib
    )
else()
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

    if(APPLE)
      target_link_libraries(matrix PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    else()
      target_link_libraries(matrix PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} lapacke)
    endif()
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
endif()


include("include/CMakeLists.txt")
include("src/CMakeLists.txt")

configure_file(${CMAKE_CURRENT_LIST_DIR}/configure_matrix.h.in ${PROJECT_BINARY_DIR}/configure_matrix.h)
message( STATUS "configure_matrix.h => ${PROJECT_BINARY_DIR}" )

include_directories(${PROJECT_BINARY_DIR})

# copy third party dependencies
if(MSVC)
  set(THIRD_PARTY_MKL_DIR "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/redist/intel64")
  set(MKL_DLLS
    "${THIRD_PARTY_MKL_DIR}/mkl_avx.1.dll"
    "${THIRD_PARTY_MKL_DIR}/mkl_core.1.dll"
    "${THIRD_PARTY_MKL_DIR}/mkl_def.1.dll"
    "${THIRD_PARTY_MKL_DIR}/mkl_sequential.1.dll"
  )
  foreach( file_i ${MKL_DLLS})
    add_custom_command(
      TARGET matrix
      POST_BUILD
      COMMAND ${CMAKE_COMMAND}
      ARGS -E copy ${file_i} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} 
    )
  endforeach( file_i )
endif()
