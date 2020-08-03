file(GLOB SRCS ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)
add_executable(benchmark_matrix ${SRCS})

target_link_libraries(
  benchmark_matrix
  matrix
  )

target_include_directories(
  benchmark_matrix
  PUBLIC
  ${PROJECT_SOURCE_DIR}/matrix_impl/include
  ${PROJECT_SOURCE_DIR}/benchmark_matrix/include
  )

target_compile_options(test_matrix PUBLIC
  -Wall)

