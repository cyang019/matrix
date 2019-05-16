file(GLOB SRCS *.cpp)
add_executable(test_matrix ${SRCS})

target_link_libraries(
  test_matrix
  matrix
  gtest_main
  )

target_include_directories(
  test_matrix
  PUBLIC
    ${PROJECT_BINARY_DIR}/googletest/googletest-src/googletest/include
  )

add_test(NAME test_matrix
         COMMAND test_matrix)
