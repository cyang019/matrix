#include(GoogleTest)
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
    ${PROJECT_SOURCE_DIR}/matrix_impl/include
  )

#gtest_add_tests(TARGET      test_matrix
#                TEST_SUFFIX .noArgs
#                TEST_LIST   noArgsTests
#                )
#set_tests_properties(${noArgsTests} PROPERTIES TIMEOUT 10)
add_test(NAME test_matrix
         COMMAND test_matrix)

