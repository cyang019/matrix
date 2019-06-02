file(GLOB SRCS *.cpp)
add_executable(test_matrix ${SRCS})

set(USE_COVERAGE OFF)

target_link_libraries(
  test_matrix
  matrix
  gtest_main
  )

if(USE_COVERAGE)
  if(APPLE)
    message("APPLE coverage...")
    target_compile_options(test_matrix PUBLIC
      -g -O0 --coverage)
  else()
    target_compile_options(test_matrix PUBLIC
      -Wall -fprofile-arcs -ftest-coverage
    )
    target_link_libraries(test_matrix gcov)
  endif()
  add_custom_command(TARGET test_matrix POST_BUILD
    COMMAND make test
    COMMAND lcov -c -d tests/CMakeFiles -o cov.info
    COMMAND genhtml cov.info -o out
    COMMAND open out/index.html
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )

else()
  target_compile_options(test_matrix PUBLIC
    -Wall)
endif()


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

