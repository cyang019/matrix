#include "matrix.h"
#include "matrix_core/common.h"
#include "gtest/gtest.h"
#include <complex>
#include <iostream>

namespace {
    using namespace std;
    using matrix::log2int;

    TEST(TestMatrix, Pow){
      auto mat_identity = matrix::identity<double>(5);
      auto mat_start = mat_identity * 2.0;

      auto mat_desired = mat_identity * 1024.0;
      
      ASSERT_TRUE(matrix::allclose(mat_desired, matrix::pow(mat_start, 10), 
            1.0e-14));
    }
}

