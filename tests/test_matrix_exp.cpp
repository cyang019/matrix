#include "Matrix.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>    // exp


namespace {
    using MatrixD = matrix::Matrix<double>; 
    using MatrixI = matrix::Matrix<int>; 

    TEST(TestMatrix, DoubleExp){
      MatrixD m1= matrix::diagonal({1.0,2.0,3.0,4.0});
      MatrixD res1 = matrix::exp(m1);
      MatrixD desired = matrix::diagonal(
          { std::exp(1.0), std::exp(2.0), std::exp(3.0), std::exp(4.0)});
      ASSERT_TRUE(matrix::allclose(desired, m1, 1.0e-14));
    }
}

