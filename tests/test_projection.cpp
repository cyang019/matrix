#include "Matrix.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>


namespace {
    using MatrixD = matrix::Matrix<double>; 

    TEST(TestMatrix, DoubleProjection){
      MatrixD d1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
      MatrixD d2 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
      double p = matrix::projection(d1, d2);
      ASSERT_DOUBLE_EQ(285.0, p);
    }

    TEST(TestMatrix, DoubleProjectionNorm){
      MatrixD d1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
      MatrixD d2 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
      double p = matrix::projectionNorm(d1, d2);
      ASSERT_DOUBLE_EQ(1.0, p);
    }
}
