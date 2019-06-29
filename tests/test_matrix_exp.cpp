#include "Matrix.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>    // exp


namespace {
    using MatrixD = matrix::Matrix<double>; 
    using MatrixI = matrix::Matrix<int>; 

    TEST(TestMatrix, areParallel){
      MatrixD mat = {
        {1, 2, -3, 4},
        {-1, 2, 3, 4},
        {1, 2, -3, 4}
      };
      
      ASSERT_EQ(true, matrix::exponential::areParallel(mat, 0, mat, 2));
      ASSERT_EQ(true, matrix::exponential::areParallel(mat, 1, mat, 3));
      ASSERT_EQ(false, matrix::exponential::areParallel(mat, 0, mat, 1));
    }

    TEST(TestMatrix, Normest){
      MatrixD m1= matrix::diagonal({1.0,2.0,3.0,4.0});
      double val = matrix::exponential::normest(m1, 1);
      ASSERT_TRUE(val < 4.0);
      val = matrix::norm1(m1);
      EXPECT_DOUBLE_EQ(4.0, val);
      val = matrix::exponential::normest(m1, 2);
      ASSERT_TRUE(val < 16.0);
    }

    // TEST(TestMatrix, DoubleExp){
    //   MatrixD m1= matrix::diagonal({1.0,2.0,3.0,4.0});
    //   MatrixD res1 = matrix::exp(m1);
    //   MatrixD desired = matrix::diagonal(
    //       { std::exp(1.0), std::exp(2.0), std::exp(3.0), std::exp(4.0)});
    //   ASSERT_TRUE(matrix::allclose(desired, m1, 1.0e-14));
    // }
}

