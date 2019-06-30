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

    TEST(TestMatrix, DoubleExp){
      MatrixD m1= matrix::diagonal({1.0,2.0,3.0,4.0});
      MatrixD res1 = matrix::exp(m1);
      MatrixD desired = matrix::diagonal(
          { std::exp(1.0), std::exp(2.0), std::exp(3.0), std::exp(4.0)});
      std::cout << std::setprecision(15);
      std::cout << "desired:\n" << desired << "\n";
      std::cout << "result:\n" << res1 << "\n";

      double eps = matrix::norm1(desired) * matrix::eps * 5;
      std::cout << "desired matrix norm1: " << matrix::norm1(desired) << "\n";

      ASSERT_TRUE(matrix::allclose(desired, res1, eps));
    }
}

