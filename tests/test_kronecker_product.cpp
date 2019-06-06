#include "Matrix.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>


namespace {
    using MatrixD = matrix::Matrix<double>; 

    TEST(TestMatrix, DoubleOuter){
      MatrixD m1 = {{1,2}, {3,4}};
      MatrixD m2 = {{0,5}, {6,7}};

      MatrixD m12 = matrix::kroneckerProduct(m1, m2);
      MatrixD m12_desired = {{0, 5, 0, 10},
                             {6, 7, 12, 14},
                             {0, 15, 0, 20},
                             {18, 21, 24, 28}};

      //std::cout << m12;
      //std::cout << m12_desired;
      ASSERT_TRUE(matrix::allclose(m12_desired, m12, 1.0e-14));
    }
}

