#include "matrix.h"
#include "matrix_core/Matrix.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>


namespace {
    using MatrixD = matrix::Matrix<double>; 

    TEST(TestMatrix, DoubleKron){
      MatrixD m1 = {{1,2}, {3,4}};
      MatrixD m2 = {{0,5}, {6,7}};

      MatrixD m12 = matrix::kroneckerProduct(m1, m2);
      MatrixD m12_desired = {{0, 5, 0, 10},
                             {6, 7, 12, 14},
                             {0, 15, 0, 20},
                             {18, 21, 24, 28}};

      std::cout << "m12:\n" << m12;
      std::cout << "m12 desired:\n" << m12_desired;
      ASSERT_TRUE(matrix::allclose(m12_desired, m12, 1.0e-14));

      MatrixD m3 = {{1, -4, 7}, {-2, 3, 3}};
      MatrixD m4 = {{8, -9, -6, 5},
                    {1, -3, -4, 7},
                    {2, 8, -8, -3},
                    {1, 2, -5, -1}};
      MatrixD m34 = matrix::kroneckerProduct(m3, m4);
      MatrixD m34_desired = {
        {8,-9,-6,5,-32,36,24,-20,56,-63,-42,35},
        {1,-3,-4,7,-4,12,16,-28,7,-21,-28,49},
        {2,8,-8,-3,-8,-32,32,12,14,56,-56,-21},
        {1,2,-5,-1,-4,-8,20,4,7,14,-35,-7},
        {-16,18,12,-10,24,-27,-18,15,24,-27,-18,15},
        {-2,6,8,-14,3,-9,-12,21,3,-9,-12,21},
        {-4,-16,16,6,6,24,-24,-9,6,24,-24,-9},
        {-2,-4,10,2,3,6,-15,-3,3,6,-15,-3}};
      std::cout << "m34:\n" << m34;
      std::cout << "m34 desired:\n" << m34_desired;
      ASSERT_TRUE(matrix::allclose(m34_desired, m34, 1.0e-24));
    }
}

