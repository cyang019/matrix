#include "matrix.h"
#include "gtest/gtest.h"
#include <limits>

namespace {
    using MatrixD = matrix::Matrix<double>; 
    using MatrixCx = matrix::Matrix<std::complex<double>>; 

    TEST(TestMatrix, DoubleMinus){
      MatrixD m1 = matrix::ones<double>(4,4);
      MatrixD m2 = {
        {1.0, 2.0, 3.0, 4.0},
        {5.0, 6.0, 7.0, 8.0},
        {9.0, 10.0, 11.0, 12.0},
        {13.0, 14.0, 15.0, 16.0}
      };

      MatrixD m_desired = {
          {0.0, 1.0, 2.0, 3.0},
          {4.0, 5.0, 6.0, 7.0},
          {8.0, 9.0, 10.0, 11.0},
          {12.0, 13.0, 14.0, 15.0}
      };

      ASSERT_TRUE(matrix::allclose(m_desired, m2-m1, 1.0e-14));
      ASSERT_TRUE(matrix::allclose(m_desired, m2-1.0, 1.0e-14));
    }
}


