#include "matrix.h"
#include "gtest/gtest.h"
#include <limits>

namespace {
    using MatrixD = matrix::Matrix<double>; 
    using MatrixCx = matrix::Matrix<std::complex<double>>; 

    TEST(TestMatrix, DoubleScale){
      MatrixD m1 = {
        {1.0, 2.0, 3.0, 4.0},
        {5.0, 6.0, 7.0, 8.0},
        {9.0, 10.0, 11.0, 12.0},
        {13.0, 14.0, 15.0, 16.0}
      };
      MatrixD m2 = m1 * 2.0;
      MatrixD m3 = 2.0 * m1;

      m1 *= 2.0;

      std::cout << "m1: " << m1 << "\n";
      std::cout << "m1 * 2.0:\n" << m2 << "\n";

      MatrixD m_desired = {
        {2.0, 4.0, 6.0, 8.0},
        {10.0, 12.0, 14.0, 16.0},
        {18.0, 20.0, 22.0, 24.0},
        {26.0, 28.0, 30.0, 32.0}
      };

      ASSERT_TRUE(matrix::allclose(m_desired, m1, 1.0e-14));
      ASSERT_TRUE(matrix::allclose(m_desired, m2, 1.0e-14));
      ASSERT_TRUE(matrix::allclose(m_desired, m3, 1.0e-14));
    }
}

