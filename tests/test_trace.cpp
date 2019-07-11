#include "matrix.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>


namespace {
    using MatrixD = matrix::Matrix<double>; 
    using MatrixI = matrix::Matrix<int>; 

    TEST(TestMatrix, DoubleTrace){
      MatrixD m12= {{0, 5, 0, 10},
                    {6, 7, 12, 14},
                    {0, 15, 0, 20},
                    {18, 21, 24, 28}};
      double t = matrix::trace(m12);
      ASSERT_DOUBLE_EQ(35, t);
    }

    TEST(TestMatrix, IntTrace){
      MatrixI m12= {{0, 5, 0, 10},
                    {6, 7, 12, 14},
                    {0, 15, 0, 20},
                    {18, 21, 24, 28}};
      int t = matrix::trace(m12);
      ASSERT_EQ(35, t);
    }
}
