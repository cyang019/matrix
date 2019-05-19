#include "Matrix.h"
#include "gtest/gtest.h"
#include <complex>

namespace {
    using MatrixD = matrix::Matrix<double>; 

    TEST(TestMatrix, DoubleZero){
      MatrixD m(2,2);
      m.setZero();
      EXPECT_EQ(0, m(0,0));
      EXPECT_EQ(0, m(0,1));
      EXPECT_EQ(0, m(1,0));
      EXPECT_EQ(0, m(1,1));
    }

    TEST(TestMatrix, DoubleOne){
      MatrixD m(2,2);
      m.setOne();
      EXPECT_EQ(1, m(0,0));
      EXPECT_EQ(1, m(0,1));
      EXPECT_EQ(1, m(1,0));
      EXPECT_EQ(1, m(1,1));
    }

    TEST(TestMatrix, DoublePlus){
      MatrixD d(2,2);

    }
}
