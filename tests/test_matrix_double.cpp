#include "Matrix.h"
#include "gtest/gtest.h"
#include <limits>

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

    TEST(TestMatrix, DoubleCopy){
      MatrixD d(2,2);
      d.setOne();
      // call ctor
      MatrixD d2(d);
      constexpr double eps = std::numeric_limits<double>::epsilon();
      ASSERT_NEAR(1.0, d2(0,0), eps);
      ASSERT_NEAR(1.0, d2(0,1), eps);
      ASSERT_NEAR(1.0, d2(1,0), eps);
      ASSERT_NEAR(1.0, d2(1,1), eps);
    }
}
