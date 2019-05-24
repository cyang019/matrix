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

    TEST(TestMatrix, DoubleInit){
      MatrixD d = {{1, 2}, {3, 4}, {5, 6}};

      EXPECT_EQ(1, d(0,0));
      EXPECT_EQ(2, d(0,1));
      EXPECT_EQ(3, d(1,0));
      EXPECT_EQ(4, d(1,1));
      EXPECT_EQ(5, d(2,0));
      EXPECT_EQ(6, d(2,1));
    }

    TEST(TestMatrix, DoubleTranspose){
      MatrixD d = {{1, 2}, {3, 4}, {5, 6}};
      EXPECT_EQ(3, d.nrows());
      EXPECT_EQ(2, d.ncols());
      auto d2 = d.t();
      EXPECT_EQ(1, d2(0,0));
      EXPECT_EQ(2, d2(1,0));
      EXPECT_EQ(3, d2(0,1));
      EXPECT_EQ(4, d2(1,1));
      EXPECT_EQ(5, d2(0,2));
      EXPECT_EQ(6, d2(1,2));

      EXPECT_EQ(2, d2.nrows());
      EXPECT_EQ(3, d2.ncols());
      d.tInplace();
      EXPECT_EQ(1, d(0,0));
      EXPECT_EQ(2, d(1,0));
      EXPECT_EQ(3, d(0,1));
      EXPECT_EQ(4, d(1,1));
      EXPECT_EQ(5, d(0,2));
      EXPECT_EQ(6, d(1,2));
      EXPECT_EQ(2, d.nrows());
      EXPECT_EQ(3, d.ncols());
    }

    TEST(TestMatrix, DoublePlus){
      MatrixD d3 = {{1, 2}, {3, 4}};
      MatrixD d4 = {{1, 2}, {3, 4}};
      MatrixD d5 = d3 + d4;
      EXPECT_EQ(2, d5(0, 0));
      EXPECT_EQ(4, d5(0, 1));
      EXPECT_EQ(6, d5(1, 0));
      EXPECT_EQ(8, d5(1, 1));
      d3 += d4;
      EXPECT_EQ(2, d3(0, 0));
      EXPECT_EQ(4, d3(0, 1));
      EXPECT_EQ(6, d3(1, 0));
      EXPECT_EQ(8, d3(1, 1));

      d3 += 2;
      EXPECT_EQ(4, d3(0, 0));
      EXPECT_EQ(6, d3(0, 1));
      EXPECT_EQ(8, d3(1, 0));
      EXPECT_EQ(10, d3(1, 1));
    }
}
