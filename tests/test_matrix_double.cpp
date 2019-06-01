#include "Matrix.h"
#include "gtest/gtest.h"
#include <limits>

namespace {
    using MatrixD = matrix::Matrix<double>; 

    TEST(TestMatrix, DoubleZero){
      MatrixD m(2,2);
      m.setZero();
      EXPECT_DOUBLE_EQ(0, m(0,0));
      EXPECT_DOUBLE_EQ(0, m(0,1));
      EXPECT_DOUBLE_EQ(0, m(1,0));
      EXPECT_DOUBLE_EQ(0, m(1,1));
    }

    TEST(TestMatrix, DoubleOne){
      MatrixD m(2,2);
      m.setOne();
      EXPECT_DOUBLE_EQ(1, m(0,0));
      EXPECT_DOUBLE_EQ(1, m(0,1));
      EXPECT_DOUBLE_EQ(1, m(1,0));
      EXPECT_DOUBLE_EQ(1, m(1,1));
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

      EXPECT_DOUBLE_EQ(1, d(0,0));
      EXPECT_DOUBLE_EQ(2, d(0,1));
      EXPECT_DOUBLE_EQ(3, d(1,0));
      EXPECT_DOUBLE_EQ(4, d(1,1));
      EXPECT_DOUBLE_EQ(5, d(2,0));
      EXPECT_DOUBLE_EQ(6, d(2,1));
    }

    TEST(TestMatrix, DoubleTranspose){
      MatrixD d = {{1, 2}, {3, 4}, {5, 6}};
      EXPECT_DOUBLE_EQ(3, d.nrows());
      EXPECT_DOUBLE_EQ(2, d.ncols());
      auto d2 = d.t();
      EXPECT_DOUBLE_EQ(1, d2(0,0));
      EXPECT_DOUBLE_EQ(2, d2(1,0));
      EXPECT_DOUBLE_EQ(3, d2(0,1));
      EXPECT_DOUBLE_EQ(4, d2(1,1));
      EXPECT_DOUBLE_EQ(5, d2(0,2));
      EXPECT_DOUBLE_EQ(6, d2(1,2));

      EXPECT_DOUBLE_EQ(2, d2.nrows());
      EXPECT_DOUBLE_EQ(3, d2.ncols());
      d.tInplace();
      EXPECT_DOUBLE_EQ(1, d(0,0));
      EXPECT_DOUBLE_EQ(2, d(1,0));
      EXPECT_DOUBLE_EQ(3, d(0,1));
      EXPECT_DOUBLE_EQ(4, d(1,1));
      EXPECT_DOUBLE_EQ(5, d(0,2));
      EXPECT_DOUBLE_EQ(6, d(1,2));
      EXPECT_DOUBLE_EQ(2, d.nrows());
      EXPECT_DOUBLE_EQ(3, d.ncols());
    }

    TEST(TestMatrix, DoublePlus){
      MatrixD d3 = {{1, 2}, {3, 4}};
      MatrixD d4 = {{1, 2}, {3, 4}};
      MatrixD d5 = d3 + d4;
      EXPECT_DOUBLE_EQ(2, d5(0, 0));
      EXPECT_DOUBLE_EQ(4, d5(0, 1));
      EXPECT_DOUBLE_EQ(6, d5(1, 0));
      EXPECT_DOUBLE_EQ(8, d5(1, 1));
      d3 += d4;
      EXPECT_DOUBLE_EQ(2, d3(0, 0));
      EXPECT_DOUBLE_EQ(4, d3(0, 1));
      EXPECT_DOUBLE_EQ(6, d3(1, 0));
      EXPECT_DOUBLE_EQ(8, d3(1, 1));

      d3 += 2;
      EXPECT_DOUBLE_EQ(4, d3(0, 0));
      EXPECT_DOUBLE_EQ(6, d3(0, 1));
      EXPECT_DOUBLE_EQ(8, d3(1, 0));
      EXPECT_DOUBLE_EQ(10, d3(1, 1));

      auto d6 = d5 + 2.0;
      EXPECT_DOUBLE_EQ(4, d6(0, 0));
      EXPECT_DOUBLE_EQ(6, d6(0, 1));
      EXPECT_DOUBLE_EQ(8, d6(1, 0));
      EXPECT_DOUBLE_EQ(10, d6(1, 1));

      auto d7 = 2.0 + d5;
      EXPECT_DOUBLE_EQ(4, d7(0, 0));
      EXPECT_DOUBLE_EQ(6, d7(0, 1));
      EXPECT_DOUBLE_EQ(8, d7(1, 0));
      EXPECT_DOUBLE_EQ(10, d7(1, 1));
    }

    TEST(TestMatrix, DoubleEqual){
      MatrixD d1 = {{1, 2}, {3, 4}};
      MatrixD d2 = d1;
      MatrixD d3 = d1 + d2;
      MatrixD d3_expected = {{2, 4}, {6, 8}};
      ASSERT_TRUE(allclose(d1, d2, 1.0e-14));
      ASSERT_TRUE(!allclose(d1, d3, 1.0e-14));
      ASSERT_TRUE(allclose(d3_expected, d3, 1.0e-14));
      ASSERT_TRUE(d1 != d3);
    }

    TEST(TestMatrix, DoubleMultiply){
      MatrixD d1 = {{1, 2, 3}};
      MatrixD d2 = {{4}, {5}, {6}};

      auto d3 = d1 * d2;
      EXPECT_DOUBLE_EQ(32, d3(0,0));

      d1 = {{1, 2, 3}, {4, 5, 6}};
      d3 = d1 * d2;
      MatrixD d3_expected = {{32.0}, {77}};
      ASSERT_TRUE(allclose(d3_expected, d3, 1.0e-14));

      MatrixD d4 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
      MatrixD d5 = {{1, 4}, {2, 5}, {3, 6}};
      auto d6 = d4 * d5;
      MatrixD d6_expected = {{14, 32}, {32, 77}, {50, 122}};
      ASSERT_TRUE(allclose(d6_expected, d6, 1.0e-14));
    }

    TEST(TestMatrix, DoubleInverse){
      MatrixD d1 = {{1.0, 0, 0}, {0, 2.0, 0}, {0, 0, 4.0}};
      auto d1_inv = d1.inverse();
      auto d2 = d1 * d1_inv;
      MatrixD d2_expected = {{1.0, 0, 0}, {0, 1.0, 0}, {0, 0, 1.0}};
      ASSERT_TRUE(allclose(d2_expected, d2, 1.0e-14));
    }
}
