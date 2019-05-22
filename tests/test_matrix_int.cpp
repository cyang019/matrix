#include "Matrix.h"
#include "gtest/gtest.h"

namespace {
    using MatrixI = matrix::Matrix<int>;
    TEST(TestMatrix, IntZero){
        MatrixI m(3,3);
        m.setZero();
        EXPECT_EQ(0, m(0,0));
        EXPECT_EQ(0, m(0,1));
        EXPECT_EQ(0, m(0,2));
        EXPECT_EQ(0, m(1,0));
        EXPECT_EQ(0, m(1,1));
        EXPECT_EQ(0, m(1,2));
        EXPECT_EQ(0, m(2,0));
        EXPECT_EQ(0, m(2,1));
        EXPECT_EQ(0, m(2,2));
    }

    TEST(TestMatrix, IntOne){
        MatrixI m(2,2);
        m.setOne();
        EXPECT_EQ(1, m(0,0));
        EXPECT_EQ(1, m(0,1));
        EXPECT_EQ(1, m(1,0));
        EXPECT_EQ(1, m(1,1));
    }

    TEST(TestMatrix, IntPlus){
        MatrixI m(2,2);
        m(0,0) = 2;
        m(0,1) = 3;
        m(1,0) = 5;
        m(1,1) = 7;

        MatrixI m2(2,2);
        m2(0,0) = 1;
        m2(0,1) = 2;
        m2(1,0) = 3;
        m2(1,1) = 4;

        MatrixI m3 = m + m2;
        EXPECT_EQ(3, m3(0,0));
        EXPECT_EQ(5, m3(0,1));
        EXPECT_EQ(8, m3(1,0));
        EXPECT_EQ(11, m3(1,1));

        m += m2;
        EXPECT_EQ(3, m(0,0));
        EXPECT_EQ(5, m(0,1));
        EXPECT_EQ(8, m(1,0));
        EXPECT_EQ(11, m(1,1));

    }

    TEST(TestMatrix, InitMat){
        MatrixI m = {{1,2,3},{4,5,6},{7,8,9}};
        EXPECT_EQ(1, m(0,0));
        EXPECT_EQ(2, m(0,1));
        EXPECT_EQ(3, m(0,2));
        EXPECT_EQ(4, m(1,0));
        EXPECT_EQ(5, m(1,1));
        EXPECT_EQ(6, m(1,2));
        EXPECT_EQ(7, m(2,0));
        EXPECT_EQ(8, m(2,1));
        EXPECT_EQ(9, m(2,2));
    }

    TEST(TestMatrix, AssignMat){
        MatrixI m1 = {{1,2,3}, {4,5,6}, {7,8,9}};
        MatrixI m2(m1.nrows(), m1.ncols());
        for(size_t i = 0; i < m2.nrows(); ++i){
            for(size_t j = 0; j < m2.ncols(); ++j){
                m2(i,j) = m1(i,j);
            }
        }
        EXPECT_EQ(1, m2(0,0));
        EXPECT_EQ(2, m2(0,1));
        EXPECT_EQ(3, m2(0,2));
        EXPECT_EQ(4, m2(1,0));
        EXPECT_EQ(5, m2(1,1));
        EXPECT_EQ(6, m2(1,2));
        EXPECT_EQ(7, m2(2,0));
        EXPECT_EQ(8, m2(2,1));
        EXPECT_EQ(9, m2(2,2));
    }
}   // namespace
