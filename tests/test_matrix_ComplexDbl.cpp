#include "Matrix.h"
#include "gtest/gtest.h"
#include <complex>

namespace {
    using MatrixCxDbl = matrix::Matrix<matrix::ComplexDbl>;

    TEST(TestMatrix, SetComplexZero){
        MatrixCxDbl m(2,2);
        m.setZero();
        ASSERT_DOUBLE_EQ(0.0, m(0,0).r);
        ASSERT_DOUBLE_EQ(0.0, m(0,1).r);
        ASSERT_DOUBLE_EQ(0.0, m(1,0).r);
        ASSERT_DOUBLE_EQ(0.0, m(1,1).r);
        ASSERT_DOUBLE_EQ(0.0, m(0,0).i);
        ASSERT_DOUBLE_EQ(0.0, m(0,1).i);
        ASSERT_DOUBLE_EQ(0.0, m(1,0).i);
        ASSERT_DOUBLE_EQ(0.0, m(1,1).i);
    }
}

