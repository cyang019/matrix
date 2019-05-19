#include "Matrix.h"
#include "gtest/gtest.h"
#include <complex>

namespace {
    using MatrixCD = matrix::Matrix<std::complex<double>>; 

    TEST(TestMatrix, SetComplexZero){
        MatrixCD m(2,2);
        m.setZero();
        ASSERT_DOUBLE_EQ(0.0, m(0,0).real());
        ASSERT_DOUBLE_EQ(0.0, m(0,1).real());
        ASSERT_DOUBLE_EQ(0.0, m(1,0).real());
        ASSERT_DOUBLE_EQ(0.0, m(1,1).real());
        ASSERT_DOUBLE_EQ(0.0, m(0,0).imag());
        ASSERT_DOUBLE_EQ(0.0, m(0,1).imag());
        ASSERT_DOUBLE_EQ(0.0, m(1,0).imag());
        ASSERT_DOUBLE_EQ(0.0, m(1,1).imag());
    }
}
