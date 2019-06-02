#include "Matrix.h"
#include "gtest/gtest.h"
#include <complex>
#include <iostream>

namespace {
    using MatrixCD = matrix::Matrix<std::complex<double>>; 
    using cxdbl = std::complex<double>;

    TEST(TestMatrix, ZheevEigenVals){
      MatrixCD m = {{cxdbl(0,0), cxdbl(0,-0.5)},
                    {cxdbl(0,0.5), cxdbl(0,0)}};
      auto vals = matrix::eigenVal<matrix::EigenMethod::zheev>(m);
      ASSERT_EQ(2, vals.nrows());
      ASSERT_DOUBLE_EQ(-0.5, vals(0,0));
      ASSERT_DOUBLE_EQ(0.5, vals(1,0));
    }
}

