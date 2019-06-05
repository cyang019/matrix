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

      MatrixCD m1 = matrix::diagnal({cxdbl(1.0,0), cxdbl(2.0,0), cxdbl(3.0,0)});
      auto vals1 = matrix::eigenVal<matrix::EigenMethod::zheev>(m1);
      ASSERT_EQ(3, vals1.nrows());
      ASSERT_DOUBLE_EQ(1.0, vals1(0,0));
      ASSERT_DOUBLE_EQ(2.0, vals1(1,0));
      ASSERT_DOUBLE_EQ(3.0, vals1(2,0));

      MatrixCD m2 = {{cxdbl(2,0), cxdbl(2,1), cxdbl(4,0)},
                     {cxdbl(2,-1), cxdbl(3,0), cxdbl(0,1)},
                     {cxdbl(4,0), cxdbl(0,-1), cxdbl(1,0)}};
      auto vals2 = matrix::eigenVal<matrix::EigenMethod::zheev>(m2);
      matrix::Matrix<double> vals2_desired = {
        {-3.1651250644074413},
        {2.8530792851693514},
        {6.312045779238087}};
      ASSERT_EQ(3, vals2.nrows());
      ASSERT_TRUE(matrix::allclose(vals2_desired, vals2, 1.0e-14));

      MatrixCD m3 = {{cxdbl(-1,0), cxdbl(1,-1), cxdbl(1,2), cxdbl(0,-1)},
                     {cxdbl(1,1), cxdbl(3,0), cxdbl(-2,0), cxdbl(3,-2)},
                     {cxdbl(1,-2), cxdbl(-2,0), cxdbl(0,0), cxdbl(4,0)},
                     {cxdbl(0,1), cxdbl(3,2), cxdbl(4,0), cxdbl(2,0)}};
      matrix::Matrix<double> vals3_desired = {
        {-5.567011257613931},
        {-1.336173855519022},
        {4.076219465400827},
        {6.826965647732123}
      };
      auto vals3 = matrix::eigenVal<matrix::EigenMethod::zheev>(m3);
      ASSERT_EQ(4, vals3.nrows());
      ASSERT_TRUE(matrix::allclose(vals3_desired, vals3, 1.0e-14));
    }

    TEST(TestMatrix, ZheevdEigenVals){
      MatrixCD m = {{cxdbl(0,0), cxdbl(0,-0.5)},
                    {cxdbl(0,0.5), cxdbl(0,0)}};
      auto vals = matrix::eigenVal<matrix::EigenMethod::zheevd>(m);
      ASSERT_EQ(2, vals.nrows());
      ASSERT_DOUBLE_EQ(-0.5, vals(0,0));
      ASSERT_DOUBLE_EQ(0.5, vals(1,0));

      MatrixCD m1 = matrix::diagnal({cxdbl(1.0,0), cxdbl(2.0,0), cxdbl(3.0,0)});
      auto vals1 = matrix::eigenVal<matrix::EigenMethod::zheevd>(m1);
      ASSERT_EQ(3, vals1.nrows());
      ASSERT_DOUBLE_EQ(1.0, vals1(0,0));
      ASSERT_DOUBLE_EQ(2.0, vals1(1,0));
      ASSERT_DOUBLE_EQ(3.0, vals1(2,0));

      MatrixCD m2 = {{cxdbl(2,0), cxdbl(2,1), cxdbl(4,0)},
                     {cxdbl(2,-1), cxdbl(3,0), cxdbl(0,1)},
                     {cxdbl(4,0), cxdbl(0,-1), cxdbl(1,0)}};
      auto vals2 = matrix::eigenVal<matrix::EigenMethod::zheevd>(m2);
      matrix::Matrix<double> vals2_desired = {
        {-3.1651250644074413},
        {2.8530792851693514},
        {6.312045779238087}};
      ASSERT_EQ(3, vals2.nrows());
      ASSERT_TRUE(matrix::allclose(vals2_desired, vals2, 1.0e-14));

      MatrixCD m3 = {{cxdbl(-1,0), cxdbl(1,-1), cxdbl(1,2), cxdbl(0,-1)},
                     {cxdbl(1,1), cxdbl(3,0), cxdbl(-2,0), cxdbl(3,-2)},
                     {cxdbl(1,-2), cxdbl(-2,0), cxdbl(0,0), cxdbl(4,0)},
                     {cxdbl(0,1), cxdbl(3,2), cxdbl(4,0), cxdbl(2,0)}};
      matrix::Matrix<double> vals3_desired = {
        {-5.567011257613931},
        {-1.336173855519022},
        {4.076219465400827},
        {6.826965647732123}
      };
      auto vals3 = matrix::eigenVal<matrix::EigenMethod::zheevd>(m3);
      ASSERT_EQ(4, vals3.nrows());
      ASSERT_TRUE(matrix::allclose(vals3_desired, vals3, 1.0e-14));
    }

    TEST(TestMatrix, ZheevrEigenVals){
      MatrixCD m = {{cxdbl(0,0), cxdbl(0,-0.5)},
                    {cxdbl(0,0.5), cxdbl(0,0)}};
      auto vals = matrix::eigenVal<matrix::EigenMethod::zheevr>(m);
      ASSERT_EQ(2, vals.nrows());
      ASSERT_DOUBLE_EQ(-0.5, vals(0,0));
      ASSERT_DOUBLE_EQ(0.5, vals(1,0));

      MatrixCD m1 = matrix::diagnal({cxdbl(1.0,0), cxdbl(2.0,0), cxdbl(3.0,0)});
      auto vals1 = matrix::eigenVal<matrix::EigenMethod::zheevr>(m1);
      ASSERT_EQ(3, vals1.nrows());
      ASSERT_DOUBLE_EQ(1.0, vals1(0,0));
      ASSERT_DOUBLE_EQ(2.0, vals1(1,0));
      ASSERT_DOUBLE_EQ(3.0, vals1(2,0));

      MatrixCD m2 = {{cxdbl(2,0), cxdbl(2,1), cxdbl(4,0)},
                     {cxdbl(2,-1), cxdbl(3,0), cxdbl(0,1)},
                     {cxdbl(4,0), cxdbl(0,-1), cxdbl(1,0)}};
      auto vals2 = matrix::eigenVal<matrix::EigenMethod::zheevr>(m2);
      matrix::Matrix<double> vals2_desired = {
        {-3.1651250644074413},
        {2.8530792851693514},
        {6.312045779238087}};
      ASSERT_EQ(3, vals2.nrows());
      ASSERT_TRUE(matrix::allclose(vals2_desired, vals2, 1.0e-14));

      MatrixCD m3 = {{cxdbl(-1,0), cxdbl(1,-1), cxdbl(1,2), cxdbl(0,-1)},
                     {cxdbl(1,1), cxdbl(3,0), cxdbl(-2,0), cxdbl(3,-2)},
                     {cxdbl(1,-2), cxdbl(-2,0), cxdbl(0,0), cxdbl(4,0)},
                     {cxdbl(0,1), cxdbl(3,2), cxdbl(4,0), cxdbl(2,0)}};
      matrix::Matrix<double> vals3_desired = {
        {-5.567011257613931},
        {-1.336173855519022},
        {4.076219465400827},
        {6.826965647732123}
      };
      auto vals3 = matrix::eigenVal<matrix::EigenMethod::zheevr>(m3);
      ASSERT_EQ(4, vals3.nrows());
      ASSERT_TRUE(matrix::allclose(vals3_desired, vals3, 1.0e-14));
    }
}

