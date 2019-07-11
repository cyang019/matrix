#include "matrix.h"
#include "gtest/gtest.h"
#include <complex>

namespace {
    using MatrixCD = matrix::Matrix<std::complex<double>>; 
    using cxdbl = std::complex<double>;

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

    TEST(TestMatrix, CastCLPKCxDbl){
      MatrixCD m = {{cxdbl(0,0), cxdbl(0,-1)},
                    {cxdbl(0, 1), cxdbl(0,0)}};

      cxdbl *ptr_data = m.data();
      ASSERT_DOUBLE_EQ(0.0, ptr_data->real());
      ASSERT_DOUBLE_EQ(0.0, ptr_data->imag());
      ++ptr_data;
      ASSERT_DOUBLE_EQ(0.0, ptr_data->real());
      ASSERT_DOUBLE_EQ(1.0, ptr_data->imag());
      ++ptr_data;
      ASSERT_DOUBLE_EQ(0.0, ptr_data->real());
      ASSERT_DOUBLE_EQ(-1.0, ptr_data->imag());
      ++ptr_data;
      ASSERT_DOUBLE_EQ(0.0, ptr_data->real());
      ASSERT_DOUBLE_EQ(0.0, ptr_data->imag());

      double *ptr_casted = reinterpret_cast<double *>(m.data());
      ASSERT_DOUBLE_EQ(0.0, ptr_casted[0]);
      ASSERT_DOUBLE_EQ(0.0, ptr_casted[1]);
      ASSERT_DOUBLE_EQ(0.0, ptr_casted[2]);
      ASSERT_DOUBLE_EQ(1.0, ptr_casted[3]);
      ASSERT_DOUBLE_EQ(0.0, ptr_casted[4]);
      ASSERT_DOUBLE_EQ(-1.0, ptr_casted[5]);
      ASSERT_DOUBLE_EQ(0.0, ptr_casted[6]);
      ASSERT_DOUBLE_EQ(0.0, ptr_casted[7]);

      struct my_cx { double r, i; };
      my_cx *ptr_casted2 = reinterpret_cast<my_cx *>(m.data());
      ASSERT_DOUBLE_EQ(0.0,  ptr_casted2[0].r);
      ASSERT_DOUBLE_EQ(0.0,  ptr_casted2[0].i);
      ASSERT_DOUBLE_EQ(0.0,  ptr_casted2[1].r);
      ASSERT_DOUBLE_EQ(1.0,  ptr_casted2[1].i);
      ASSERT_DOUBLE_EQ(0.0,  ptr_casted2[2].r);
      ASSERT_DOUBLE_EQ(-1.0, ptr_casted2[2].i);
      ASSERT_DOUBLE_EQ(0.0,  ptr_casted2[3].r);
      ASSERT_DOUBLE_EQ(0.0,  ptr_casted2[3].i);

      MatrixCD m_2 = {{cxdbl(1, 2), cxdbl(3,-4)},
                     {cxdbl(5, 6), cxdbl(7, 8)}};
      my_cx *ptr_casted_2 = reinterpret_cast<my_cx *>(m_2.data());
      ASSERT_DOUBLE_EQ(1.0,  ptr_casted_2[0].r);
      ASSERT_DOUBLE_EQ(2.0,  ptr_casted_2[0].i);
      ASSERT_DOUBLE_EQ(5.0,  ptr_casted_2[1].r);
      ASSERT_DOUBLE_EQ(6.0,  ptr_casted_2[1].i);
      ASSERT_DOUBLE_EQ(3.0,  ptr_casted_2[2].r);
      ASSERT_DOUBLE_EQ(-4.0, ptr_casted_2[2].i);
      ASSERT_DOUBLE_EQ(7.0,  ptr_casted_2[3].r);
      ASSERT_DOUBLE_EQ(8.0,  ptr_casted_2[3].i);
    }
}
