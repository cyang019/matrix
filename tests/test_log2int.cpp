#include "matrix.h"
#include "gtest/gtest.h"
#include <complex>
#include <iostream>

namespace {
    using namespace std;
    using matrix::log2int;

    TEST(TestMatrix, Log2int){
      ASSERT_EQ(3, log2int(4));
      ASSERT_EQ(3, log2int(5));
      ASSERT_EQ(4, log2int(8));
      ASSERT_EQ(4, log2int(9));
      ASSERT_EQ(5, log2int(16));
      ASSERT_EQ(6, log2int(33));
    }
}
