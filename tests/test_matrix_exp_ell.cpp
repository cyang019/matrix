#include "Matrix.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>    // exp


namespace {
    using MatrixD = matrix::Matrix<double>; 
    using MatrixCx = matrix::Matrix<std::complex<double>>; 

    TEST(TestMatrix, ellDouble){
      MatrixD mat = matrix::ones<double>(4, 4);
      std::cout << "matrix:\n" << mat << "\n";

      int s = matrix::exponential::ell(mat, 3);
      std::cout << "s of ell(mat,3): " << s << "\n";

      s = matrix::exponential::ell(mat, 5);
      std::cout << "s of ell(mat,5): " << s << "\n";

      s = matrix::exponential::ell(mat, 7);
      std::cout << "s of ell(mat,7): " << s << "\n";

      s = matrix::exponential::ell(mat, 9);
      std::cout << "s of ell(mat,9): " << s << "\n";

      s = matrix::exponential::ell(mat, 11);
      std::cout << "s of ell(mat,11): " << s << "\n";

      s = matrix::exponential::ell(mat, 13);
      std::cout << "s of ell(mat,13): " << s << "\n";
    }
}

