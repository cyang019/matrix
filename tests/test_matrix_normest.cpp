#include "Matrix.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>    // exp


namespace {
    using MatrixD = matrix::Matrix<double>; 
    using MatrixCx = matrix::Matrix<std::complex<double>>; 

    TEST(TestMatrix, normestDouble){
      constexpr auto evalMat = [](const MatrixD &mat){
        double norm = matrix::norm1(mat);
        double est = matrix::exponential::normest(mat, 1);

#ifdef VERBOSE
        std::cout << mat << "\n";
#endif
        std::cout << "1-norm: " << norm << "\n";
        std::cout << "1-norm estimation: " << est << "\n\n";
        ASSERT_TRUE(est <= norm);
      };

      MatrixD mat1 = matrix::ones<double>(4, 4);
      evalMat(mat1);

      MatrixD mat2 = matrix::diagonal({0.5, -0.5});
      evalMat(mat2);

      MatrixD mat3 = {{0, 0.5}, {0.5, 0}};
      evalMat(mat3);

      MatrixD mat4 = matrix::diagonal({-1.0, 0.0, 1.0});
      evalMat(mat4);

      MatrixD mat5 = {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}};
      evalMat(mat5);

      MatrixD mat6 = {{0, 1, 0}, {-1, 0, 1}, {0, -1, 0}};
      evalMat(mat6);

      MatrixD mat7 = {
        {0, 2, 0, 0, 0},
        {2, 0, std::sqrt(6.0), 0, 0},
        {0, std::sqrt(6.0), 0, std::sqrt(6.0), 0},
        {0, 0, std::sqrt(6.0), 0, 2.0},
        {0, 0, 0, 2.0, 0}
      };
      evalMat(mat7);

      MatrixD mat8 = matrix::diagonal({2.5, 1.5, 0.5, -0.5, -1.5, -2.5});
      evalMat(mat8);

      MatrixD mat9 = {
        {0, std::sqrt(5), 0, 0, 0, 0},
        {std::sqrt(5), 0, std::sqrt(8), 0, 0, 0},
        {0, std::sqrt(8), 0, 3, 0, 0},
        {0, 0, 3, 0, std::sqrt(8), 0},
        {0, 0, 0, std::sqrt(8), 0, std::sqrt(5)},
        {0, 0, 0, 0, std::sqrt(5), 0}
      };
      evalMat(mat9);

      MatrixD mat10(16, 16);
      for(size_t i = 0; i < 5; ++i){
        mat10.setRandom();
        evalMat(mat10);
      }
    }

    TEST(TestMatrix, normestCxDbl){
      constexpr auto evalMat = [](const MatrixCx &mat){
        double norm = matrix::norm1(mat);
        double est = matrix::exponential::normest(mat, 1);

#ifdef VERBOSE
        std::cout << mat << "\n";
#endif
        std::cout << "1-norm: " << norm << "\n";
        std::cout << "1-norm estimation: " << est << "\n\n";
        ASSERT_TRUE(est <= norm);
      };

      MatrixCx mat1 = matrix::ones<std::complex<double>>(4, 4);
      evalMat(mat1);

      MatrixCx mat2 = matrix::diagonal({0.5, -0.5});
      evalMat(mat2);

      MatrixCx mat3 = {{0, 0.5}, {0.5, 0}};
      evalMat(mat3);

      MatrixCx mat4 = matrix::diagonal({-1.0, 0.0, 1.0});
      evalMat(mat4);

      MatrixCx mat5 = {{0.0, 1.0, 0.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 0.0}};
      evalMat(mat5);

      MatrixCx mat6 = {{0.0, 1.0, 0.0}, {-1.0, 0.0, 1.0}, {0.0, -1.0, 0.0}};
      evalMat(mat6);

      MatrixCx mat7 = {
        {0, 2, 0, 0, 0},
        {2, 0, std::sqrt(6.0), 0, 0},
        {0, std::sqrt(6.0), 0, std::sqrt(6.0), 0},
        {0, 0, std::sqrt(6.0), 0, 2.0},
        {0, 0, 0, 2.0, 0}
      };
      evalMat(mat7);

      MatrixCx mat8 = matrix::diagonal({2.5, 1.5, 0.5, -0.5, -1.5, -2.5});
      evalMat(mat8);

      MatrixCx mat9 = {
        {0, std::sqrt(5), 0, 0, 0, 0},
        {std::sqrt(5), 0, std::sqrt(8), 0, 0, 0},
        {0, std::sqrt(8), 0, 3, 0, 0},
        {0, 0, 3, 0, std::sqrt(8), 0},
        {0, 0, 0, std::sqrt(8), 0, std::sqrt(5)},
        {0, 0, 0, 0, std::sqrt(5), 0}
      };
      evalMat(mat9);

      MatrixCx mat10(16, 16);
      for(size_t i = 0; i < 5; ++i){
        mat10.setRandom();
        evalMat(mat10);
      }
    }
}


