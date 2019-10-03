#include "matrix.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>    // exp


namespace {
    using MatrixD = matrix::Matrix<double>; 
    using MatrixCD = matrix::Matrix<std::complex<double>>; 
    using MatrixI = matrix::Matrix<int>; 

    TEST(TestMatrix, areParallel){
      MatrixD mat = {
        {1, 2, -3, 4},
        {-1, 2, 3, 4},
        {1, 2, -3, 4}
      };
      
      ASSERT_EQ(true, matrix::exponential::areParallel(mat, 0, mat, 2));
      ASSERT_EQ(true, matrix::exponential::areParallel(mat, 1, mat, 3));
      ASSERT_EQ(false, matrix::exponential::areParallel(mat, 0, mat, 1));
    }

    TEST(TestMatrix, DoubleExp){
      MatrixD m1= matrix::diagonal({1.0,2.0,3.0,4.0});
      MatrixD res1 = matrix::exp(m1);
      MatrixD desired = matrix::diagonal(
          { std::exp(1.0), std::exp(2.0), std::exp(3.0), std::exp(4.0)});
      std::cout << std::setprecision(15);
      std::cout << "desired:\n" << desired << "\n";
      std::cout << "result:\n" << res1 << "\n";

      double eps = matrix::norm1(desired) * matrix::eps * 5;
      std::cout << "desired matrix norm1: " << matrix::norm1(desired) << "\n";

      ASSERT_TRUE(matrix::allclose(desired, res1, eps));

      MatrixCD m2 = {{1.0, std::complex<double>(0.0, 2.0), 10.0}, 
                     {std::complex<double>(0.0, -2.0), 0.0, 0.3}, 
                     {10.0, 0.3, -0.1}};
      auto [eigen_vals, eigen_vecs] = matrix::eigenSys(m2);
      MatrixCD m2_desired = matrix::zeros<std::complex<double>>(m2.nrows(), m2.ncols());
      for(size_t i = 0; i < m2.nrows(); ++i){
        m2_desired(i,i) = std::exp(eigen_vals(i,0));
      }
      m2_desired = eigen_vecs.adjoint() * m2_desired * eigen_vecs;
      auto m2_exp = matrix::exp(m2);
      std::cout << "exp by expansion:\n" << m2_exp << "\n";
      std::cout << "exp by eigen:\n" << m2_desired << std::endl;

      MatrixCD m3 = {{0.0, std::complex<double>(0.0, -0.5)},
                     {std::complex<double>(0.0, 0.5), 0.0}};
      auto [eigen_vals3, eigen_vecs3] = matrix::eigenSys(m3);
      MatrixCD m3_desired = matrix::zeros<std::complex<double>>(m3.nrows(), m3.ncols());
      for(size_t i = 0; i < m3.nrows(); ++i){
        m3_desired(i,i) = std::exp(eigen_vals(i,0));
      }
      m3_desired = eigen_vecs3.adjoint() * m3_desired * eigen_vecs3;
      auto m3_exp = matrix::exp(m3);
      std::cout << "exp by expansion:\n" << m3_exp << "\n";
      std::cout << "exp by eigen:\n" << m3_desired << std::endl;
    }

    TEST(TestMatrix, ExpMinusOne){
      MatrixCD m = {{-0.5, 0}, {0, 0.5}};
      MatrixCD desired = {{std::exp(-0.5) - 1.0, 0.0}, {0.0, std::exp(0.5) - 1.0}};
      MatrixCD res = matrix::expMinusIdentity(m);

      double eps = matrix::norm1(desired) * matrix::eps * 5;
      std::cout << "desired matrix norm1: " << matrix::norm1(desired) << "\n";
      ASSERT_TRUE(matrix::allclose(desired, res, eps));
    }
}

