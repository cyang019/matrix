#include "cxdbl_mat_op.h"
#include <chrono>
#include <iostream>

using namespace std;

void mat_mat_mul(size_t n, size_t cnt)
{
  using cxdbl = std::complex<double>;
  matrix::Matrix<cxdbl> mat1(n, n);
  matrix::Matrix<cxdbl> mat2(n, n);
  mat1.setRandom();
  mat2.setRandom();


  cout << "matrix multiplication [" << n << " x " << n << "] repeat "
       << cnt << " times took: \n";
  auto t1 = std::chrono::high_resolution_clock::now();
  for(size_t i = 0; i < cnt; ++i) {
    auto res = mat1 * mat2;
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  cout << duration << " microseconds." << endl;
}
