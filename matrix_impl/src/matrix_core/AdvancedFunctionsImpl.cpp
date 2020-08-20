#include "configure_matrix.h"
#include "matrix_core/Matrix.h"
#include "blas_wrapper/linsol.h"
#include <memory>
#include <utility>
#include <cmath>
#include <iostream>
#include <chrono>
#include <ctime>

using namespace std;

namespace matrix {
  inline namespace v1 {
    // solve for A*x = B by minimizing 2-norm(abs(b - A*x))
    std::pair<Matrix<cxdbl>, int> lstsq(
        const Matrix<cxdbl> &a, const Matrix<cxdbl> &b,
        double rcond)
    {
#ifndef NDEBUG
      if(a.nrows() != b.nrows()){
         throw MatrixSizeMismatchError("A *  X = B A.nrows() needs to be the same as B.nrows()");
      }
#endif
      const size_t m = a.nrows();
      const size_t n = a.ncols();     
      const size_t nrhs = b.ncols();
      size_t lda =  m > 1 ? m : 1;
      size_t ldb = n > lda ? n : lda;
      Matrix<cxdbl> a_local = a;
      auto b_mat = make_unique<cxdbl[]>(ldb * nrhs);
      // b was m x nrhs
      for(size_t i = 0; i < nrhs; ++i){
        for(size_t j = 0; j < m; ++j){
          size_t idx = i * m + j;
          b_mat[idx] = b(j, i);
        }
      }

      auto jpvt = make_unique<int[]>(n);
      for(size_t i = 0; i < n; ++i){
        jpvt[i] = 0;
      }
      int rank = 0;
      int info = 0;

      // least square problem
      int status = mat_zgelsy(m, n, nrhs, 
          a_local.data(), lda, b_mat.get(), ldb, jpvt.get(), &rcond,
          &rank, &info);

      Matrix<cxdbl> x(n, nrhs);
      for(size_t i = 0; i < nrhs; ++i){
        for(size_t j = 0; j < n; ++j){
          x(j, i) = b_mat[j + i * m];
        }
      }
      status = status || info;
      return make_pair(std::move(x), status);
    }
  } // namespace v1
} // namespace matrix

