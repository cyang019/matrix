#include "matrix_core/Matrix.h"
#include "blas_wrapper/linsol.h"
#include "clapack.h"
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
    std::pair<Matrix<cxdbl>, int> lstsq(const Matrix<cxdbl> &a, const Matrix<cxdbl> &b,
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

      auto s = make_unique<double[]>(std::min(m,n));
      int rank = 0;
      int lwork = 1;
      auto work = make_unique<cxdbl[]>(lwork);
      auto rwork = make_unique<double[]>(1);
      auto iwork = make_unique<int[]>(1);
      int info = 0;

      // query sizes
#ifndef NDEBUG
      auto t1 = std::clock();
#endif
      int res = mat_zgelsd(m, n,nrhs, 
          a_local.data(), lda, b_mat.get(), ldb, s.get(), &rcond,
          &rank, work.get(), -1, rwork.get(), iwork.get(), &info);
#ifndef NDEBUG
      auto t2 = std::clock();
      cout << "Query cpu time: " << 1000.0 * (t2 - t1)/CLOCKS_PER_SEC << endl;
#endif

      lwork = static_cast<int>(work[0].real());
      int lrwork = static_cast<int>(rwork[0]);
      int liwork = static_cast<int>(iwork[0]);
#ifndef NDEBUG
      cout << "lwork: " << lwork << "\t"
           << "lrwork: " << lrwork << "\t"
           << "lrwork: " << liwork << endl;
#endif
      
      work = make_unique<cxdbl[]>(lwork);
      rwork = make_unique<double[]>(lrwork);
      iwork = make_unique<int[]>(liwork);

#ifndef NDEBUG
      auto t3 = std::clock();
#endif
      res = mat_zgelsd(m, n, nrhs, 
          a_local.data(), lda, b_mat.get(), ldb, s.get(), &rcond,
          &rank, work.get(), lwork, rwork.get(), iwork.get(), &info);
#ifndef NDEBUG
      auto t4 = std::clock();
      cout << "compute cpu time: " << 1000.0 * (t4 - t3)/CLOCKS_PER_SEC << endl;
#endif
      Matrix<cxdbl> x(n, nrhs);
      for(size_t i = 0; i < nrhs; ++i){
        for(size_t j = 0; j < n; ++j){
          x(j, i) = b_mat[j + i * m];
        }
      }
      return make_pair(std::move(x), info);
    }

  } // namespace v1
} // namespace matrix
