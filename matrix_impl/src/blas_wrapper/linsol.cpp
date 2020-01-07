#include "configure_matrix.h"
#ifdef HAVE_APPLE_LAPACK
  #include "clapack.h"
#elif defined HAVE_LAPACKE
  #include "lapacke.h"
#else
#endif
#include "blas_wrapper/linsol.h"
#include "matrix_core/common.h"
#include "matrix_core/errors.h"
#include <memory>
#include <chrono>
#include <ctime>
#include <iostream>


namespace matrix {
  inline namespace v1 {
    int mat_zgesv(
        size_t n, size_t nrhs, cxdbl *a, size_t lda,
        int *ipiv, cxdbl *b, size_t ldb, int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max || ldb > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      int dim1 = (int)n;
      int dim2 = (int)nrhs;
      int lower_dim_a = (lda > 0) ? (int)lda : 1;
      int lower_dim_b = (ldb > 0) ? (int)ldb : 1;

      int res = 0;
#ifdef HAVE_APPLE_LAPACK
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_b = reinterpret_cast<__CLPK_doublecomplex *>(b);
      res = zgesv_(&dim1, &dim2, clpk_a, &lower_dim_a,
          ipiv, clpk_b, &lower_dim_b, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_zgesv(LAPACK_COL_MAJOR,
          dim1, dim2, a, lower_dim_a,
          ipiv, b, lower_dim_b);
      res = *info;
#endif
      return res;
    }

    int mat_dgesv(
        size_t n, size_t nrhs, double *a, size_t lda,
        int *ipiv, double *b, size_t ldb, int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max || ldb > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      int dim1 = (int)n;
      int dim2 = (int)nrhs;
      int lower_dim_a = (lda > 0) ? (int)lda : 1;
      int lower_dim_b = (ldb > 0) ? (int)ldb : 1;

      int res = 0;
#ifdef HAVE_APPLE_LAPACK
      res = dgesv_(&dim1, &dim2, a, &lower_dim_a,
          ipiv, b, &lower_dim_b, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_dgesv(LAPACK_COL_MAJOR, dim1, dim2, 
          a, lower_dim_a,
          ipiv, b, lower_dim_b);
      res = *info;
#endif
      return res;
    }

    int mat_zgelsd(
        size_t m, size_t n, size_t nrhs,
        cxdbl *a, size_t lda, cxdbl *b, size_t ldb,
        double *s, double *rcond,
        int *rank, int *info)
    {
#ifndef NDEBUG
      if(m > (size_t)int_max 
          || n > (size_t)int_max
          || nrhs > (size_t)int_max
          || lda > (size_t)int_max 
          || ldb > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      int nrows = (int)m;
      int ncols = (int)n;
      int ncols_b = (int)nrhs;
      int lower_dim_a = (lda > 0) ? (int)lda : 1;
      int lower_dim_b = (ldb > 0) ? (int)ldb : 1;

      int res = 0;
#ifdef HAVE_APPLE_LAPACK
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_b = reinterpret_cast<__CLPK_doublecomplex *>(b);

      int lwork = 1;
      auto work = std::make_unique<cxdbl[]>(lwork);
      auto rwork = std::make_unique<double[]>(1);
      auto iwork = std::make_unique<int[]>(1);

      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
#ifdef VERBOSE
      auto t1 = std::clock();
#endif
      // query size first
      res = zgelsd_(&nrows, &ncols, &ncols_b,
          clpk_a, &lower_dim_a, clpk_b, &lower_dim_b,
          s, rcond, rank, clpk_work, -1,
          rwork.get(), iwork.get(), info);
#ifdef VERBOSE
      auto t2 = std::clock();
      std::cout << "Query cpu time: " 
                << 1000.0 * (t2 - t1)/CLOCKS_PER_SEC << std::endl;
#endif

      lwork = static_cast<int>(work[0].real());
      int lrwork = static_cast<int>(rwork[0]);
      int liwork = static_cast<int>(iwork[0]);
#ifdef VERBOSE
      cout << "lwork: " << lwork << "\t"
           << "lrwork: " << lrwork << "\t"
           << "lrwork: " << liwork << endl;
#endif
      work = std::make_unique<cxdbl[]>(lwork);
      rwork = std::make_unique<double[]>(lrwork);
      iwork = std::make_unique<int[]>(liwork);
      clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
#ifdef VERBOSE
      auto t3 = std::clock();
#endif
      // the real calculation
      res = zgelsd_(&nrows, &ncols, &ncols_b,
          clpk_a, &lower_dim_a, clpk_b, &lower_dim_b,
          s, rcond, rank, clpk_work, &lwork,
          rwork.get(), iwork.get(), info);
#ifdef VERBOSE
      auto t4 = std::clock();
      cout << "compute cpu time: " << 1000.0 * (t4 - t3)/CLOCKS_PER_SEC << endl;
#endif
#elif defined HAVE_LAPACKE
      *info = LAPACKE_zgelsd(LAPACK_COL_MAJOR, nrows, ncols, 
          ncols_b, a, lower_dim_a, b, lower_dim_b,
          s, rcond, rank);
      res = *info;
#endif
      return res;
    }
  } // namespace v1
} // namespace matrix
