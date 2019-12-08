#include "configure_matrix.h"
#ifdef HAVE_APPLE_LAPACK
#include "clapack.h"
#elif defined HAVE_CLAPACK
#include "clapack.h"
#elif defined HAVE_LAPACKE
#include "lapacke.h"
#else
#endif

#include "blas_wrapper/linsol.h"
#include "matrix_core/errors.h"


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

      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_b = reinterpret_cast<__CLPK_doublecomplex *>(b);
      int res = zgesv_(&dim1, &dim2, clpk_a, &lower_dim_a,
          ipiv, clpk_b, &lower_dim_b, info);
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

      int res = dgesv_(&dim1, &dim2, a, &lower_dim_a,
          ipiv, b, &lower_dim_b, info);
      return res;
    }

    int mat_zgelsd(
        size_t m, size_t n, size_t nrhs,
        cxdbl *a, size_t lda, cxdbl *b, size_t ldb,
        double *s, double *rcond,
        int *rank, cxdbl *work,
        int lwork, double *rwork,
        int *iwork, int *info)
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
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_b = reinterpret_cast<__CLPK_doublecomplex *>(b);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work);

      int res = zgelsd_(&nrows, &ncols, &ncols_b,
          clpk_a, &lower_dim_a, clpk_b, &lower_dim_b,
          s, rcond, rank, clpk_work, &lwork,
          rwork, iwork, info);
      return res;
    }
  } // namespace v1
} // namespace matrix
