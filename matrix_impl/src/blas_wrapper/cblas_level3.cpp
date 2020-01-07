#include "configure_matrix.h"
#include "blas_wrapper/cblas_level3.h"
#include "matrix_core/errors.h"
#include <complex>

#include "cblas.h"

#ifdef HAVE_APPLE_LAPACK
  #include "clapack.h"
#elif defined HAVE_LAPACKE
  #include "lapacke.h"
#else
#endif


using namespace std;


namespace matrix {
  inline namespace v1 {
    void lvl3_dgemm(
        const CblasOrder &layout,
        const CblasTranspose &TransA,
        const CblasTranspose &TransB,
        size_t m, size_t n, size_t k,
        double alpha, const double *A, int lda, const double *B, int ldb,
        double beta, double *C, int ldc)
    {
#ifndef NDEBUG
      if (m >= (size_t)int_max || n >= (size_t)int_max || k >= (size_t)int_max) {
        throw IndexOutOfBound("Matrices dimensions need to be smaller than INT_MAX.");
      }
#endif
      cblas_dgemm((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
        (int)m, (int)n, (int)k,
        alpha, A, lda, B, ldb, beta, C, ldc);
    }

    void lvl3_zgemm(
        const CblasOrder &layout,
        const CblasTranspose &TransA,
        const CblasTranspose &TransB,
        size_t m, size_t n, size_t k,
        cxdbl alpha, const cxdbl *A, int lda, const cxdbl *B, int ldb,
        cxdbl beta, cxdbl *C, int ldc)
    {
#ifndef NDEBUG
      if (m >= (size_t)int_max || n >= (size_t)int_max || k >= (size_t)int_max) {
        throw IndexOutOfBound("Matrices dimensions need to be smaller than INT_MAX.");
      }
#endif
      cblas_zgemm((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
          (int)m, (int)n, (int)k,
          &alpha, A, lda, B, ldb, &beta, C, ldc);
    }

    int mat_dgetrf(size_t m, size_t n, double *A, size_t lda, int *ipiv, int *info)
    {
#ifndef NDEBUG
      if (m >= (size_t)int_max || n >= (size_t)int_max || lda >= (size_t)int_max) {
        throw IndexOutOfBound("Matrices dimensions need to be smaller than INT_MAX.");
      }
#endif

      int res = 0;
#ifdef HAVE_APPLE_LAPACK
      int i_m = (int)m;
      int i_n = (int)n;
      int i_lda = (int)lda;
<<<<<<< HEAD
#ifdef HAVE_APPLE_LAPACK
      int res = dgetrf_(&i_m, &i_n, A, &i_lda, ipiv, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, i_m, i_n, A, i_lda, ipiv);
      int res = *info;
#else
=======
      res = dgetrf_(&i_m, &i_n, A, &i_lda, ipiv, info);
>>>>>>> bfea27faacc032e19f2a464d8791b07a71ad6670
#endif
      return res;
    }

    int mat_dgetri(size_t n, double *A, size_t lda, int * ipiv, 
        double *work, size_t lwork, int *info)
    {
#ifndef NDEBUG
      if (n >= (size_t)int_max || lda >= (size_t)int_max || lwork >= (size_t)int_max) {
        throw IndexOutOfBound("Matrices dimensions need to be smaller than INT_MAX.");
      }
#endif
      int res = 0;
#ifdef HAVE_APPLE_LAPACK
      int i_n = (int)n;
      int i_lda = (int)lda;
#ifdef HAVE_APPLE_LAPACK
      int i_lwork = (int)lwork;
      int res = dgetri_(&i_n, A, &i_lda, ipiv, work, &i_lwork, info);
#elif defined HAVE_CLAPACK
      int i_lwork = (int)lwork;
<<<<<<< HEAD
      int res = dgetri_(&i_n, A, &i_lda, ipiv, work, &i_lwork, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_dgetri(LAPACK_COL_MAJOR, i_n, A, i_lda, ipiv);
      int res = *info;
#else
=======
      res = dgetri_(&i_n, A, &i_lda, ipiv, work, &i_lwork, info);
>>>>>>> bfea27faacc032e19f2a464d8791b07a71ad6670
#endif
      return res;
    }
  } // namespace v1
} // namespace matrix
