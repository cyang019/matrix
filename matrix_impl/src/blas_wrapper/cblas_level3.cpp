#include "configure_matrix.h"
#include "blas_wrapper/cblas_level3.h"
#include "matrix_core/errors.h"
#include <complex>
#include <memory>

#ifdef MSVC
#include "mkl_cblas.h"
#else
#include "cblas.h"
#endif

#ifdef MSVC
#include "mkl_lapacke.h"
#elif defined(HAVE_APPLE_LAPACK) || defined(HAVE_CLAPACK)
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
      int i_m = (int)m;
      int i_n = (int)n;
      int i_lda = (int)lda;
#if defined(HAVE_APPLE_LAPACK) || defined(HAVE_CLAPACK)
      res = dgetrf_(&i_m, &i_n, A, &i_lda, ipiv, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, i_m, i_n, A, i_lda, ipiv);
      res = *info;
#endif
      return res;
    }

    int mat_dgetri(size_t n, double *A, size_t lda, int * ipiv, 
        int *info)
    {
#ifndef NDEBUG
      if (n >= (size_t)int_max || lda >= (size_t)int_max) {
        throw IndexOutOfBound("Matrices dimensions need to be smaller than INT_MAX.");
      }
#endif
      int res = 0;
      int i_n = (int)n;
      int i_lda = (int)lda;
#if defined(HAVE_APPLE_LAPACK) || defined(HAVE_CLAPACK)
      // query lwork
      int lwork = -1;
      auto work = std::make_unique<double[]>(1);
      res = dgetri_(&i_n, A, &i_lda, ipiv, work.get(), &lwork, info);
      lwork = static_cast<int>(work[0]);
      work = std::make_unique<double[]>(lwork);
      // actual calculation
      // inv(A)*L = inv(U)
      res = dgetri_(&i_n, A, &i_lda, ipiv, work.get(), &lwork, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_dgetri(LAPACK_COL_MAJOR, i_n, A, i_lda, ipiv);
      res = *info;
#endif
      return res;
    }

    int mat_zgetrf(size_t m, size_t n,
        cxdbl *A, size_t lda,
        int *ipiv, int *info)
    {
#ifndef NDEBUG
      if (m >= (size_t)int_max || n >= (size_t)int_max || lda >= (size_t)int_max) {
        throw IndexOutOfBound("Matrices dimensions need to be smaller than INT_MAX.");
      }
#endif

      int res = 0;
      int i_m = (int)m;
      int i_n = (int)n;
      int i_lda = (int)lda;
#if defined(HAVE_APPLE_LAPACK) || defined(HAVE_CLAPACK)
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(A);
      res = zgetrf_(&i_m, &i_n, clpk_a, &i_lda, ipiv, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, i_m, i_n, A, i_lda, ipiv);
      res = *info;
#endif
      return res;
    }

    int mat_zgetri(size_t n, cxdbl *A, size_t lda, int * ipiv, 
        int *info)
    {
#ifndef NDEBUG
      if (n >= (size_t)int_max || lda >= (size_t)int_max) {
        throw IndexOutOfBound("Matrices dimensions need to be smaller than INT_MAX.");
      }
#endif
      int res = 0;
      int i_n = (int)n;
      int i_lda = (int)lda;
#if defined(HAVE_APPLE_LAPACK) || defined(HAVE_CLAPACK)
      // query lwork
      int lwork = -1;
      auto work = std::make_unique<cxdbl[]>(1);
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(A);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
      res = zgetri_(&i_n, clpk_a, &i_lda, ipiv, clpk_work, &lwork, info);
      lwork = static_cast<int>(work[0].real());
      work = std::make_unique<cxdbl[]>(lwork);
      clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
      // actual calculation
      // inv(A)*L = inv(U)
      res = zgetri_(&i_n, clpk_a, &i_lda, ipiv, clpk_work, &lwork, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_zgetri(LAPACK_COL_MAJOR, i_n, A, i_lda, ipiv);
      res = *info;
#endif
      return res;
    }
  } // namespace v1
} // namespace matrix
