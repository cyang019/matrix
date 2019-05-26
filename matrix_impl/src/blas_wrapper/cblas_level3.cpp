#include "cblas.h"
#include "blas_wrapper/cblas_level3.h"
#include "errors.h"


namespace matrix {
  inline namespace v1 {
    void lvl3_cblas_dgemm(
        const CblasOrder &layout,
        const CblasTranspose &TransA,
        const CblasTranspose &TransB,
        size_t m, size_t n, size_t k,
        double alpha, const double *A, int lda, const double *B, int ldb,
        double beta, double *C, int ldc)
    {
#ifndef NDEBUG
      if (m >= int_max || n >= int_max || k >= int_max) {
        throw IndexOutOfBound("Matrices dimensions need to be smaller than INT_MAX.");
      }
#endif
      cblas_dgemm((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
        (int)m, (int)n, (int)k,
        alpha, A, lda, B, ldb, beta, C, ldc);
    }

    void lvl3_cblas_zgemm(
        const CblasOrder &layout,
        const CblasTranspose &TransA,
        const CblasTranspose &TransB,
        size_t m, size_t n, size_t k,
        cxdbl alpha, const cxdbl *A, int lda, const cxdbl *B, int ldb,
        cxdbl beta, cxdbl *C, int ldc)
    {
#ifndef NDEBUG
      if (m >= int_max || n >= int_max || k >= int_max) {
        throw IndexOutOfBound("Matrices dimensions need to be smaller than INT_MAX.");
      }
#endif
      cblas_zgemm((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
          (int)m, (int)n, (int)k,
          &alpha, A, lda, B, ldb, &beta, C, ldc);
    }
  } // namespace v1
} // namespace matrix
