#ifndef _MATRIX_BLAS_LEVEL3_H
#define _MATRIX_BLAS_LEVEL3_H

#include "blas_wrapper/cblas_common.h"
#include "common.h"


namespace matrix {
  inline namespace v1 {
    void lvl3_cblas_dgemm(
        const CblasOrder &layout,
        const CblasTranspose &TransA,
        const CblasTranspose &TransB,
        size_t m, size_t n, size_t k,
        double alpha, const double *A, int lda, const double *B, int ldb,
        double beta, double *C, int ldc);

    void lvl3_cblas_zgemm(
        const CblasOrder &layout,
        const CblasTranspose &TransA,
        const CblasTranspose &TransB,
        size_t m, size_t n, size_t k,
        cxdbl alpha, const cxdbl *A, int lda, const cxdbl *B, int ldb,
        cxdbl beta, cxdbl *C, int ldc);
  } // namespace v1
} // namespace matrix

#endif

