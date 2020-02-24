#ifndef _MATRIX_BLAS_LEVEL2_H
#define _MATRIX_BLAS_LEVEL2_H

#include "blas_wrapper/cblas_common.h"
#include "matrix_core/common.h"


namespace matrix {
  inline namespace v1 {
    /// m & n cannot exceed int_max
    void lvl2_dgemv(
        const CblasOrder &layout, const CblasTranspose &TransA,
        size_t m, size_t n, double alpha, const double *A, size_t lda,
        const double *x, int incx, double beta, double *y, int incy);

    void lvl2_zgemv(
        const CblasOrder &layout, const CblasTranspose &TransA,
        size_t m, size_t n, cxdbl alpha, const cxdbl *A, size_t lda,
        const cxdbl *x, int incx, cxdbl beta, cxdbl *y, int incy);

    void lvl2_dger(
        const CblasOrder &layout,
        size_t m, size_t n, double alpha, const double *x, size_t inc_x,
        const double *y, size_t inc_y, double *a, size_t lda);

    void lvl2_zgeru(
        const CblasOrder &layout,
        size_t m, size_t n, cxdbl alpha, const cxdbl *x, size_t inc_x,
        const cxdbl *y, size_t inc_y, cxdbl *a, size_t lda);

    void lvl2_zgerc(
        const CblasOrder &layout,
        size_t m, size_t n, cxdbl alpha, const cxdbl *x, size_t inc_x,
        const cxdbl *y, size_t inc_y, cxdbl *a, size_t lda);
  } // namespace v1
} // namespace matrix

#endif

