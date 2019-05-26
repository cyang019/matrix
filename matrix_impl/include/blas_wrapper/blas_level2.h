#ifndef _MATRIX_BLAS_LEVEL2_H
#define _MATRIX_BLAS_LEVEL2_H

#include "cblas.h"
#include "common.h"


namespace matrix {
  inline namespace v1 {
    /// m & n cannot exceed int_max
    void lvl2_cblas_dgemv(
        const CBLAS_ORDER &layout, const CBLAS_TRANSPOSE &TransA,
        size_t m, size_t n, double alpha, const double *A, size_t lda,
        const double *x, int incx, double beta, double *y, int incy);

    void lvl2_cblas_zgemv(
        const CBLAS_ORDER &layout, const CBLAS_TRANSPOSE &TransA,
        size_t m, size_t n, cxdbl alpha, const cxdbl *A, size_t lda,
        const cxdbl *x, int incx, cxdbl beta, cxdbl *y, int incy);
                
  } // namespace v1
} // namespace matrix

#endif

