#include "blas_wrapper/blas_level2.h"
#include "errors.h"


namespace matrix {
  inline namespace v1 {
    void lvl2_cblas_dgemv(
        const CBLAS_ORDER &layout, const CBLAS_TRANSPOSE &TransA,
        size_t m, size_t n, double alpha, const double *A, size_t lda,
        const double *x, int incx, double beta, double *y, int incy)
    {
#ifndef NDEBUG
      if(m >= int_max || n >= int_max || lda >= int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      cblas_dgemv(layout, TransA, (int)m, (int)n, alpha, A, (int)lda, x, incx, beta, y, incy);
    }

    void lvl2_cblas_zgemv(
        const CBLAS_ORDER &layout, const CBLAS_TRANSPOSE &TransA,
        size_t m, size_t n, cxdbl alpha, const cxdbl *A, size_t lda,
        const cxdbl *x, int incx, cxdbl beta, cxdbl *y, int incy)
    {
#ifndef NDEBUG
      if(m >= int_max || n >= int_max || lda >= int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      cblas_zgemv(layout, TransA, (int)m, (int)n, &alpha, A, (int)lda, x, incx, &beta, y, incy);
    }

  } // namespace v1
} // namespace matrix
