#include "configure_matrix.h"
#include "cblas.h"
#include "blas_wrapper/cblas_level2.h"
#include "matrix_core/errors.h"
#include <stdexcept>
#include <cmath>


namespace matrix {
  inline namespace v1 {
    void lvl2_dgemv(
        const CblasOrder &layout, const CblasTranspose &TransA,
        size_t m, size_t n, double alpha, const double *A, size_t lda,
        const double *x, int incx, double beta, double *y, int incy)
    {
#ifndef NDEBUG
      if(m >= int_max || n >= int_max || lda >= int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      cblas_dgemv((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA,
          (int)m, (int)n, alpha, A, (int)lda, x, incx, beta, y, incy);
    }

    void lvl2_zgemv(
        const CblasOrder &layout, const CblasTranspose &TransA,
        size_t m, size_t n, cxdbl alpha, const cxdbl *A, size_t lda,
        const cxdbl *x, int incx, cxdbl beta, cxdbl *y, int incy)
    {
#ifndef NDEBUG
      if(m >= int_max || n >= int_max || lda >= int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      cblas_zgemv((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)TransA,
          (int)m, (int)n, &alpha, A, (int)lda, x, incx, &beta, y, incy);
    }
    void lvl2_dger(
        const CblasOrder &layout,
        size_t m, size_t n, double alpha, const double *x, size_t inc_x,
        const double *y, size_t inc_y, double *a, size_t lda)
    {
#ifndef NDEBUG
      if(m >= int_max || n >= int_max || lda >= int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      cblas_dger((CBLAS_ORDER)layout,
          (int)m, (int)n, alpha, x, (int)inc_x, y, (int)inc_y, 
          a, std::max(1, (int)lda));
    }

    void lvl2_zgeru(
        const CblasOrder &layout,
        size_t m, size_t n, cxdbl alpha, const cxdbl *x, size_t inc_x,
        const cxdbl *y, size_t inc_y, cxdbl *a, size_t lda)
    {
#ifndef NDEBUG
      if(m >= int_max || n >= int_max || lda >= int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      cblas_zgeru((CBLAS_ORDER)layout,
          (int)m, (int)n, &alpha, x, (int)inc_x, y, (int)inc_y,
          a, std::max(1, (int)lda));
    }

    void lvl2_zgerc(
        const CblasOrder &layout,
        size_t m, size_t n, cxdbl alpha, const cxdbl *x, size_t inc_x,
        const cxdbl *y, size_t inc_y, cxdbl *a, size_t lda)
    {
#ifndef NDEBUG
      if(m >= int_max || n >= int_max || lda >= int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      cblas_zgerc((CBLAS_ORDER)layout,
          (int)m, (int)n, &alpha, x, (int)inc_x, y, (int)inc_y, 
          a, std::max(1, (int)lda));
    }
  } // namespace v1
} // namespace matrix
