#ifndef _MATRIX_BLAS_LEVEL1_H
#define _MATRIX_BLAS_LEVEL1_H
#endif

#include "cblas.h"
#include "common.h"

namespace matrix {
  inline namespace v1 {
    void matrix_dscal(size_t n, double alpha, double *x, int incx);
    void matrix_zscal(size_t n, cxdbl alpha, cxdbl *x, int incx);

    void matrix_dcopy(size_t n, const double *x, int inc_x, double *y, int inc_y);
    void matrix_zcopy(size_t n, const cxdbl *x, int inc_x, cxdbl *y, int inc_y);

    void matrix_daxpy(size_t n, double alpha, const double *x, int inc_x, double *y, int inc_y);
    void matrix_zaxpy(size_t n, cxdbl alpha, const cxdbl *x, int inc_x, cxdbl *y, int inc_y);
  }
}


