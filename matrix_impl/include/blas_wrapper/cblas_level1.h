#ifndef _MATRIX_BLAS_LEVEL1_H
#define _MATRIX_BLAS_LEVEL1_H

#include "blas_wrapper/cblas_common.h"
#include "common.h"

namespace matrix {
  inline namespace v1 {
    void lvl1_dscal(size_t n, double alpha, double *x, int incx);
    void lvl1_zscal(size_t n, cxdbl alpha, cxdbl *x, int incx);
    //void lvl1_zscal(size_t n, ComplexDbl alpha, ComplexDbl *x, int incx);
    void lvl1_zdscal(size_t n, double alpha, cxdbl *x, int incx);
    //void lvl1_zdscal(size_t n, double alpha, ComplexDbl *x, int incx);

    void lvl1_dcopy(size_t n, const double *x, int inc_x, double *y, int inc_y);
    void lvl1_zcopy(size_t n, const cxdbl *x, int inc_x, cxdbl *y, int inc_y);
    //void lvl1_zcopy(size_t n, const ComplexDbl *x, int inc_x, ComplexDbl *y, int inc_y);

    double lvl1_ddot(size_t n, const double *x, int inc_x, const double *y, int inc_y);
    cxdbl lvl1_zdot(size_t n, const cxdbl *x, int inc_x, const cxdbl *y, int inc_y);
    //ComplexDbl lvl1_zdot(size_t n, const ComplexDbl *x, int inc_x, const ComplexDbl *y, int inc_y);

    void lvl1_daxpy(size_t n, double alpha, const double *x, int inc_x, double *y, int inc_y);
    void lvl1_zaxpy(size_t n, cxdbl alpha, const cxdbl *x, int inc_x, cxdbl *y, int inc_y);
    //void lvl1_zaxpy(size_t n, ComplexDbl alpha, const ComplexDbl *x, int inc_x, ComplexDbl *y, int inc_y);
  }
}


#endif
