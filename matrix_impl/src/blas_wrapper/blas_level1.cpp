#include "blas_wrapper/blas_level1.h"


namespace matrix {
  inline namespace v1 {
    void lvl1_dscal(size_t n, double alpha, double *x, int incx)
    {
      while(n > (size_t)int_max){
        cblas_dscal(int_max, alpha, x, incx);

        n -= (size_t)int_max;
        x += int_max;
      }
      cblas_dscal(n, alpha, x, incx);
    }

    void lvl1_zscal(size_t n, cxdbl alpha, cxdbl *x, int incx)
    {
      while(n > (size_t)int_max){
        cblas_zscal(int_max, &alpha, x, incx);

        n -= (size_t)int_max;
        x += int_max;
      }
      cblas_zscal(n, &alpha, x, incx);
    }

    void lvl1_zdscal(size_t n, double alpha, cxdbl *x, int incx)
    {
      while(n > (size_t)int_max){
        cblas_zdscal(int_max, alpha, x, incx);

        n -= (size_t)int_max;
        x += int_max;
      }
      cblas_zdscal(n, alpha, x, incx);
    }

    void lvl1_dcopy(size_t n, const double *x, int inc_x, double *y, int inc_y)
    {
      while(n > (size_t)int_max){
        cblas_dcopy(int_max, x, inc_x, y, inc_y);

        n -= (size_t)int_max;
        x += int_max;
        y += int_max;
      }
      cblas_dcopy(n, x, inc_x, y, inc_y);
    }

    double lvl1_ddot(size_t n, const double *x, int inc_x, const double *y, int inc_y)
    {
      double res = 0;
      while(n > (size_t)int_max){
        res += cblas_ddot(int_max, x, inc_x, y, inc_y);

        n -= (size_t)int_max;
        x += int_max;
        y += int_max;
      }
      res += cblas_ddot(n, x, inc_x, y, inc_y);
      return res;
    }

    cxdbl lvl1_zdot(size_t n, const cxdbl *x, int inc_x, const cxdbl *y, int inc_y)
    {
      cxdbl res = 0;
      while(n > (size_t)int_max){
        cxdbl temp = 0;
        cblas_zdotc_sub(int_max, x, inc_x, y, inc_y, &temp);

        n -= (size_t)int_max;
        x += int_max;
        y += int_max;
        res += temp;
      }
      cxdbl temp = 0;
      cblas_zdotc_sub(n, x, inc_x, y, inc_y, &temp);
      res += temp;
      return res;
    }

    void lvl1_zcopy(size_t n, const cxdbl *x, int inc_x, cxdbl *y, int inc_y)
    {
      while(n > (size_t)int_max){
        cblas_zcopy(int_max, x, inc_x, y, inc_y);

        n -= (size_t)int_max;
        x += int_max;
        y += int_max;
      }
      cblas_zcopy(n, x, inc_x, y, inc_y);
    }

    void lvl1_daxpy(size_t n, double alpha, const double *x, int inc_x, double *y, int inc_y)
    {
      while(n > (size_t)int_max){
        cblas_daxpy(int_max, alpha, x, inc_x, y, inc_y);

        n -= (size_t)int_max;
        x += int_max;
        y += int_max;
      }
      cblas_daxpy(n, alpha, x, inc_x, y, inc_y);
    }

    void lvl1_zaxpy(size_t n, cxdbl alpha, const cxdbl *x, int inc_x, cxdbl *y, int inc_y)
    {
      while(n > (size_t)int_max){
        cblas_zaxpy(int_max, &alpha, x, inc_x, y, inc_y);

        n -= (size_t)int_max;
        x += int_max;
        y += int_max;
      }
      cblas_zaxpy(n, &alpha, x, inc_x, y, inc_y);
    }
  }
}

