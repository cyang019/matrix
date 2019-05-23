#include "blas_wrapper/blas_level1.h"


namespace matrix {
  inline namespace v1 {

    void matrix_dcopy(size_t n, const double *x, int inc_x, double *y, int inc_y)
    {
      while(n > (size_t)int_max){
        cblas_dcopy(int_max, x, inc_x, y, inc_y);

        n -= int_max;
        x += int_max;
        y += int_max;
      }
      cblas_dcopy(n, x, inc_x, y, inc_y);
    }

    void matrix_zcopy(size_t n, const cxdbl *x, int inc_x, cxdbl *y, int inc_y)
    {
      while(n > (size_t)int_max){
        cblas_zcopy(int_max, x, inc_x, y, inc_y);

        n -= int_max;
        x += int_max;
        y += int_max;
      }
      cblas_zcopy(n, x, inc_x, y, inc_y);
    }
  }
}

