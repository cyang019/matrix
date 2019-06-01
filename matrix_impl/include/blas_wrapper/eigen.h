#ifndef _MATRIX_EIGEN_H
#define _MATRIX_EIGEN_H

#include "blas_wrapper/cblas_common.h"
#include "common.h"


namespace matrix {
  inline namespace v1 {
    enum class EigenMethod : int {
      // hermitian
      zheev         = 1,
      zheev_2stage  = 2,
      zheevd        = 3,
      zheevd_2stage = 4,
      zheevr        = 5,
      zheevr_2stage = 6,
      zheevx        = 7,
      zheevx_2stage = 8
    };

    void mat_zheev();
    void mat_zheev_2stage();
    void mat_zheevd();
    void mat_zheevd_2stage();
    void mat_zheevr();
    void mat_zheevr_2stage();
    void mat_zheevx();
    void mat_zheevx_2stage();

  } // namespace v1
} // namespace matrix

#endif
