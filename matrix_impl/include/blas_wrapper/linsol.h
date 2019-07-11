#ifndef _MATRIX_LINSOL_H
#define _MATRIX_LINSOL_H

#include "blas_wrapper/cblas_common.h"
#include "matrix_core/common.h"


namespace matrix {
  inline namespace v1 {
    int mat_zgesv(
        size_t n, size_t nrhs, cxdbl *a, size_t lda,
        int *ipiv, cxdbl *b, size_t ldb, int *info);

    int mat_dgesv(
        size_t n, size_t nrhs, double *a, size_t lda,
        int *ipiv, double *b, size_t ldb, int *info);
  }   // namespace v1
}   // namespace matrix

#endif
