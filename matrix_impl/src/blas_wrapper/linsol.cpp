#include "clapack.h"
#include "blas_wrapper/linsol.h"
#include "errors.h"


namespace matrix {
  inline namespace v1 {
    int mat_zgesv(
        size_t n, size_t nrhs, cxdbl *a, size_t lda,
        int *ipiv, cxdbl *b, size_t ldb, int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max || ldb > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      int dim1 = (int)n;
      int dim2 = (int)nrhs;
      int lower_dim_a = (lda > 0) ? (int)lda : 1;
      int lower_dim_b = (ldb > 0) ? (int)ldb : 1;

      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_b = reinterpret_cast<__CLPK_doublecomplex *>(b);
      int res = zgesv_(&dim1, &dim2, clpk_a, &lower_dim_a,
          ipiv, clpk_b, &lower_dim_b, info);
      return res;
    }

    int mat_dgesv(
        size_t n, size_t nrhs, double *a, size_t lda,
        int *ipiv, double *b, size_t ldb, int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max || ldb > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      int dim1 = (int)n;
      int dim2 = (int)nrhs;
      int lower_dim_a = (lda > 0) ? (int)lda : 1;
      int lower_dim_b = (ldb > 0) ? (int)ldb : 1;

      int res = dgesv_(&dim1, &dim2, a, &lower_dim_a,
          ipiv, b, &lower_dim_b, info);
      return res;
    }
  } // namespace v1
} // namespace matrix
