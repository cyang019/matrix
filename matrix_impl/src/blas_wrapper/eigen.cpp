#include "clapack.h"
#include "blas_wrapper/eigen.h"
#include "errors.h"


namespace matrix {
  inline namespace v1 {
    int mat_zheev(char jobz, char uplo, size_t n, ComplexDbl *a,
        size_t lda, 
        double *w, ComplexDbl *work, size_t lwork, 
        double *rwork, int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }

      if(lwork > (size_t)int_max){
        throw IndexOutOfBound("lwork cannot exceed INT_MAX.");
      }
#endif
      int dim = (int)n;
      int lower_dim = (int)lda;
      int work_length = (int)lwork;
      auto clpk_a = static_cast<__CLPK_doublecomplex *>(a);
      auto clpk_work = static_cast<__CLPK_doublecomplex *>(work);
      int res = zheev_(&jobz, &uplo, &dim, clpk_a, &lower_dim, 
          w, clpk_work, &work_length, rwork, info);
      return res;
    }

  } // namespace v1
} // namespace matrix

