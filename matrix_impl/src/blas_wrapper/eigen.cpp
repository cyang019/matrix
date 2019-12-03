#include "configure_matrix.h"
#ifdef HAVE_CLAPACK
  #include "clapack.h"
#elif defined HAVE_LAPACKE
  #include "lapacke.h"
#else
#endif
#include "blas_wrapper/eigen.h"
#include "matrix_core/errors.h"


namespace matrix {
  inline namespace v1 {
    int mat_zheev(char jobz, char uplo, size_t n, cxdbl *a,
        size_t lda, 
        double *w, cxdbl *work, size_t lwork, 
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
      /// unique_ptr<cxdbl> a; a.get() directly
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work);
      int res = zheev_(&jobz, &uplo, &dim, clpk_a, &lower_dim, 
          w, clpk_work, &work_length, rwork, info);
      return res;
    }

    int mat_zheevd(char jobz, char uplo, size_t n, cxdbl *a,
        size_t lda,
        double *w, cxdbl *work, size_t lwork,
        double *rwork, size_t lrwork, int *iwork, size_t liwork,
        int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }

      if(lwork > (size_t)int_max || lrwork > (size_t)int_max || liwork > (size_t)int_max){
        throw IndexOutOfBound("lwork cannot exceed INT_MAX.");
      }
#endif
      int dim = (int)n;
      int lower_dim = (int)lda;
      int work_length = (int)lwork;
      int rwork_length = (int)lrwork;
      int iwork_length = (int)liwork;

      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work);

      int res = zheevd_(&jobz, &uplo, &dim, clpk_a, &lower_dim, w,
          clpk_work, &work_length, rwork, &rwork_length,
          iwork, &iwork_length, info);
      return res;
    }

    int mat_zheevr(char jobz, char range, char uplo, size_t n, cxdbl *a,
        size_t lda, double vl, double vu, size_t il, size_t iu, double abstol,
        size_t m, double *w, cxdbl *z, size_t ldz, int *isuppz,
        cxdbl *work, size_t lwork, double *rwork, size_t lrwork, int *iwork, size_t liwork,
        int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }

      if(lwork > (size_t)int_max || lrwork > (size_t)int_max || liwork > (size_t)int_max){
        throw IndexOutOfBound("lwork cannot exceed INT_MAX.");
      }

      if(m > (size_t)int_max) 
        throw IndexOutOfBound("number of eigen values cannot exceed matrix dimension");

      if(il > (size_t)int_max || iu > (size_t)int_max){
        throw IndexOutOfBound("bound(s) of eigen value count cannot exceed matrix dimension");
      }
#endif
      int dim = (int)n;
      int lower_dim = (int)lda;
      int work_length = (int)lwork;
      int rwork_length = (int)lrwork;
      int iwork_length = (int)liwork;
      int n_eigen = (int)m;
      int lead_dim_z = (int)ldz;
      int lower_idx = (int)il;
      int upper_idx = (int)iu;

      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work);
      auto clpk_z = reinterpret_cast<__CLPK_doublecomplex *>(z);

      int res = zheevr_(&jobz, &range, &uplo, &dim, clpk_a, &lower_dim,
          &vl, &vu, &lower_idx, &upper_idx, &abstol, &n_eigen, w, clpk_z,
          &lead_dim_z, isuppz, clpk_work, &work_length, rwork, &rwork_length,
          iwork, &iwork_length, info);
      return res;
    }

    // out: a, d, e, tau, work, info
    int zhetrd(char uplo, size_t n, cxdbl *a, size_t lda,
        double *d, double *e, cxdbl *tau,
        cxdbl *work, size_t lwork, int *info)
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
      /// unique_ptr<cxdbl> a; a.get() directly
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_tau = reinterpret_cast<__CLPK_doublecomplex *>(tau);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work);
      int res = zhetrd_(&uplo, &dim, clpk_a, &lower_dim, d, e,
          clpk_tau, clpk_work, &work_length, info);
      return res;
    }
  } // namespace v1
} // namespace matrix

