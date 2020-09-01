#include "configure_matrix.h"
#ifdef HAVE_APPLE_LAPACK
#include "clapack.h"
#elif defined HAVE_CLAPACK
#include "clapack.h"
#elif defined HAVE_LAPACKE
#include "lapacke.h"
#else
#endif

#include "blas_wrapper/eigen.h"
#include "matrix_core/errors.h"
#include <memory>


namespace matrix {
  inline namespace v1 {
    int mat_zheev(char jobz, char uplo, size_t n, 
        cxdbl *a, size_t lda, 
        double *w, int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif
      int dim = (int)n;
      int lower_dim = (int)lda;

      int res = 0;
#if defined(HAVE_APPLE_LAPACK) || defined(HAVE_CLAPACK)
      const size_t lwork = n > 0 ? (n+1) * n : 1;
      auto work = std::make_unique<cxdbl[]>(lwork);
      const size_t rwork_size = n > 0 ? 3 * n - 2 : 1;
      auto rwork = std::make_unique<double[]>(rwork_size);
      int work_length = (int)lwork;
      /// unique_ptr<cxdbl> a; a.get() directly
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
      res = zheev_(&jobz, &uplo, &dim, clpk_a, &lower_dim, 
          w, clpk_work, &work_length, rwork.get(), info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_zheev(LAPACK_COL_MAJOR, jobz, uplo, dim,
          a, lower_dim, w);
      res = *info;
#endif
      return res;
    }

    int mat_zheevd(char jobz, char uplo, size_t n, cxdbl *a,
        size_t lda, double *w, int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif

      int res = 0;
#if defined(HAVE_APPLE_LAPACK) || defined(HAVE_CLAPACK)
      int dim = (int)n;
      int lower_dim = (int)lda;
      // query sizes
      int lwork = -1;
      int lrwork = -1;
      int liwork = -1;
      auto work = std::make_unique<cxdbl[]>(1);
      auto rwork = std::make_unique<double[]>(1);
      auto iwork = std::make_unique<int[]>(1);

      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
      res = zheevd_(&jobz, &uplo, &dim, clpk_a, &lower_dim, w,
          clpk_work, &lwork, rwork.get(), &lrwork,
          iwork.get(), &liwork, info);

      lwork = static_cast<int>(work[0].real());
      lrwork = static_cast<int>(rwork[0]);
      liwork = iwork[0];
      work = std::make_unique<cxdbl[]>(lwork);
      rwork = std::make_unique<double[]>(lrwork);
      iwork = std::make_unique<int[]>(liwork);
      clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
      // actual calculation
      res = zheevd_(&jobz, &uplo, &dim, clpk_a, &lower_dim, w,
          clpk_work, &lwork, rwork.get(), &lrwork,
          iwork.get(), &liwork, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_zheevd(LAPACK_COL_MAJOR, jobz, uplo, (int)n,
          a, (int)lda, w);
      res = *info;
#endif
      return res;
    }

    int mat_zheevr(char jobz, char range, char uplo, size_t n, cxdbl *a,
        size_t lda, double vl, double vu, size_t il, size_t iu, double abstol,
        size_t m, double *w, cxdbl *z, size_t ldz, int *isuppz,
        int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }

      if(m > (size_t)int_max) 
        throw IndexOutOfBound("number of eigen values cannot exceed matrix dimension");

      if(il > (size_t)int_max || iu > (size_t)int_max){
        throw IndexOutOfBound("bound(s) of eigen value count cannot exceed matrix dimension");
      }
#endif
      int dim = (int)n;
      int lower_dim = (int)lda;
      int n_eigen = (int)m;
      int lead_dim_z = (int)ldz;
      int lower_idx = (int)il;
      int upper_idx = (int)iu;

      int res = 0;
#if defined(HAVE_APPLE_LAPACK) || defined(HAVE_CLAPACK)
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_z = reinterpret_cast<__CLPK_doublecomplex *>(z);

      // query
      int lwork = -1;
      int lrwork = -1;
      int liwork = -1;
      auto work = std::make_unique<std::complex<double>[]>(1);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
      auto rwork = std::make_unique<double[]>(1);
      auto iwork = std::make_unique<int[]>(1);
      res = zheevr_(&jobz, &range, &uplo, &dim, clpk_a, &lower_dim,
          &vl, &vu, &lower_idx, &upper_idx, 
          &abstol, &n_eigen, w, clpk_z,
          &lead_dim_z, isuppz, 
          clpk_work, &lwork, rwork.get(), &lrwork,
          iwork.get(), &liwork, info);
      
      // calculate
      lwork = static_cast<int>(work[0].real());
      lrwork = static_cast<int>(rwork[0]);
      liwork = iwork[0];
      work = std::make_unique<cxdbl[]>(lwork);
      rwork = std::make_unique<double[]>(lrwork);
      iwork = std::make_unique<int[]>(liwork);
      clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());

      res = zheevr_(&jobz, &range, &uplo, &dim, clpk_a, &lower_dim,
          &vl, &vu, &lower_idx, &upper_idx, &abstol, &n_eigen, w, clpk_z,
          &lead_dim_z, isuppz, clpk_work, &lwork, rwork.get(), &lrwork,
          iwork.get(), &liwork, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_zheevr(
          LAPACK_COL_MAJOR, jobz, range, uplo, 
          dim, a, 
          lower_dim, vl, vu, lower_idx, 
          upper_idx, abstol, &n_eigen,
          w, z, lead_dim_z, isuppz);
      res = *info;
#endif
      return res;
    }

    // out: a, d, e, tau, work, info
    int zhetrd(char uplo, size_t n, cxdbl *a, size_t lda,
        double *d, double *e, cxdbl *tau, int *info)
    {
#ifndef NDEBUG
      if(n > (size_t)int_max || lda > (size_t)int_max){
        throw IndexOutOfBound("matrix dimension cannot exceed INT_MAX.");
      }
#endif

      int dim = (int)n;
      int lower_dim = (int)lda;
      int res = 0;
#if defined(HAVE_APPLE_LAPACK) || defined(HAVE_CLAPACK)
      int lwork = -1;
      /// unique_ptr<cxdbl> a; a.get() directly
      auto clpk_a = reinterpret_cast<__CLPK_doublecomplex *>(a);
      auto clpk_tau = reinterpret_cast<__CLPK_doublecomplex *>(tau);

      // query
      auto work = std::make_unique<cxdbl[]>(1);
      auto clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
      res = zhetrd_(&uplo, &dim, clpk_a, &lower_dim, d, e,
          clpk_tau, clpk_work, &lwork, info);

      // calculate
      lwork = static_cast<int>(work[0].real());
      work = std::make_unique<cxdbl[]>(lwork);
      clpk_work = reinterpret_cast<__CLPK_doublecomplex *>(work.get());
      res = zhetrd_(&uplo, &dim, clpk_a, &lower_dim, d, e,
          clpk_tau, clpk_work, &lwork, info);
#elif defined HAVE_LAPACKE
      *info = LAPACKE_zhetrd(
          LAPACK_COL_MAJOR, uplo, dim,
          a, lower_dim, d, e, tau);
      res = *info;
#endif
      return res;
    }
  } // namespace v1
} // namespace matrix

