#ifndef _MATRIX_EIGEN_H
#define _MATRIX_EIGEN_H

#include "blas_wrapper/cblas_common.h"
#include "matrix_core/common.h"


namespace matrix {
  inline namespace v1 {
    enum class EigenMethod : int {
      // hermitian
      zheev         = 1,
      //zheev_2stage  = 2,
      zheevd        = 3,
      //zheevd_2stage = 4,
      zheevr        = 5,
      //zheevr_2stage = 6,
      zheevx        = 7,
      //zheevx_2stage = 8
    };

    int mat_zheev(char JOBZ, char UPLO, size_t n, cxdbl *a,
        size_t lda, 
        double *w, cxdbl *work, size_t lwork, 
        double *rwork, int *info);

    int mat_zheevd(char jobz, char uplo, size_t n, cxdbl *a,
        size_t lda,
        double *w, cxdbl *work, size_t lwork,
        double *rwork, size_t lrwork, int*iwork, size_t liwork,
        int *info);

    int mat_zheevr(char jobz, char range, char uplo, size_t n, cxdbl *a,
        size_t lda, double vl, double vu, size_t il, size_t iu, double abstol,
        size_t m, double *w, cxdbl *z, size_t ldz, int *isuppz,
        cxdbl *work, size_t lwork, double *rwork, size_t lrwork, int *iwork, size_t liwork,
        int *info);
    // not implemented
    //int mat_zheev_2stage(char JOBZ, char UPLO, size_t n, cxdbl *a,
    //    size_t lda, 
    //    double *w, cxdbl *work, size_t lwork, 
    //    double *rwork, int *info);
    //int mat_zheevd_2stage();
    //void mat_zheevr_2stage();
    void mat_zheevx();
    //void mat_zheevx_2stage();

    // helpers
    // make hermitian matrix upper or lower triangular
    int zhetrd(char uplo, size_t n, cxdbl *a, size_t lda,
        double *d, double *e, cxdbl *tau,
        cxdbl *work, size_t lwork, int *info);

  } // namespace v1
} // namespace matrix

#endif
