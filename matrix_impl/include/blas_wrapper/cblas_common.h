#ifndef _MATRIX_CBLAS_COMMON_H
#define _MATRIX_CBLAS_COMMON_H


namespace matrix {
  inline namespace v1 {
    enum class CblasOrder : int {
      CblasRowMajor=101,
      CblasColMajor=102
    };
    
    enum class CblasTranspose : int {
      CblasNoTrans = 111,
      CblasTrans=112,
      CblasConjTrans=113
    };

  }
}

#endif
