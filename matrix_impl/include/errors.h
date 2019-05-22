#ifndef _MATRIX_ERRORS_H
#define _MATRIX_ERRORS_H

#include <exception>


namespace matrix {
  inline namespace v1 {
    class MatrixColsSizeUnevenError : public std::length_error {
    public:
      using std::length_error::length_error;
    };
  }
}

#endif
