#ifndef _MATRIX_ERRORS_H
#define _MATRIX_ERRORS_H

#include <exception>


namespace matrix {
  inline namespace v1 {
    class MatrixColsSizeUnevenError : public std::length_error {
    public:
      using std::length_error::length_error;
    };

    class MatrixSizeMismatchError : public std::length_error {
    public:
      using std::length_error::length_error;
    };

    class IndexOutOfBound : public std::out_of_range {
    public:
      using std::out_of_range::out_of_range;
    };

    class MultiplicationError : public std::length_error {
    public:
      using std::length_error::length_error;
    };

    class NotImplementedError : public std::logic_error {
    public:
      using std::logic_error::logic_error;
    };

    class NonInvertibleMatrix : public std::logic_error {
    public:
      using std::logic_error::logic_error;
    };

    class InvalidEigenValue : public std::logic_error {
      public:
        using std::logic_error::logic_error;
    };
  }
}

#endif
