#ifndef _MATRIX_MATRIX_H
#define _MATRIX_MATRIX_H

#include <memory>   // unique_ptr
#include <utility>
#include <exception>
#include <type_traits>
#include <initializer_list>
#include <complex>
#include <random>
#include <cassert>
#include <cmath>    // abs
#include "common.h"
#include "errors.h"
#include "blas_wrapper/cblas_functions.h"

// temperary
#include <iostream>


namespace matrix {
    inline namespace v1 {
        template<typename T>
        class Matrix;

        //basic operations
        template<typename T>
        bool operator==(const Matrix<T> &, const Matrix<T> &);

        template<typename T>
        bool operator!=(const Matrix<T> &, const Matrix<T> &);

        bool allclose(const Matrix<double> &, const Matrix<double> &, double);

        template<typename T>
        Matrix<T> operator+(const Matrix<T> &, const Matrix<T> &);

        template<typename T>
        Matrix<T> operator+(const Matrix<T> &, const T &);

        template<typename T>
        Matrix<T> operator+(const T&, const Matrix<T> &);

        template<typename T>
        Matrix<T> operator-(const Matrix<T> &, const Matrix<T> &);

        template<typename T>
        Matrix<T> operator-(const Matrix<T> &, const T &);

        template<typename T>
        Matrix<T> operator-(const T&, const Matrix<T> &);

        template<typename T>
        Matrix<T> operator*(const Matrix<T> &, const Matrix<T> &);

        template<typename T>
        Matrix<T> operator*(const Matrix<T> &, const T &);

        template<typename T>
        Matrix<T> operator*(const T&, const Matrix<T> &);

        template<typename T>
        Matrix<T> operator/(const Matrix<T> &, const T &);

        Matrix<cxdbl> exp(const Matrix<cxdbl> &);

        // @brief: column major matrix type.
        template<typename T>
        class Matrix{
        public:
          friend bool operator==<T>(const Matrix<T> &, const Matrix<T> &);
          friend bool operator!=<T>(const Matrix<T> &, const Matrix<T> &);
          friend Matrix<T> operator+<T>(const Matrix<T> &t_m1, const Matrix<T> &t_m2);
          friend Matrix<T> operator+<T>(const T &t_v1, const Matrix<T> &t_m2);
          friend Matrix<T> operator+<T>(const Matrix<T> &t_m1, const T &t_v2);
          friend Matrix<T> operator-<T>(const Matrix<T> &t_m1, const Matrix<T> &t_m2);
          friend Matrix<T> operator-<T>(const T &t_v1, const Matrix<T> &t_m2);
          friend Matrix<T> operator-<T>(const Matrix<T> &t_m1, const T &t_v2);
          friend Matrix<T> operator*<T>(const Matrix<T> &t_m1, const Matrix<T> &t_m2);
          friend Matrix<T> operator*<T>(const T &t_v1, const Matrix<T> &t_m2);
          friend Matrix<T> operator*<T>(const Matrix<T> &t_m1, const T &t_v2);
          friend Matrix<T> operator/<T>(const Matrix<T> &t_m1, const T &t_v2);

          Matrix();
          Matrix(size_t, size_t);
          Matrix(size_t);     ///< square matrix
          Matrix(std::initializer_list<std::initializer_list<T>> il);
          Matrix(size_t n_total, const T *raw_data);

          // copy ctor
          Matrix(const Matrix<T> &);

          // move ctor
          Matrix(Matrix<T> &&) noexcept;

          // copy ctor
          Matrix<T>& operator=(const Matrix<T> &);

          // move ctor
          Matrix<T>& operator=(Matrix<T> &&) noexcept;
          ~Matrix();
          
          Matrix<T>& swapCols(size_t, size_t);
          Matrix<T>& swapRows(size_t, size_t);
          Matrix<T>& operator+=(const Matrix<T> &);
          Matrix<T>& operator+=(const T &);
          Matrix<T>& operator-=(const Matrix<T> &);
          Matrix<T>& operator-=(const T &);

          Matrix<T>& operator*=(const T &);
          Matrix<T>& operator/=(const T &);
          Matrix<T>& setZero();
          Matrix<T>& setOne();
          Matrix<T>& setRandom();

          T& operator()(size_t, size_t) const;
          T& operator()(size_t, size_t);

          // Transpose
          Matrix<T> t() const;
          Matrix<T>& tInplace();

          // Adjoint
          Matrix<T> adjoint() const;
          Matrix<T>& adjointInplace();

          T* data();
          const T* data() const;

          size_t nrows() const;
          size_t ncols() const;

          std::tuple<size_t, size_t> shape() const;
        private:
          size_t m_nrows;
          size_t m_ncols;
          std::unique_ptr<T[]> m_data; ///< serialized matrix
        };
    }
}

#include "MatrixImpl.hpp"

#include "OperatorsImpl.hpp"


#endif
