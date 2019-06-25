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
#include <cmath>    // abs, max, min
#include <vector>
#include <unordered_set>
#include <algorithm>  // sort
#include <functional> // less
#include <limits>
#include <iostream>
#include <sstream>
#include <cstdint>   // uint32_t
#include <cstdlib>    // rand, srand
#include "common.h"
#include "errors.h"
#include "blas_wrapper/cblas_functions.h"


namespace matrix {
    inline namespace v1 {
        template<typename T>
        class Matrix;

        //basic operations
        template<typename T>
        std::ostream& operator<<(std::ostream &os, const Matrix<T> &m);

        template<typename T>
        bool operator==(const Matrix<T> &, const Matrix<T> &);

        template<typename T>
        bool operator!=(const Matrix<T> &, const Matrix<T> &);

        bool allclose(const Matrix<double> &, const Matrix<double> &, double);
        bool allclose(const Matrix<cxdbl> &, const Matrix<cxdbl> &, double);

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

        // -------------------------
        /// advanced operations
        // -------------------------
        template<typename T>
        Matrix<T> abs(const Matrix<T> &);

        template<typename T>
        Matrix<T> exp(const Matrix<T> &);

        /// solves for A * X = B
        /// A is a square matrix
        /// A & B will be destroyed after this function
        template<typename T>
        Matrix<T> linearSolveSq(Matrix<T> &, Matrix<T> &);

        template<EigenMethod em>
        std::tuple<Matrix<double>, Matrix<cxdbl>> eigenSys(const Matrix<cxdbl> &);

        template<EigenMethod em>
        Matrix<double> eigenVal(const Matrix<cxdbl> &);

        template<typename T>
        Matrix<T> kroneckerProduct(const Matrix<T> &, const Matrix<T> &);

        /// trace of A^T. B or A^H . B
        /// based on xdot & xdotc
        template<typename T>
        T projection(const Matrix<T> &, const Matrix<T> &);

        template<typename T>
        T projectionNorm(const Matrix<T> &, const Matrix<T> &);

        template<typename T>
        T trace(const Matrix<T> &);

        template<typename T>
        Matrix<T> pow(const Matrix<T> &, std::uint64_t);

        // the maximum absolute column sum of a matrix
        template<typename T>
        T norm1(const Matrix<T> &);

        // @brief: column major matrix type.
        template<typename T>
        class Matrix{
        public:
          friend std::ostream& operator<<<T>(std::ostream &os, const Matrix<T> &m);
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
          Matrix<T>& operator*=(const Matrix<T> &);
          Matrix<T>& operator/=(const T &);
          Matrix<T>& setZero();
          Matrix<T>& setOne();
          Matrix<T>& setRandom();

          Matrix<T>& row(size_t);
          Matrix<T>& col(size_t);

          T& operator()(size_t, size_t) const;
          T& operator()(size_t, size_t);

          // Transpose
          Matrix<T> t() const;
          Matrix<T>& tInplace();

          // Adjoint
          Matrix<T> adjoint() const;
          Matrix<T>& adjointInplace();

          // Inverse
          Matrix<T> inverse() const;
          Matrix<T>& inverseInplace();

          T* data();
          const T* data() const;

          size_t nrows() const;
          size_t ncols() const;
          size_t nelements() const { return m_nrows * m_ncols; }

          T trace() const;

          std::tuple<size_t, size_t> shape() const;
        private:
          size_t m_nrows;
          size_t m_ncols;
          std::unique_ptr<T[]> m_data; ///< serialized matrix
        };

        template<typename T>
        Matrix<T> diagonal(std::initializer_list<T> vals);

        template<typename T>
        Matrix<T> identity(size_t n);

        template<typename T>
        Matrix<T> zeros(size_t nrows, size_t ncols);

        template<typename T>
        Matrix<T> ones(size_t nrows, size_t ncols);

        /// for double mat, either 1 or -1
        /// for complex<double> mat, a_ij/|a_ij|, sign(0) = 1
        template<typename T>
        Matrix<T> sign(const Matrix<T> &);
    }
}

#include "MatrixImpl.hpp"

#include "OperatorsImpl.hpp"

#include "AdvancedFunctionsImpl.hpp"

#include "MatrixExponential.hpp"


#endif
