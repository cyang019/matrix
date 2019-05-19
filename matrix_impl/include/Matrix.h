#ifndef _MATRIX_MATRIX_H
#define _MATRIX_MATRIX_H

#include <memory>   // unique_ptr
#include <exception>
#include <type_traits>
#include <initializer_list>
#include <complex>
#include <random>
#include "blas_wrapper/blas_functions.h"


namespace matrix {
    inline namespace v1 {
        template<typename T> struct is_complex : std::false_type {};
        template<typename T> struct is_complex<std::complex<T>> 
          : std::true_type {};

        template<typename T> struct is_double : std::false_type {};
        template<>
        struct is_double<double> : std::true_type {};

        template<typename T> struct is_int : std::false_type {};
        template<>
        struct is_int<int> : std::true_type {};

        template<typename T>
        class Matrix;

        // basic operations
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

        Matrix<double> exp(const Matrix<double> &, const double &);

        /// @brief: column major matrix type.
        template<typename T>
        class Matrix{
        public:
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
            Matrix(const Matrix<T> &);
            Matrix(Matrix<T> &&) noexcept;
            Matrix<T>& operator=(const Matrix<T> &);
            Matrix<T>& operator=(Matrix<T> &&) noexcept;
            ~Matrix();
            
            Matrix<T>& swapCols(size_t, size_t);
            Matrix<T>& swapRows(size_t, size_t);
            Matrix<T>& operator+=(const Matrix<T> &);
            Matrix<T>& operator-=(const Matrix<T> &);
            Matrix<T>& operator*=(const Matrix<T> &);
            Matrix<T>& setZero();
            Matrix<T>& setOne();
            Matrix<T>& setRandom();

            T& operator()(size_t, size_t) const;
            T& operator()(size_t, size_t);

            // Transpose
            Matrix<T>& t();

            // Adjoint
            Matrix<T>& adjoint();

            size_t rows() const;
            size_t cols() const;
        private:
            std::unique_ptr<T[]> m_data; ///< serialized matrix
            size_t m_nrows;
            size_t m_ncols;
        };
    }
}

#include "MatrixImpl.hpp"

#endif
