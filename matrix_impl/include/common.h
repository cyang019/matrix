#ifndef _MATRIX_COMMON_H
#define _MATRIX_COMMON_H

#include <type_traits>
#include <complex>
#include <limits>

namespace matrix {
  inline namespace v1 {
    using cxdbl = std::complex<double>;

    constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr int int_max = std::numeric_limits<int>::max();
    constexpr cxdbl cx_zero = 0;

    // complex double
    // To be compatible with lapack
    //using ComplexDbl = __CLPK_doublecomplex; ///< struct { double r, i; };

    //ComplexDbl operator+(const ComplexDbl &lhs, const ComplexDbl &rhs);
    //ComplexDbl operator-(const ComplexDbl &lhs, const ComplexDbl &rhs);
    //ComplexDbl operator*(const ComplexDbl &lhs, const ComplexDbl &rhs);
    //ComplexDbl operator/(const ComplexDbl &lhs, const ComplexDbl &rhs);
    //ComplexDbl conj(const ComplexDbl &num);
    //ComplexDbl getComplexDblZero();
    //ComplexDbl getComplexDblOne();
    //ComplexDbl getComplexDblNegOne();

    template <typename T> using static_not = std::integral_constant<bool, !T::value>;

    template <typename T> struct is_default : std::true_type {};

    template<>
    struct is_default<double> : std::false_type {};

    template<>
    struct is_default<cxdbl> : std::false_type {};

    template<typename T> struct is_complex : std::false_type {};
    template<typename T> struct is_complex<std::complex<T>> 
      : std::true_type {};

    template<typename T> struct is_double : std::false_type {};
    template<>
    struct is_double<double> : std::true_type {};

    template<typename T> struct is_int : std::false_type {};
    template<>
    struct is_int<int> : std::true_type {};

    constexpr long long int factorial(int n){
      long long int res = 1;
      while(n > 1){
        res *= n;
        --n;
      }
      return res;
    }
  }
}

#endif
