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

    template<typename T>
    bool approxEqual(T num1, T num2, double e)
    {
      if constexpr(is_complex<T>::value){
        if(std::abs(num1.real() - num2.real()) > e) return false;
        if(std::abs(num1.imag() - num2.imag()) > e) return false;
        return true;
      } else {
        if(std::abs(num1 - num2) < e) return true;
        return false;
      }
    }
  }
}

#endif
