#ifndef _MATRIX_COMMON_H
#define _MATRIX_COMMON_H

#include <type_traits>
#include <complex>

namespace matrix {
  inline namespace v1 {
    using cxdbl = std::complex<double>;

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

  }
}

#endif
