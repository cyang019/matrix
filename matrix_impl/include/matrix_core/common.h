#ifndef _MATRIX_COMMON_H
#define _MATRIX_COMMON_H

#include <type_traits>
#include <complex>
#include <limits>
#include <cstdint>
#include <ctime>

namespace matrix {
  inline namespace v1 {
    using cxdbl = std::complex<double>;

    constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr int int_max = std::numeric_limits<int>::max();
    constexpr cxdbl cx_zero = 0;
    inline std::int64_t log2int(std::int64_t value) 
    { std::int64_t target = 0; 
      while(value > 0) {
        ++target;
        value >>= 1;
      }
      return target;
    }

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
    bool approxEqual(const T &num1, const T &num2, double atol, double rtol=1.0e-14)
    {
      return (std::abs(num1 - num2) < atol, + rtol*std::abs(num2));
    }

    inline long double factorial(std::size_t n)
    {
      switch(n){
        case 0:
          return 1.0;
        case 1:
          return 1.0;
        case 2:
          return 2.0;
        case 3:
          return 6.0;
        case 4:
          return 24.0;
        case 5:
          return 120.0;
        case 6:
          return 720.0;
        case 7:
          return 5040.0;
        case 8:
          return 40320.0;
        case 9:
          return 362880.0;
        case 10:
          return 3628800.0;
        case 11:
          return 39916800.0;
        case 12:
          return 479001600.0;
        case 13:
          return 6227020800.0;
        case 14:
          return 87178291200.0;
        case 15:
          return 1.3076744e+12;
        default:
          {
            long double res = 1.3076744e+12;
            for(std::size_t i = 16; i <= n; ++i){
              res *= static_cast<long double>(i);
            }
            return res;
          }
          break;
      }
      return -1.0;
    }
  }   // namespace v1
}   // namespace matrix

#endif
