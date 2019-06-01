#include "common.h"

namespace matrix {
  inline namespace v1 {
    ComplexDbl operator+(const ComplexDbl &lhs, const ComplexDbl &rhs)
    {
      auto res = lhs;
      res.r += rhs.r;
      res.i += rhs.i;
      return res;
    }

    ComplexDbl operator-(const ComplexDbl &lhs, const ComplexDbl &rhs)
    {
      auto res = lhs;
      res.r -= rhs.r;
      res.i -= rhs.i;
      return res;
    }

    ComplexDbl operator*(const ComplexDbl &lhs, const ComplexDbl &rhs)
    {
      ComplexDbl res;
      res.r = lhs.r * rhs.r - lhs.i * rhs.i;
      res.i = lhs.r * rhs.i + lhs.i * rhs.r;
      return res;
    }

    ComplexDbl operator/(const ComplexDbl &lhs, const ComplexDbl &rhs)
    {
      ComplexDbl res;
      const double denominator = rhs.r * rhs.r + rhs.i * rhs.i;
      res.r = (rhs.r * lhs.r + rhs.i * lhs.i)/denominator;
      res.i = (rhs.i * lhs.r - rhs.r * lhs.i)/denominator;
      return res;
    }

    ComplexDbl conj(const ComplexDbl &num)
    {
      ComplexDbl res;
      res.r = num.r;
      res.i = -num.i;
      return res;
    }

    ComplexDbl getComplexDblZero()
    {
      ComplexDbl res;
      res.r = 0;
      res.i = 0;
      return res;
    }

    ComplexDbl getComplexDblOne()
    {
      ComplexDbl res;
      res.r = 1;
      res.i = 0;
      return res;
    }

    ComplexDbl getComplexDblNegOne()
    {
      ComplexDbl res;
      res.r = -1;
      res.i = 0;
      return res;
    }
  } // namespace v1
} // namespace matrix
