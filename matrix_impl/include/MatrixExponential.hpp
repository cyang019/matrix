#include "matrix_exp_impl/indices.hpp"
#include "matrix_exp_impl/norm_estimate.hpp"

namespace matrix { inline namespace v1 {
  namespace exponential {
    // Al-Mohy, Awad H., and Nicholas J. Higham. "A new scaling and squaring 
    // algorithm for the matrix exponential." SIAM Journal on Matrix Analysis 
    // and Applications 31.3 (2009): 970-989.
    template<typename T>
    Matrix<T> exp(const Matrix<T> &mat)
    {
      const auto A0 = identity(mat.nrows());
      const auto A2 = mat * mat;
      double d6 = std::pow(normest(A2, 3), 1.0/6);
      const double eta_1 = std::max(std::pow(normest(A2, 2), 0.25), d6);
      if(eta_1 < theta_3 + eps && ell(mat, 3) == 0){
        Matrix<T> u3 = mat * (b(3, 1) * A0 + b(3, 3) * A2);
        Matrix<T> v3 = b(3, 0) * A0 + b(3, 2) * A2;
        Matrix<T> p3 = u3 + v3;
        Matrix<T> q3 = v3 - u3;
        Matrix<T> r3 = linearSolveSq(q3, p3);
        return r3;
      }
      const auto A4 = A2 * A2;
      const double d4 = std::pow(norm1(A4), 0.25);
      const double eta_2 = std::max(d4, d6);
      if(eta_2 < theta_5 + eps && ell(mat, 5) == 0){
        Matrix<T> u5 = mat * (b(5, 1) * A0 + b(5, 3) * A2 + b(5,5) * A4);
        Matrix<T> v5 = b(5, 0) * A0 + b(5, 2) * A2 + b(5, 4) * A4;
        Matrix<T> p5 = u5 + v5;
        Matrix<T> q5 = v5 - u5;
        Matrix<T> r5 = linearSolveSq(q5, p5);
        return r5;
      }
      const auto A6 = A4 * A2;
      d6 = std::pow(norm1(A6), 1.0/6);
      const double d8 = std::pow(normest(A4, 2), 0.125);
      const double eta_3 = std::max(d6, d8);
      if(eta_3 < theta_7 + eps && ell(mat, 7) == 0){
        Matrix<T> u7 = mat * (b(7, 1) * A0 + b(7, 3) * A2 + b(7, 5) * A4 + b(7, 7) * A6);
        Matrix<T> v7 = b(7, 0) * A0 + b(7, 2) * A2 + b(7, 4) * A4 + b(7, 6) * A6;
        Matrix<T> p7 = u7 + v7;
        Matrix<T> q7 = v7 - u7;
        Matrix<T> r7 = linearSolveSq(q7, p7);
        return r7;
      }
      if(eta_3 < theta_9 + eps && ell(mat, 9) == 0){
        const auto A8 = A4 * A4;
        Matrix<T> u9 = mat * (b(9, 1) * A0 + b(9, 3) * A2 + b(9, 5) * A4 + b(9, 7) * A6 + b(9, 9) * A8);
        Matrix<T> v9 = b(9, 0) * A0 + b(9, 2) * A2 + b(9, 4) * A4 + b(9, 6) * A6 + b(9, 8) * A8;
        Matrix<T> p9 = u9 + v9;
        Matrix<T> q9 = v9 - u9;
        Matrix<T> r9 = linearSolveSq(q9, p9);
        return r9;
      }

      const double eta_4 = std::max(d8, std::pow(norm1(A4 * A6), 0.1));
      const double eta_5 = std::min(eta_3, eta_4);

      const double val = std::log2(eta_5/theta_13);
      int s = (int)std::max(std::ceil(val), 0.0);
      s += ell(std::pow(2, -s) * mat, 13);
      mat *= std::pow(0.5, s);
      A2 *= std::pow(0.5, 2*s); 
      A4 *= std::pow(0.5, 4*s);
      A6 *= std::pow(0.5, 6*s);

      Matrix<T> u13 = mat * (A6 * (b(13, 13) * A6 + b(13, 11) * A4 + b(13, 9) * A2) 
                           + b(13, 7) * A6 + b(13, 5) * A4 + b(13, 3) * A2
                           + b(13, 1) * A0);
      Matrix<T> v13 = A6 * (b(13, 12) * A6 + b(13, 10) * A4 + b(13, 8) * A2)
                    + b(13, 6) * A6 + b(13, 4) * A4 + b(13, 2) * A2 
                    + b(13, 0) * A0;
      Matrix<T> p13 = u13 + v13;
      Matrix<T> q13 = v13 - u13;
      Matrix<T> r13 = linearSolveSq(q13, p13);
      Matrix<T> res = pow(r13, std::pow(2, s));
      return res;
    }
  }
  // squaring and scaling & Pade Approximation
  template<typename T>
  Matrix<T> exp(const Matrix<T> &mat)
  {
    return exponential::exp(mat);
  }

}  // namespace v1
} // namespace matrix

