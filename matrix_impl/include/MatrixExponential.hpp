#include "matrix_exp_impl/norm_estimate.hpp"
#include "matrix_exp_impl/indices.hpp"

namespace matrix { inline namespace v1 {
  namespace exponential {
    // Al-Mohy, Awad H., and Nicholas J. Higham. "A new scaling and squaring 
    // algorithm for the matrix exponential." SIAM Journal on Matrix Analysis 
    // and Applications 31.3 (2009): 970-989.
    Matrix<T> exp(const Matrix<T> &mat)
    {
      const auto A0 = identity(mat.nrows());
      const auto A2 = mat * mat;
      const double d_6 = std::pow(normest(A2, 3), 1.0/6);
      const double eta_1 = std::max(std::pow(normest(A2, 2), 0.25), d_6);
      if(eta_1 < theta_3 + eps && ell(mat, 3) == 0){
        Matrix<T> u3 = mat * (b(3, 1) * A0 + b(3, 3) * A2);
        Matrix<T> v3 = b(3, 0) * A0 + b(3, 2) * A2;
        Matrix<T> p3 = u3 + v3;
        Matrix<T> q3 = v3 - u3;
        Matrix<T> r3 = linearSolveSq(q3, p3);
        return r3;
      }
      const auto A4 = A2 * A2;
      const double d_4 = std::pow(norm1(A4), 0.25);
      const double eta_2 = std::max(d4, d6);
      if(eta2 < theta_3 + eps && ell(mat, 5) == 0){
        Matrix<T> u5 = mat * (b(5, 1) * A0 + b(5, 3) * A2 + b(5,5) * A4);
        Matrix<T> v5 = b(5, 0) * A0 + b(5, 2) * A2 + b(5, 4) * A4;
        Matrix<T> p5 = u5 + v5;
        Matrix<T> q5 = v5 - u5;
        Matrix<T> r5 = linearSolveSq(q5, p5);
        return r5;
      }
      const auto A6 = A4 * A2;
      const double d6 = std::pow(norm1(A6), 1.0/6);
      const double d8 = std::pow(normest(A4, 2), 0.125);
      const double eta_3 = std::max(d6, d8);
    }
  }
  // squaring and scaling & Pade Approximation
  template<typename T>
  Matrix<T> exp(const Matrix<T> &mat)
  {
    const auto A2 = mat * mat;
    const auto A4 = A2 * A2;
    const auto A6 = A4 * A2;
    const double d_6 = std::pow(norm1(A6), 1.0/6);
    const double eta_1 =  std::max(std::pow(norm1(A6), 0.25), d_6);
    if (eta_1 <= theta_3 && true){
    }

    Matrix<T> res(mat.nrows(), mat.ncols());
    return res;
  }

}  // namespace v1
} // namespace matrix

