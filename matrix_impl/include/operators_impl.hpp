namespace matrix {
  inline namespace v1 {
    template<typename T>
    Matrix<T> operator+(const Matrix<T> &t_m1, const Matrix<T> &t_m2)
    {
      Matrix<T> m = t_m1;
      m += t_m2;
      return m;
    }

    template<typename T>
    Matrix<T> operator-(const Matrix<T> &t_m1, const Matrix<T> &t_m2)
    {
      Matrix<T> m = t_m1;
      m -= t_m2;
      return m;
    }
  }
}
