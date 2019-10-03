namespace matrix {
  inline namespace v1 {
    template<typename T>
    std::ostream& operator<<(std::ostream &os, const Matrix<T> &m)
    {
      for(size_t i = 0; i < m.nrows(); ++i){
        for(size_t j = 0; j < m.ncols(); ++j){
          os << m(i,j) << "\t";
        }
        os << "\n";
      }
      return os;
    }

    template<typename T>
    bool operator==(const Matrix<T> &lhs, const Matrix<T> &rhs)
    {
      if(lhs.m_ncols != rhs.m_ncols || lhs.m_nrows != rhs.m_nrows)
        return false;

      const size_t n_total = lhs.m_ncols * lhs.m_nrows;
      for(size_t i = 0; i < n_total; ++i){
        if(lhs.m_data[i] != rhs.m_data[i]){
          return false;
        }
      }

      return true;
    }

    template<typename T>
    bool operator!=(const Matrix<T> &lhs, const Matrix<T> &rhs)
    {
      return !operator==<T>(lhs, rhs);
    }

    inline bool allclose(const Matrix<double> &lhs, const Matrix<double> &rhs, double eps)
    {
      if(lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
        return false;

      const size_t n_total = lhs.ncols() * lhs.nrows();
      auto d1 = lhs.data();
      auto d2 = rhs.data();
      for(size_t i = 0; i < n_total; ++i){
        if(std::abs(d1[i] - d2[i]) > eps){
          return false;
        }
      }

      return true;
    }

    inline bool allclose(const Matrix<cxdbl> &lhs, const Matrix<cxdbl> &rhs, double eps)
    {
      if(lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
        return false;

      const size_t n_total = lhs.ncols() * lhs.nrows();
      auto d1 = lhs.data();
      auto d2 = rhs.data();
      for(size_t i = 0; i < n_total; ++i){
        if(std::abs(d1[i].real() - d2[i].real()) > eps){
          return false;
        }
        if(std::abs(d1[i].imag() - d2[i].imag()) > eps){
          return false;
        }
      }

      return true;
    }

    template<typename T>
    Matrix<T> operator+(const Matrix<T> &t_m1, const Matrix<T> &t_m2)
    {
      Matrix<T> m = t_m1;
      m += t_m2;
      return m;
    }

    template<typename T>
    Matrix<T> operator+(const Matrix<T> &t_m1, const T &t_val2)
    {
      Matrix<T> m = t_m1;
      m += t_val2;
      return m;
    }

    template<typename T>
    Matrix<T> operator+(const T &t_val1, const Matrix<T> &t_m2)
    {
      Matrix<T> m = t_m2;
      m += t_val1;
      return m;
    }

    template<typename T>
    Matrix<T> operator-(const Matrix<T> &t_m1, const Matrix<T> &t_m2)
    {
      Matrix<T> m = t_m1;
      m -= t_m2;
      return m;
    }

    template<typename T>
    Matrix<T> operator-(const Matrix<T> &t_m1, const T &t_val2)
    {
      Matrix<T> m = t_m1;
      m -= t_val2;
      return m;
    }

    template<typename T>
    Matrix<T> operator-(const T &t_val1, const Matrix<T> &t_m2)
    {
      Matrix<T> m = t_val1 * identity<T>(t_m2.nrows(), t_m2.ncols());
      m -= t_m2;
      return m;
    }

    template<typename T>
    Matrix<T> operator-(const Matrix<T> &mat)
    {
      auto res = -1.0 * mat;
      return res;
    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T> &lhs, const Matrix<T> &rhs)
    {
#ifndef NDEBUG
      if(lhs.m_ncols != rhs.m_nrows){
        throw MultiplicationError(
            "A of size a1 x a2 times B of size b1 x b2: a2 should equal b1.");
      }
#endif
      Matrix<T> res(lhs.m_nrows, rhs.m_ncols);
      for(size_t c2 = 0; c2 < rhs.m_ncols; ++c2){
        for(size_t r1 = 0; r1 < lhs.m_nrows; ++r1){
          // the dot product to produce (r1, c2) element
          res(r1, c2) = 0;
          for(size_t c1 = 0; c1 < lhs.m_ncols; ++c1){
            size_t pos2 = rhs.m_nrows * c2;
            for(size_t r2 = 0; r2 < rhs.m_nrows; ++r2) {  ///< column major layout
              res(r1, c2) += lhs(r1, c1) * rhs.m_data[pos2++];
            }
          }
        }
      }
      return res;
    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T> &lhs, const T &rhs)
    {
      auto res = lhs;
      res *= rhs;
      return res;
    }

    template<typename T>
    Matrix<T> operator*(Matrix<T> &&lhs, const T &rhs)
    {
      auto res = std::move(lhs);
      res *= rhs;
      return res;
    }

    template<typename T>
    Matrix<T> operator*(const T &lhs, const Matrix<T> &rhs)
    {
      auto res = rhs;
      res *= lhs;
      return res;
    }

    template<typename T>
    Matrix<T> operator*(const T &lhs, Matrix<T> &&rhs)
    {
      auto res = std::move(rhs);
      res *= lhs;
      return res;
    }

    template<>
    inline
    Matrix<double> operator*(const Matrix<double> &lhs, const Matrix<double> &rhs)
    {
#ifndef NDEBUG
      if(lhs.m_ncols != rhs.m_nrows){
        throw MultiplicationError(
            "A of size a1 x a2 times B of size b1 x b2: a2 should equal b1.");
      }
#endif
      Matrix<double> res(lhs.m_nrows, rhs.m_ncols);
      if(lhs.m_nrows == 1 && rhs.m_ncols == 1) {
        res.m_data[0] = lvl1_ddot(lhs.m_ncols, lhs.m_data.get(), 1, rhs.m_data.get(), 1);
      } else if(rhs.m_ncols == 1) {
        lvl2_dgemv(CblasOrder::CblasColMajor, CblasTranspose::CblasNoTrans,
            lhs.m_nrows, lhs.m_ncols, 
            1.0, lhs.m_data.get(), lhs.m_nrows, rhs.m_data.get(), 1,
            0.0, res.m_data.get(), 1);
      } else if(lhs.m_nrows == 1) {
        lvl2_dgemv(CblasOrder::CblasColMajor, CblasTranspose::CblasTrans,
            rhs.m_nrows, rhs.m_ncols,
            1.0, rhs.m_data.get(), rhs.m_nrows, lhs.m_data.get(), 1,
            0.0, res.m_data.get(), 1);
      } else {
        lvl3_dgemm(CblasOrder::CblasColMajor, CblasTranspose::CblasNoTrans, CblasTranspose::CblasNoTrans,
            lhs.m_nrows, rhs.m_ncols, lhs.m_ncols,
            1.0, lhs.m_data.get(), lhs.m_nrows, rhs.m_data.get(), rhs.m_nrows,
            0.0, res.m_data.get(), res.m_nrows);
      }
      return res;
    }

    template<>
    inline
    Matrix<cxdbl> operator*(const Matrix<cxdbl> &lhs, const Matrix<cxdbl> &rhs)
    {
#ifndef NDEBUG
      if(lhs.m_ncols != rhs.m_nrows){
        throw MultiplicationError(
            "A of size a1 x a2 times B of size b1 x b2: a2 should equal b1.");
      }
#endif
      Matrix<cxdbl> res(lhs.m_nrows, rhs.m_ncols);
      if(rhs.m_ncols == 1) {
        lvl2_zgemv(CblasOrder::CblasColMajor, CblasTranspose::CblasNoTrans,
            lhs.m_nrows, lhs.m_ncols, 
            cxdbl(1.0, 0.0), lhs.m_data.get(), lhs.m_nrows, rhs.m_data.get(), 1,
            cxdbl(0.0, 0.0), res.m_data.get(), 1);
      } else if(lhs.m_nrows == 1) {
        lvl2_zgemv(CblasOrder::CblasColMajor, CblasTranspose::CblasConjTrans,
            rhs.m_nrows, rhs.m_ncols,
            cxdbl(1.0, 0.0), rhs.m_data.get(), rhs.m_nrows, lhs.m_data.get(), 1,
            cxdbl(0.0, 0.0), res.m_data.get(), 1);
      } else {
        lvl3_zgemm(CblasOrder::CblasColMajor, CblasTranspose::CblasNoTrans, CblasTranspose::CblasNoTrans,
            lhs.m_nrows, rhs.m_ncols, lhs.m_ncols,
            cxdbl(1.0, 0.0), lhs.m_data.get(), lhs.m_nrows, rhs.m_data.get(), rhs.m_nrows,
            cxdbl(0.0, 0.0), res.m_data.get(), res.m_nrows);
      }
      return res;
    }

    template<typename T>
    Matrix<T> flatten(const Matrix<T> &mat, char c)
    {
      if(c != 'C' && c != 'c' && c != 'R' && c != 'r'){
        throw IllegalValue("flatten() second parameter can only be c or r.");
      }
      if(c == 'C' || c == 'c'){
        Matrix<T> res(mat.m_nrows * mat.m_ncols, 1);
        
        size_t idx = 0;
        for(size_t i = 0; i < mat.m_nrows; ++i){
          for(size_t j = 0; j < mat.m_ncols; ++j){
            res.m_data[idx++] = mat(i, j);
          }
        }
        return res;
      } else {
        Matrix<T> res(1, mat.m_nrows * mat.m_ncols);
        
        size_t idx = 0;
        for(size_t i = 0; i < mat.m_nrows; ++i){
          for(size_t j = 0; j < mat.m_ncols; ++j){
            res.m_data[idx++] = mat(i, j);
          }
        }
        return res;
      }
    }
  } // namespace v1
} // namespace matrix
