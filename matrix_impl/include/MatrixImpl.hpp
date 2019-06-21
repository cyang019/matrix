namespace matrix { inline namespace v1 {
    template<typename T>
    Matrix<T>::Matrix()
        : m_nrows(0), m_ncols(0), m_data(nullptr)
    {}

    template<typename T>
    Matrix<T>::Matrix(size_t nrows, size_t ncols)
        : m_nrows(nrows), m_ncols(ncols),
        m_data(std::make_unique<T[]>(nrows * ncols))
    {}

    template<typename T>
    Matrix<T>::Matrix(size_t n)
        : m_nrows(n), m_ncols(n),
        m_data(std::make_unique<T[]>(n*n))
    {}

    template<typename T>
    Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> il)
        : m_nrows(il.size()), m_ncols(0)
    {
      auto r_it = il.begin();
      for(; r_it != il.end(); ++r_it){
          if(m_ncols == 0){
            m_ncols = r_it->size();
          }
          else if(r_it->size() != m_ncols){
            throw(MatrixColsSizeUnevenError("Initializer size error."));
          }
      }   ///< find number of columns
      m_data = std::make_unique<T[]>(m_nrows * m_ncols);

      size_t at_row = 0;
      r_it = il.begin();
      while(r_it!=il.end()){
          size_t at_col = 0;
          auto c_it = r_it->begin(); 
          while(c_it != r_it->end()){
              const size_t idx = at_col * m_nrows + at_row;
              
              //if constexpr (is_complex<T>::value) {
              //  std::cout << "index: " << idx << " " 
              //            << "r:" << at_row << " c:" << at_col << "\n";
              //}

              m_data[idx] = *c_it;
              ++c_it;
              ++at_col;
          }
          ++r_it;
          ++at_row;
      }
    }

    template<typename T>
    Matrix<T>::Matrix(size_t n_total, const T *raw_data)
    {
      for(size_t i = 0; i < n_total; ++i){
        m_data[i] = raw_data[i];
      }
    }

    template<>
    inline
    Matrix<double>::Matrix(size_t n_total, const double *raw_data)
    {
      lvl1_dcopy(n_total, raw_data, 1, m_data.get(), 1);
    }

    template<>
    inline
    Matrix<cxdbl>::Matrix(size_t n_total, const cxdbl *raw_data)
    {
      lvl1_zcopy(n_total, raw_data, 1, m_data.get(), 1);
    }

    template<typename T>
    Matrix<T>::Matrix(const Matrix<T> &rhs)
        : m_nrows(rhs.m_nrows), m_ncols(rhs.m_ncols),
        m_data(std::make_unique<T[]>(rhs.m_nrows * rhs.m_ncols))
    {
      for(size_t i = 0; i < m_nrows * m_ncols; ++i){
          m_data[i] = rhs.m_data[i];
      }
    }

    template<>
    inline
    Matrix<double>::Matrix(const Matrix<double> &rhs)
        : m_nrows(rhs.m_nrows), m_ncols(rhs.m_ncols),
        m_data(std::make_unique<double[]>(rhs.m_nrows * rhs.m_ncols))
    {
      const size_t n_total = rhs.m_nrows * rhs.m_ncols;
      lvl1_dcopy(n_total, rhs.m_data.get(), 1, m_data.get(), 1);
    }

    template<>
    inline
    Matrix<cxdbl>::Matrix(const Matrix<cxdbl> &rhs)
        : m_nrows(rhs.m_nrows), m_ncols(rhs.m_ncols),
        m_data(std::make_unique<cxdbl[]>(rhs.m_nrows * rhs.m_ncols))
    {
      const size_t n_total = rhs.m_nrows * rhs.m_ncols;
      lvl1_zcopy(n_total, rhs.m_data.get(), 1, m_data.get(), 1);
    }

    template<typename T>
    Matrix<T>::Matrix(Matrix<T> &&rhs) noexcept
        : m_nrows(rhs.m_nrows), m_ncols(rhs.m_ncols),
        m_data(std::move(rhs.m_data))
    {}

    template<typename T>
    Matrix<T>& Matrix<T>::operator=(const Matrix<T> &rhs)
    {
      m_nrows = rhs.m_nrows;
      m_ncols = rhs.m_ncols;
      m_data = std::make_unique<T[]>(rhs.m_nrows * rhs.m_ncols);
      for(size_t i = 0; i < m_nrows * m_ncols; ++i){
          m_data[i] = rhs.m_data[i];
      }
      return *this;
    }

    template<>
    inline
    Matrix<double>& Matrix<double>::operator=(const Matrix<double> &rhs)
    {
      m_nrows = rhs.m_nrows;
      m_ncols = rhs.m_ncols;
      m_data = std::make_unique<double[]>(rhs.m_nrows * rhs.m_ncols);
      lvl1_dcopy(m_nrows * m_ncols, rhs.m_data.get(), 1, m_data.get(), 1);
      return *this;
    }

    template<>
    inline
    Matrix<cxdbl>& Matrix<cxdbl>::operator=(const Matrix<cxdbl> &rhs)
    {
      m_nrows = rhs.m_nrows;
      m_ncols = rhs.m_ncols;
      m_data = std::make_unique<cxdbl[]>(rhs.m_nrows * rhs.m_ncols);
      lvl1_zcopy(m_nrows * m_ncols, rhs.m_data.get(), 1, m_data.get(), 1);
      return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator=(Matrix<T> &&rhs) noexcept
    {
        m_nrows = rhs.m_nrows;
        m_ncols = rhs.m_ncols;
        m_data = std::move(rhs.m_data);
        return *this;
    }

    template<typename T>
    Matrix<T>::~Matrix()
    {}

    /// column major ordering
    template<typename T>
    Matrix<T>& Matrix<T>::swapCols(size_t t_i, size_t t_j)
    {
        if(t_i >= m_ncols || t_j >= m_ncols){
            throw std::out_of_range("column index out of bound.");
        }
        for(size_t idx = 0; idx < m_nrows; ++idx){
            std::swap(m_data[t_i * m_nrows + idx], m_data[t_j * m_nrows + idx]);
        }

        return *this;
    }

    /// column major ordering
    template<typename T>
    Matrix<T>& Matrix<T>::swapRows(size_t t_i, size_t t_j)
    {
#ifndef NDEBUG
      if(t_i >= m_nrows || t_j >= m_nrows){
          throw IndexOutOfBound("row index out of bound.");
      }
#endif
      for(size_t idx = 0; idx < m_ncols; idx+=m_nrows){
          std::swap(m_data[t_i + idx * m_ncols], m_data[t_j + idx * m_ncols]);
      }

      return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator+=(const Matrix<T> &rhs)
    {
#ifndef NDEBUG
      if(rhs.m_nrows != m_nrows || rhs.m_ncols != m_ncols)
          throw MatrixSizeMismatchError("cannot add matrices of different dimensions.");
#endif
      for(size_t i = 0; i < m_ncols * m_nrows; ++i){
          m_data[i] += rhs.m_data[i];
      }

      return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator+=(const T &rhs)
    {
      for(size_t i = 0; i < m_ncols * m_nrows; ++i){
          m_data[i] += rhs;
      }

      return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator-=(const T &rhs)
    {
      for(size_t i = 0; i < m_ncols * m_nrows; ++i){
          m_data[i] -= rhs;
      }

      return *this;
    }

    template<>
    inline
    Matrix<double>& Matrix<double>::operator+=(const Matrix<double> &rhs)
    {
#ifndef NDEBUG
      if(rhs.m_nrows != m_nrows || rhs.m_ncols != m_ncols)
          throw MatrixSizeMismatchError("cannot add matrices of different dimensions.");
#endif

      lvl1_daxpy(m_ncols * m_nrows, 1.0, rhs.m_data.get(), 1, m_data.get(), 1);
      return *this;
    }

    template<>
    inline
    Matrix<cxdbl>& Matrix<cxdbl>::operator+=(const Matrix<cxdbl> &rhs)
    {
#ifndef NDEBUG
      if(rhs.m_nrows != m_nrows || rhs.m_ncols != m_ncols)
          throw MatrixSizeMismatchError("cannot add matrices of different dimensions.");
#endif
      assert(rhs.m_nrows == m_nrows && rhs.m_ncols == m_ncols);
      lvl1_zaxpy(m_ncols * m_nrows, 1.0, rhs.m_data.get(), 1, m_data.get(), 1);
      return *this;
    }
    
    template<typename T>
    Matrix<T>& Matrix<T>::operator-=(const Matrix<T> &rhs)
    {
#ifndef NDEBUG
      if(rhs.m_nrows != m_nrows || rhs.m_ncols != m_ncols)
          throw MatrixSizeMismatchError("cannot add matrices of different dimensions.");
#endif

      for(size_t i = 0; i < m_ncols * m_nrows; ++i){
          m_data[i] -= rhs.m_data[i];
      }

      return *this;
    }

    template<>
    inline
    Matrix<double>& Matrix<double>::operator-=(const Matrix<double> &rhs)
    {
#ifndef NDEBUG
      if(rhs.m_nrows != m_nrows || rhs.m_ncols != m_ncols)
          throw MatrixSizeMismatchError("cannot add matrices of different dimensions.");
#endif
      lvl1_daxpy(m_ncols * m_nrows, -1.0, rhs.m_data.get(), 1, m_data.get(), 1);
      return *this;
    }

    template<>
    inline
    Matrix<cxdbl>& Matrix<cxdbl>::operator-=(const Matrix<cxdbl> &rhs)
    {
#ifndef NDEBUG
      if(rhs.m_nrows != m_nrows || rhs.m_ncols != m_ncols)
          throw MatrixSizeMismatchError("cannot add matrices of different dimensions.");
#endif
      lvl1_zaxpy(m_ncols * m_nrows, -1.0, rhs.m_data.get(), 1, m_data.get(), 1);
      return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator*=(const T &rhs)
    {
      for(size_t i = 0; i < m_ncols * m_nrows; ++i){
        m_data[i] *= rhs;
      }
      return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator*=(const Matrix<T> &rhs)
    {
#ifndef NDEBUG
      if(m_ncols != rhs.m_nrows){
        throw MultiplicationError(
            "A of size a1 x a2 times B of size b1 x b2: a2 should equal b1.");
      }
#endif
      auto new_data = std::make_unique<T[]>(m_nrows * rhs.m_ncols);
      for(size_t c2 = 0; c2 < rhs.m_ncols; ++c2){
        for(size_t r1 = 0; r1 < m_nrows; ++r1){
          // the dot product to produce (r1, c2) element
          new_data[c2 * m_nrows + r1] = 0;
          for(size_t c1 = 0; c1 < m_ncols; ++c1){
            size_t pos2 = rhs.m_nrows * c2;
            for(size_t r2 = 0; r2 < rhs.m_nrows; ++r2) {  ///< column major layout
              new_data[c2 * m_nrows + r1] += m_data[c1 * m_nrows + r1] * rhs.m_data[pos2];
              ++pos2;
            }
          }
        }
      }
      m_data = std::move(new_data);
      return *this;
    }

    template<>
    inline
    Matrix<double>& Matrix<double>::operator*=(const Matrix<double> &rhs)
    {
#ifndef NDEBUG
      if(m_ncols != rhs.m_nrows){
        throw MultiplicationError(
            "A of size a1 x a2 times B of size b1 x b2: a2 should equal b1.");
      }
#endif
      auto new_data = std::make_unique<double[]>(m_nrows * rhs.m_ncols);
      lvl3_dgemm(CblasOrder::CblasColMajor,
          CblasTranspose::CblasNoTrans,
          CblasTranspose::CblasNoTrans,
          m_nrows, rhs.m_ncols, m_ncols,
          1.0, m_data.get(), m_nrows, rhs.data(), rhs.m_nrows,
          0.0, new_data.get(), m_nrows);

      m_data = std::move(new_data);
      return *this;
    }


    template<>
    inline
    Matrix<double>& Matrix<double>::operator*=(const double &rhs)
    {
      lvl1_dscal(m_ncols * m_nrows, rhs, m_data.get(), 1);
      return *this;
    }

    template<>
    inline
    Matrix<cxdbl>& Matrix<cxdbl>::operator*=(const cxdbl &rhs)
    {
      lvl1_zscal(m_ncols * m_nrows, rhs, m_data.get(), 1);
      return *this;
    }

    inline
    Matrix<cxdbl>& operator*=(Matrix<cxdbl> &lhs, const double &rhs)
    {
      lvl1_zdscal(lhs.ncols() * lhs.nrows(), rhs, lhs.data(), 1);
      return lhs;
    }


    template<typename T>
    Matrix<T>& Matrix<T>::operator/=(const T &rhs)
    {
      for(size_t i = 0; i < m_ncols * m_nrows; ++i){
        m_data[i] /= rhs;
      }
      return *this;
    }

    template<>
    inline
    Matrix<double>& Matrix<double>::operator/=(const double &rhs)
    {
      lvl1_dscal(m_nrows * m_ncols, 1.0/rhs, m_data.get(), 1);
      return *this;
    }

    template<>
    inline
    Matrix<cxdbl>& Matrix<cxdbl>::operator/=(const cxdbl &rhs)
    {
      lvl1_zscal(m_nrows * m_ncols, 1.0/rhs, m_data.get(), 1);
      return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::setZero()
    {
      if constexpr (std::is_floating_point<T>::value
              || std::is_integral<T>::value){
          // real type
          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = 0;
          }
      } else if constexpr (is_complex<T>::value) {
          // complex type
          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = cx_zero;
          }
      }else {
          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = T();
          }
      }

      return *this;
    }

    template<>
    inline
    Matrix<double>& Matrix<double>::setZero()
    {
      lvl1_dscal(m_ncols*m_nrows, 0.0, m_data.get(), 1);
      return *this;
    }

    template<>
    inline
    Matrix<cxdbl>& Matrix<cxdbl>::setZero()
    {
      lvl1_zscal(m_ncols*m_nrows, cx_zero, m_data.get(), 1);
      return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::setOne()
    {
      if constexpr (std::is_floating_point<T>::value
              || std::is_integral<T>::value){
          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = 1.0;
          }
      } else if constexpr(is_complex<T>::value){
          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = T(1, 0);
          }
      } else {
          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = T();
          }
      }

      return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::row(size_t r)
    {
#ifndef NDEBUG
      if(r >= m_nrows)
        throw IndexOutOfBound("row number out of bound.");
#endif
      Matrix<T> res(1, m_ncols);
      for(size_t i = 0; i < m_ncols; ++i){
        res.m_data[i] = m_data[i*m_nrows + r]; 
      }
      return res;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::col(size_t c)
    {
#ifndef NDEBUG
      if(c >= m_ncols)
        throw IndexOutOfBound("column number out of bound.");
#endif
      Matrix<T> res(m_nrows, 1);
      for(size_t i = 0; i < m_nrows; ++i){
        res.m_data[i] = m_data[c*m_nrows + i]; 
      }
      return res;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::setRandom()
    {
      std::random_device rd;
      std::mt19937 gen(rd());
      if constexpr(std::is_integral<T>::value){
          std::uniform_int_distribution<> dist(1, 10);

          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = dist(gen);
          }
          return *this;
      } else if constexpr(std::is_floating_point<T>::value){
          std::uniform_real_distribution<> dist(0.0, 1.0);

          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = dist(gen);
          }
          return *this;
      } else if constexpr(is_complex<T>::value){
          std::uniform_real_distribution<> dist(0.0, 1.0);

          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i].real() = dist(gen);
              m_data[i].imag() = dist(gen);
          }
          return *this;
      } else{
          throw std::logic_error("non numeric type.");
      }
    }

    template<typename T>
    T& Matrix<T>::operator()(size_t t_row, size_t t_col) const
    {
#ifndef NDEBUG
      if (t_row >= m_nrows || t_col >= m_ncols)
        throw IndexOutOfBound("index out of bound.");
#endif
      const size_t idx = m_nrows * t_col + t_row;
      return m_data[idx];
    }

    template<typename T>
    T& Matrix<T>::operator()(size_t t_row, size_t t_col)
    {
#ifndef NDEBUG
      if (t_row >= m_nrows || t_col >= m_ncols)
        throw IndexOutOfBound("index out of bound.");
#endif
      const size_t idx = m_nrows * t_col + t_row;
      return m_data[idx];
    }

    template<typename T>
    Matrix<T> Matrix<T>::t() const
    {
      auto res = Matrix<T>(m_ncols, m_nrows);
      for(size_t i = 0; i < m_ncols; ++i){
        for(size_t j = 0; j < m_nrows; ++j){
          const auto orig_pos = i * m_nrows + j;
          const auto dest_pos = j * m_ncols + i;
          res.m_data[dest_pos] = m_data[orig_pos];
        }
      }
      return res;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::tInplace()
    {
      if(m_nrows != 1 && m_ncols != 1) {
        auto transposed = std::make_unique<double[]>(m_ncols * m_nrows);
        for(size_t i = 0; i < m_ncols; ++i){
          for(size_t j = 0; j < m_nrows; ++j){
            const auto orig_pos = i * m_nrows + j;
            const auto dest_pos = j * m_ncols + i;
            transposed[dest_pos] = m_data[orig_pos];
          }
        }
        m_data = std::move(transposed);
      }
      std::swap(m_nrows, m_ncols);
      return *this;
    }

    template<typename T>
    Matrix<T> Matrix<T>::adjoint() const
    {
      if constexpr (is_complex<T>::value){
        auto res = Matrix<T>(m_ncols, m_nrows);
        for(size_t i = 0; i < m_ncols; ++i){
          for(size_t j = 0; j < m_nrows; ++j){
            const auto orig_pos = i * m_nrows + j;
            const auto dest_pos = j * m_ncols + i;
            res.m_data[dest_pos] = std::conj(m_data[orig_pos]);
          }
        }
        return res;
      }
      else{
        return Matrix<T>::t();
      }
    }

    template<typename T>
    Matrix<T>& Matrix<T>::adjointInplace()
    {
      if constexpr (is_complex<T>::value){
        auto transposed = std::make_unique<cxdbl[]>(m_ncols * m_nrows);
        for(size_t i = 0; i < m_ncols; ++i){
          for(size_t j = 0; j < m_nrows; ++j){
            const auto orig_pos = i * m_nrows + j;
            const auto dest_pos = j * m_ncols + i;
            transposed[dest_pos] = std::conj(m_data[orig_pos]);
          }
        }
        std::swap(m_nrows, m_ncols);
        m_data = std::move(transposed);
        return *this;
      }
      else{
        return Matrix<T>::tInplace();
      }
    }

    template<typename T> 
    Matrix<T> Matrix<T>::inverse() const
    {
      throw NotImplementedError("Matrix Inverse not implemented.");
    }

    template<>
    inline
    Matrix<double> Matrix<double>::inverse() const
    {
#ifndef NDEBUG
      if(m_nrows != m_ncols){
        throw MatrixSizeMismatchError("Row and column numbers need to be the same.");
      }
      if(m_nrows == 0){
        throw MatrixSizeMismatchError("Row and column numbers cannot be zero.");
      }
#endif
      auto ptr_ipiv = std::make_unique<int[]>(m_nrows);
      int errorHandler;
      auto ptr_lwork = std::make_unique<double[]>(
          m_nrows * m_ncols
          );

      Matrix<double> result = *this;
      mat_dgetrf(m_nrows, m_ncols, result.m_data.get(), m_nrows, ptr_ipiv.get(), &errorHandler);
      if(errorHandler > 0){
        throw NonInvertibleMatrix("Matrix non-invertible.");
      }

      mat_dgetri(m_nrows, result.m_data.get(), m_nrows, ptr_ipiv.get(), ptr_lwork.get(), m_nrows * 2, &errorHandler);

      return result;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::inverseInplace()
    {
      throw NotImplementedError("Matrix Inverse not implemented.");
    }

    template<>
    inline
    Matrix<double>& Matrix<double>::inverseInplace()
    {
#ifndef NDEBUG
      if(m_nrows != m_ncols){
        throw MatrixSizeMismatchError("Row and column numbers need to be the same.");
      }
      if(m_nrows == 0){
        throw MatrixSizeMismatchError("Row and column numbers cannot be zero.");
      }
#endif
      auto ptr_ipiv = std::make_unique<int[]>(m_nrows);
      int errorHandler;
      auto ptr_lwork = std::make_unique<double[]>(
          m_nrows * m_ncols);

      mat_dgetrf(m_nrows, m_ncols, m_data.get(), m_nrows, ptr_ipiv.get(), &errorHandler);
      if(errorHandler > 0){
        throw NonInvertibleMatrix("Matrix non-invertible.");
      }

      mat_dgetri(m_nrows, m_data.get(), m_nrows, ptr_ipiv.get(), ptr_lwork.get(), m_nrows * 2, &errorHandler);

      return *this;
    }
  
    template<typename T>
    T* Matrix<T>::data()
    { return m_data.get(); }

    template<typename T>
    const T* Matrix<T>::data() const
    { return m_data.get(); }

    template<typename T>
    size_t Matrix<T>::nrows() const
    { return m_nrows; }

    template<typename T>
    size_t Matrix<T>::ncols() const
    { return m_ncols; }

    template<typename T>
    T Matrix<T>::trace() const
    {
      T result = 0;
      const size_t n_total = std::min(m_nrows, m_ncols);
      for(size_t i = 0; i < n_total; ++i){
        const size_t pos = i * m_nrows + i;
        result += m_data[pos];
      }
      return result;
    }

    template<typename T>
    std::tuple<size_t, size_t> Matrix<T>::shape() const
    { return std::make_tuple(m_nrows, m_ncols); }

    template<typename T>
    Matrix<T> diagnal(std::initializer_list<T> vals)
    {
      Matrix<T> res(vals.size(), vals.size());
      auto iter = vals.begin();
      for(size_t i = 0; i < vals.size(); ++i, ++iter){
        for(size_t j = 0; j < vals.size(); ++j){
          if(i != j){
            res(j, i) = 0;
          }
          else{
            res(i,i) = *iter;
          }
        }
      }
      return res;
    }

    template<typename T>
    Matrix<T> identity(size_t n)
    {
      auto res = Matrix<T>(n, n);
      for(size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j){
          if(i != j){
            res(j, i) = 0;
          } else {
            res(i, j) = 1;
          }
        }
      }
      return res;
    }

    template<typename T>
    Matrix<T> zeros(size_t nrows, size_t ncols)
    {
      auto res = Matrix<T>(nrows, ncols);
      for(size_t i = 0; i < ncols; ++i){
        for(size_t j = 0; j < nrows; ++j){
          res(j, i) = 0;
        }
      }
      return res;
    }

    template<typename T>
    Matrix<T> ones(size_t nrows, size_t ncols)
    {
      auto res = Matrix<T>(nrows, ncols);
      for(size_t i = 0; i < ncols; ++i){
        for(size_t j = 0; j < nrows; ++j){
          res(j, i) = 1;
        }
      }
      return res;
    }

    template<typename T>
    Matrix<T> sign(const Matrix<T> &mat)
    {
      if constexpr(is_complex<T>::value){
        auto res = mat;
        auto ptr_res = res.data();
        for(size_t i = 0; i < res.nelements(); ++i){
          const auto abs_val = std::abs(*ptr_res);
          if(abs_val < eps){
            *ptr_res = 1;
          } else {
            *ptr_res /= abs_val;
          }
          ++ptr_res;
        }
        return res;
      } else {
        auto res = Matrix<T>(mat.nrows(), mat.ncols());
        auto ptr_res = res.data();
        auto ptr_mat = mat.data();
        for(size_t i = 0; i < res.nelements(); ++i){
          if(*ptr_mat > -eps) *ptr_res = 1;
          else *ptr_res = -1;
          ++ptr_res;
          ++ptr_mat;
        }

        return res;
      }
    }
} ///< inline namespace v1
} ///< namespace matrix
