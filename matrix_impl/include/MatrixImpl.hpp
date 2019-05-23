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
      size_t offset = 0;
      const auto &beg = m_data.get();
      const double *raw_pos = raw_data;
      while(n_total > static_cast<size_t>(int_max)){
        cblas_dcopy(int_max, raw_pos, 1, std::next(beg, offset), 1);

        n_total -= (size_t)int_max;
        offset += (size_t)int_max;
        raw_pos += int_max;
      }
      cblas_dcopy(n_total, raw_pos, 1, std::next(beg, offset), 1);
    }

    template<>
    inline
    Matrix<cxdbl>::Matrix(size_t n_total, const cxdbl *raw_data)
    {
      size_t offset = 0;
      const auto &beg = m_data.get();
      const cxdbl *raw_pos = raw_data;
      while(n_total > static_cast<size_t>(int_max)){
        cblas_zcopy(int_max, raw_pos, 1, std::next(beg, offset), 1);

        n_total -= (size_t)int_max;
        offset += (size_t)int_max;
        raw_pos += int_max;
      }
      cblas_zcopy(n_total, raw_pos, 1, std::next(beg, offset), 1);
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
      size_t n_total = rhs.m_nrows * rhs.m_ncols;
      size_t offset = 0;
      const auto &beg = m_data.get();
      const auto &raw_beg = rhs.m_data.get();
      while(n_total > static_cast<size_t>(int_max)){
        cblas_dcopy(int_max, std::next(raw_beg, offset), 1, std::next(beg, offset), 1);

        n_total -= (size_t)int_max;
        offset += (size_t)int_max;
      }
      cblas_dcopy(int_max, std::next(raw_beg, offset), 1, std::next(beg, offset), 1);
    }

    template<>
    inline
    Matrix<cxdbl>::Matrix(const Matrix<cxdbl> &rhs)
        : m_nrows(rhs.m_nrows), m_ncols(rhs.m_ncols),
        m_data(std::make_unique<cxdbl[]>(rhs.m_nrows * rhs.m_ncols))
    {
      size_t n_total = rhs.m_nrows * rhs.m_ncols;
      size_t offset = 0;
      const auto &beg = m_data.get();
      const auto &raw_beg = rhs.m_data.get();
      while(n_total > static_cast<size_t>(int_max)){
        cblas_zcopy(int_max, std::next(raw_beg, offset), 1, std::next(beg, offset), 1);

        n_total -= (size_t)int_max;
        offset += (size_t)int_max;
      }
      cblas_zcopy(int_max, std::next(raw_beg, offset), 1, std::next(beg, offset), 1);
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
        if(t_i >= m_nrows || t_j >= m_nrows){
            throw std::out_of_range("row index out of bound.");
        }
        for(size_t idx = 0; idx < m_ncols; idx+=m_nrows){
            std::swap(m_data[t_i + idx * m_ncols], m_data[t_j + idx * m_ncols]);
        }

        return *this;
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator+=(const Matrix<T> &rhs)
    {
      if(rhs.m_nrows != m_nrows || rhs.m_ncols != m_ncols)
          throw std::length_error("cannot add matrices of different dimensions.");

      for(size_t i = 0; i < m_ncols * m_nrows; ++i){
          m_data[i] += rhs.m_data[i];
      }

      return *this;
    }

    template<>
    inline
    Matrix<double>& Matrix<double>::operator+=(const Matrix<double> &rhs)
    {
    }

    template<typename T>
    Matrix<T>& Matrix<T>::operator-=(const Matrix<T> &rhs)
    {
        if(rhs.m_nrows != m_nrows || rhs.m_ncols != m_ncols)
            throw std::length_error("cannot subtract matrices of different dimensions.");

        for(size_t i = 0; i < m_ncols * m_nrows; ++i){
            m_data[i] -= rhs.m_data[i];
        }

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
              m_data[i] = T(0.0,0.0);
          }
      }else {
          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = T();
          }
      }

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
              m_data[i] = T(1.0,0);
          }
      } else {
          for(size_t i = 0; i < m_nrows * m_ncols; ++i){
              m_data[i] = T();
          }
      }

      return *this;
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
        const size_t idx = m_nrows * t_col + t_row;
        return m_data[idx];
    }

    template<typename T>
    T& Matrix<T>::operator()(size_t t_row, size_t t_col)
    {
      const size_t idx = m_nrows * t_col + t_row;
      return m_data[idx];
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


} ///< inline namespace v1
} ///< namespace matrix
