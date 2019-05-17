amespace practice { inline namespace v1 {
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
        : m_nrows(il.size())
    {
        auto r_it = il.begin();
        size_t ncols = 0;
        for(; r_it != il.end(); ++r_it){
            if(r_it->size() > ncols){
                ncols = r_it->size();
            }
        }   ///< find number of columns
        m_ncols = ncols;
        m_data = std::make_unique<T[]>(il.size() * ncols);

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
    Matrix<T>::Matrix(const Matrix<T> &rhs)
        : m_nrows(rhs.m_nrows), m_ncols(rhs.m_ncols),
        m_data(std::make_unique<T[]>(rhs.m_nrows * rhs.m_ncols))
    {
        for(size_t i = 0; i < m_nrows * m_ncols; ++i){
            m_data[i] = rhs.m_data[i];
        }
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
                m_data[i] = 1;
            }
        } else if constexpr(is_complex<T>::value){
            for(size_t i = 0; i < m_nrows * m_ncols; ++i){
                m_data[i] = T(0,0);
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
    size_t Matrix<T>::rows() const
    { return m_nrows; }

    template<typename T>
    size_t Matrix<T>::cols() const
    { return m_ncols; }

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
