namespace matrix { inline namespace v1 {
        template<EigenMethod em=EigenMethod::zheevd>
        std::tuple<Matrix<double>, Matrix<cxdbl>> eigenSys(const Matrix<cxdbl> &mat)
        {
#ifndef NDEBUG
          if(mat.nrows() != mat.ncols()){
            throw MatrixSizeMismatchError("number of rows need to be the same as number of columns.");
          }
          if(mat.nrows() < 1){
            throw MatrixSizeMismatchError("matrix need to have at least 1 rows.");
          }
#endif
          const size_t n = mat.nrows();
          constexpr char uplo = 'L';
          constexpr char jobz = 'V';
          const size_t lda = n > 0 ? n : 1;

          Matrix<cxdbl> a = mat;
          Matrix<double> eigen_vals(n, 1);

          if constexpr(em == EigenMethod::zheev){
            const size_t lwork = n > 0 ? (n+1)*n : 1;
            auto work = std::make_unique<cxdbl[]>(lwork);

            const size_t rwork_size = n > 0 ? 3*n-2 : 1;
            auto rwork = std::make_unique<double[]>(rwork_size);

            int info = 0;
            mat_zheev(jobz, uplo, n, a.data(), lda,
                eigen_vals.data(), work.get(), lwork, rwork.get(), &info);

            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          else if constexpr(em == EigenMethod::zheevd){
            const size_t lwork = 2 * n + n * n;
            auto work = std::make_unique<cxdbl[]>(lwork);

            const size_t lrwork = 1 + 5*n + 2*n*n;
            auto rwork = std::make_unique<double[]>(lrwork);

            const size_t liwork = 3 + 5*n;
            auto iwork = std::make_unique<int[]>(liwork);

            int info = 0;
            mat_zheevd(jobz, uplo, n, a.data(), lda,
                eigen_vals.data(), work.get(), lwork, rwork.get(), lrwork,
                iwork.get(), liwork, &info);
            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          else if constexpr(em == EigenMethod::zheevr){
            constexpr char range = 'A';
            constexpr double vl = -1; ///< not used
            constexpr double vu = 1; ///< not used
            constexpr size_t il = 1; ///< not used
            size_t iu = n; ///< not used
            constexpr double abstol = std::numeric_limits<double>::epsilon();
            const size_t m = n;
            const size_t ldz = m;
            const size_t z_dim = (m > 0) ? ldz * m : ldz * 1;
            auto z = std::make_unique<cxdbl[]>(z_dim);
            const size_t isuppz_size = (m > 0) ? 2 * m : 2 * 1;
            auto isuppz = std::make_unique<int[]>(isuppz_size);

            const size_t lwork = (n + 1) * n;
            auto work = std::make_unique<cxdbl[]>(lwork);

            const size_t lrwork = (m > 0) ? 24 * n : 1;
            auto rwork = std::make_unique<double[]>(lrwork);

            const size_t liwork = (m > 0) ? 10 * n : 1;
            auto iwork = std::make_unique<int[]>(liwork);

            int info = 0;
            mat_zheevr(jobz, range, uplo, n, a.data(), lda,
                vl, vu, il, iu, abstol, m,
                eigen_vals.data(), z.get(), ldz, isuppz.get(), 
                work.get(), lwork, rwork.get(), lrwork,
                iwork.get(), liwork, &info);
            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          return std::make_tuple(std::move(eigen_vals), std::move(a));
        }

        template<EigenMethod em=EigenMethod::zheevd>
        Matrix<double> eigenVal(const Matrix<cxdbl> &mat)
        {
          const size_t n = mat.nrows();
          constexpr char uplo = 'L';
          constexpr char jobz = 'N';
          const size_t lda = n > 0 ? n : 1;

          Matrix<cxdbl> a = mat;
          Matrix<double> eigen_vals(n, 1);

          if constexpr(em == EigenMethod::zheev){
            const size_t lwork = n > 0 ? (n+1)*n : 1;
            auto work = std::make_unique<cxdbl[]>(lwork);

            const size_t rwork_size = n > 0 ? 3*n-2 : 1;
            auto rwork = std::make_unique<double[]>(rwork_size);

            int info = 0;
            mat_zheev(jobz, uplo, mat.nrows(), a.data(), lda,
                eigen_vals.data(), work.get(), lwork, rwork.get(), &info);

            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          else if constexpr(em == EigenMethod::zheevd){
            const size_t lwork = n + 1;
            auto work = std::make_unique<cxdbl[]>(lwork);

            const size_t lrwork = n;
            auto rwork = std::make_unique<double[]>(lrwork);

            const size_t liwork = 1;
            auto iwork = std::make_unique<int[]>(liwork);

            int info = 0;
            mat_zheevd(jobz, uplo, n, a.data(), lda,
                eigen_vals.data(), work.get(), lwork, rwork.get(), lrwork,
                iwork.get(), liwork, &info);
            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          else if constexpr(em == EigenMethod::zheevr){
            constexpr char range = 'A';
            constexpr double vl = -1; ///< not used
            constexpr double vu = 1; ///< not used
            constexpr size_t il = 1; ///< not used
            size_t iu = n; ///< not used
            constexpr double abstol = std::numeric_limits<double>::epsilon();
            const size_t m = n;
            const size_t ldz = m;
            const size_t z_dim = (m > 0) ? ldz * m : ldz * 1;
            auto z = std::make_unique<cxdbl[]>(z_dim);
            const size_t isuppz_size = (m > 0) ? 2 * m : 2 * 1;
            auto isuppz = std::make_unique<int[]>(isuppz_size);

            const size_t lwork = (n + 1) * n;
            auto work = std::make_unique<cxdbl[]>(lwork);

            const size_t lrwork = (m > 0) ? 24 * n : 1;
            auto rwork = std::make_unique<double[]>(lrwork);

            const size_t liwork = (m > 0) ? 10 * n : 1;
            auto iwork = std::make_unique<int[]>(liwork);

            int info = 0;
            mat_zheevr(jobz, range, uplo, n, a.data(), lda,
                vl, vu, il, iu, abstol, m,
                eigen_vals.data(), z.get(), ldz, isuppz.get(), 
                work.get(), lwork, rwork.get(), lrwork,
                iwork.get(), liwork, &info);
            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          return eigen_vals;
        }   // eigenVal()

        template<typename T>
        Matrix<T> kroneckerProduct(const Matrix<T> &m1, const Matrix<T> &m2)
        {
          size_t nrows = m1.nrows() * m2.nrows();
          size_t ncols = m1.ncols() * m2.ncols();
          Matrix<T> res(nrows, ncols);

          for(size_t c1 = 0; c1 < m1.ncols(); ++c1){
            for(size_t r1 = 0; r1 < m1.nrows(); ++r1){
              for(size_t c2 = 0; c2 < m2.ncols(); ++c2){
                for(size_t r2 = 0; r2 < m2.nrows(); ++r2){
                  const size_t r = r1*m2.nrows() + r2;
                  const size_t c = c1*m2.ncols() + c2;
                  res(r,c) = m1(r1,c1) * m2(r2, c2);
                }
              }
            }
          }
          return res;
        }

        template<typename T>
        T trace(const Matrix<T> &m)
        {
          const size_t min_size = std::min(m.nrows(), m.ncols());
          T result = 0;
          for(size_t i = 0; i < min_size; ++i){
            result += m(i,i);
          }
          return result;
        }

        template<typename T>
        T projection(const Matrix<T> &lhs, const Matrix<T> &rhs)
        {
#ifndef NDEBUG
          if(lhs.nrows() != rhs.nrows() || lhs.ncols() != rhs.ncols())
            throw MatrixSizeMismatchError("matrices need to have the same size to project.");
          if(lhs.nrows() != lhs.ncols())
            throw MatrixSizeMismatchError("Projection needs to be done on a square matrix.");
#endif
          if constexpr(is_complex<T>::value){
            T res = lvl1_zdotc(lhs.nelements(), lhs.data(), 1, rhs.data(), 1);
            return res;
          } else if constexpr(is_double<T>::value){
            T res = lvl1_ddot(lhs.nelements(), lhs.data(), 1, rhs.data(), 1);
            return res;
          } else {
            T res = 0;
            for(size_t i = 0; i < lhs.ncols(); ++i){
              T temp = 0;
              for(size_t j = 0; j < lhs.nrows(); ++j){
                temp += lhs(j,i) * rhs(i,j);
              }
              res += temp;
            }
            return res;
          }
        }   // projection()

        template<typename T>
        T projectionNorm(const Matrix<T> &lhs, const Matrix<T> &rhs)
        {
          auto val = projection<T>(lhs, rhs);
          auto norm = projection<T>(rhs, rhs);
          return val/norm;
        }   // projectionNorm

        
        template<typename T>
        Matrix<T> pow(const Matrix<T> &mat, size_t n)
        {
#ifndef NDEBUG
          if(mat.nrows() != mat.ncols()){
            throw MatrixIsNotSquare("Matrix need to be square for pow().");
          }
#endif
          if (n == 0){
            return identity<T>(mat.nrows());
          } else if (n % 2 == 0) {
            auto res = pow<T>(mat, n/2) * pow<T>(mat, n/2);
            return res;
          } else if (n == 1){
            return mat;
          } else {
            auto res = mat * pow<T>(mat, n-1);
            return res;
          }
        }

        template<typename T>
        T norm1(const Matrix<T> &mat)
        {
          T val = 0;
          for(size_t c = 0; c < mat.ncols(); ++c){
            T col_sum = 0;
            for(size_t r = 0; r < mat.nrows(); ++r){
              col_sum += std::abs(mat(r, c));
            }
            if(col_sum > val){
              val = col_sum;
            }
          }
          return val;
        }
} // namespace v1
} // namespace matrix

