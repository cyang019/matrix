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
        }
} // namespace v1
} // namespace matrix

