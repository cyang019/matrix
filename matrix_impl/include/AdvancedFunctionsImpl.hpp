namespace matrix { inline namespace v1 {
        template<EigenMethod em>
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
          return std::make_tuple(std::move(eigen_vals), std::move(a));
        }

        template<EigenMethod em>
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
          return eigen_vals;
        }
} // namespace v1
} // namespace matrix

