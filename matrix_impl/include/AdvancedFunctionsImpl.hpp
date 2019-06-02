namespace matrix { inline namespace v1 {
        template<EigenMethod em>
        std::tuple<Matrix<double>, Matrix<cxdbl>> eigenSys(const Matrix<cxdbl> &mat)
        {
#ifndef NDEBUG
          if(mat.nrows() != mat.ncols()){
            throw MatrixSizeMismatchError("number of rows need to be the same as number of columns.");
          }
#endif
          const size_t n = mat.nrows();
          Matrix<double> eigen_vals(n, 1);
          Matrix<cxdbl> eigen_vecs = mat;
          if constexpr(em == EigenMethod::zheev){
            const size_t lda = n > 0 ? n : 1;
            const size_t lwork = n > 0 ? (n+1)*n : 1;
            auto work = std::make_unique<cxdbl[]>(lwork);

            const size_t rwork_size = n > 0 ? 3*n-2 : 1;
            auto rwork = std::make_unique<double[]>(rwork_size);

            int info = 0;
            mat_zheev('V', 'U', mat.nrows(), eigen_vecs.data(), lda,
                eigen_vals.data(), work.get(), lwork, rwork.get(), &info);

            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }

          }
          return std::make_tuple(std::move(eigen_vals), std::move(eigen_vecs));
        }

        template<EigenMethod em>
        Matrix<double> eigenVal(const Matrix<cxdbl> &mat)
        {
          const size_t n = mat.nrows();
          Matrix<double> vals(n, 1);
          Matrix<cxdbl> vecs = mat;
          if constexpr(em == EigenMethod::zheev){
            const size_t lda = n > 0 ? n : 1;
            const size_t lwork = n > 0 ? (n+1)*n : 1;
            auto work = std::make_unique<cxdbl[]>(lwork);

            const size_t rwork_size = n > 0 ? 3*n-2 : 1;
            auto rwork = std::make_unique<double[]>(rwork_size);

            int info = 0;
            mat_zheev('N', 'U', mat.nrows(), vecs.data(), lda,
                vals.data(), work.get(), lwork, rwork.get(), &info);

            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          return vals;
        }
} // namespace v1
} // namespace matrix

