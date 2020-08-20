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

            int info = 0;
            mat_zheev(jobz, uplo, n, a.data(), lda,
                eigen_vals.data(), &info);

            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          else if constexpr(em == EigenMethod::zheevd){
            int info = 0;
            mat_zheevd(jobz, uplo, n, a.data(), lda,
                eigen_vals.data(), &info);
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

            int info = 0;
            mat_zheevr(jobz, range, uplo, n, a.data(), lda,
                vl, vu, il, iu, abstol, m,
                eigen_vals.data(), z.get(), ldz, isuppz.get(), 
                &info);
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
            int info = 0;
            mat_zheev(jobz, uplo, mat.nrows(), a.data(), lda,
                eigen_vals.data(), &info);

            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          else if constexpr(em == EigenMethod::zheevd){
            int info = 0;
            mat_zheevd(jobz, uplo, n, a.data(), lda,
                eigen_vals.data(), &info);
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

            int info = 0;
            mat_zheevr(jobz, range, uplo, n, a.data(), lda,
                vl, vu, il, iu, abstol, m,
                eigen_vals.data(), z.get(), ldz, isuppz.get(), 
                &info);
            if(info != 0){
              throw InvalidEigenValue("eigenSys calculation invalid.");
            }
          }
          return eigen_vals;
        }   // eigenVal()

        template<typename T>
        Matrix<T> kroneckerProduct(const Matrix<T> &m1, const Matrix<T> &m2)
        {
          if(m1.nrows() == 0 || m1.ncols() == 0) return m2;
          if(m2.nrows() == 0 || m2.ncols() == 0) return m1;
          const size_t sz1 = m1.nelements();
          const size_t sz2 = m2.nelements();
          const size_t r1 = m1.nrows();
          const size_t c1 = m1.ncols();
          const size_t r2 = m2.nrows();
          const size_t c2 = m2.ncols();
          size_t nrows = r1 * r2;
          size_t ncols = c1 * c2;
          Matrix<T> tmp = zeros<T>(sz2, sz1);
          Matrix<T> res(nrows, ncols);
          if constexpr(is_double<T>::value){
            // outer product
            lvl2_dger(CblasOrder::CblasColMajor,
                sz2, sz1, 1.0,
                m2.data(), 1,
                m1.data(), 1,
                tmp.data(), sz2);
            // rearrange elements
            for(size_t idx = 0; idx < sz1; ++idx){
              for(size_t offset = 0; offset < c2; ++offset){
                const size_t orig_r = offset * r2; 
                const size_t orig_offset = orig_r + idx * sz2;
                const size_t dest_c = (idx / r1) * c2 + offset;
                const size_t dest_r = (idx % r1) * r2;
                const size_t dest_offset = dest_r + dest_c * nrows;
                lvl1_dcopy(
                    r2, 
                    tmp.data() + orig_offset, 1, 
                    res.data() + dest_offset, 1);      
              }
            }
          }
          else if constexpr(is_complex<T>::value){
            // outer product
            lvl2_zgeru(CblasOrder::CblasColMajor,
                sz2, sz1, cxdbl(1.0, 0.0),
                m2.data(), 1,
                m1.data(), 1,
                tmp.data(), sz2);
            // rearrange elements
            for(size_t idx = 0; idx < m1.nelements(); ++idx){
              for(size_t offset = 0; offset < c2; ++offset){
                const size_t orig_r = offset * r2; 
                const size_t orig_offset = orig_r + idx * sz2;
                const size_t dest_c = (idx / r1) * c2 + offset;
                const size_t dest_r = (idx % r1) * r2;
                const size_t dest_offset = dest_r + dest_c * nrows;
                lvl1_zcopy(
                    r2, 
                    tmp.data() + orig_offset, 1, 
                    res.data() + dest_offset, 1);      
              }
            }
          }
          else {    // default other types
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
          }

          return res;
        }

        template<typename T>
        Matrix<T> kroneckerProduct(const std::vector<Matrix<T>> &ms)
        {
          if(ms.size() == 0) return Matrix<T>();

          if(ms.size() == 1) return ms[0];

          auto res = ms[0];
          for(size_t i = 1; i < ms.size(); ++i){
            res = kroneckerProduct(res, ms[i]);
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
          //if(lhs.nrows() != lhs.ncols())
          //  throw MatrixSizeMismatchError("Projection needs to be done on a square matrix.");
#endif
          if constexpr(is_complex<T>::value){
            T res = lvl1_zdotc(lhs.nelements(), lhs.data(), 1, rhs.data(), 1);
            return res;
          } else if constexpr(is_double<T>::value){
            T res = lvl1_ddot(lhs.nelements(), lhs.data(), 1, rhs.data(), 1);
            return res;
          } else {  /// Trace(lhs * rhs.T)
            T res = 0;
            for(size_t i = 0; i < lhs.ncols(); ++i){
              T temp = 0;
              for(size_t j = 0; j < lhs.nrows(); ++j){
                temp += lhs(j,i) * rhs(j,i);
              }
              res += temp;
            }
            return res;
          }
        }   // projection()

        /// @ project lhs onto rhs. rhs will be scaled to norm.
        template<typename T>
        T projectionNorm(const Matrix<T> &lhs, const Matrix<T> &rhs)
        {
          auto val = projection<T>(lhs, rhs);
          auto norm = projection<T>(rhs, rhs);
          if(std::abs(norm) > eps){
            val /= norm;
          }
          return val;
        }   // projectionNorm

        
        template<typename T>
        Matrix<T> pow(const Matrix<T> &mat, std::uint64_t n)
        {
#ifndef NDEBUG
          if(mat.nrows() != mat.ncols()){
            throw MatrixIsNotSquare("Matrix need to be square for pow().");
          }
#endif
          auto result = identity<T>(mat.nrows());
          auto temp = mat;
          while(n > 0){
            if(n & 0x1){
              result *= temp;
            }
            temp *= temp;
            n >>= 1;
          }
          return result;
        }

        template<typename T>
        double norm1(const Matrix<T> &mat)
        {
          double val = 0;
          for(size_t c = 0; c < mat.ncols(); ++c){
            double col_sum = 0;
            for(size_t r = 0; r < mat.nrows(); ++r){
              col_sum += std::abs(mat(r, c));
            }
            if(col_sum > val){
              val = col_sum;
            }
          }
          return val;
        }


        template<typename T>
        Matrix<T> abs(const Matrix<T> &mat)
        {
          auto res = mat;
          for(size_t i = 0; i < mat.ncols(); ++i){
            for(size_t j = 0; j < mat.nrows(); ++j){
              res(j, i) = std::abs(res(j, i));
            }
          }
          return res;
        }

        // solve for A*x = B by minimizing 2-norm(abs(b - A*x))
        std::pair<Matrix<cxdbl>, int> lstsq(const Matrix<cxdbl> &a, const Matrix<cxdbl> &b, double rcond=1.0e-12);

        // solves for A * X = B
        template<typename T>
        Matrix<T> linearSolveSq(Matrix<T> &A, Matrix<T> &B)
        {
#ifndef NDEBUG
          if(A.nrows() != B.nrows()){
            throw MatrixSizeMismatchError("A *  X = B A.nrows() needs to be the same as B.nrows()");
          }
          if(A.nrows() != A.ncols()){
            throw MatrixSizeMismatchError("A needs to be square.");
          }
#endif
          const size_t n = A.nrows();
          const size_t nrhs = B.ncols();
          auto ptr_ipiv = std::make_unique<int[]>(n);
          int info = 0;
          if constexpr(is_complex<T>::value){
            mat_zgesv(n, nrhs, A.data(), n, ptr_ipiv.get(), B.data(), n, &info);
          } else {
            mat_dgesv(n, nrhs, A.data(), n, ptr_ipiv.get(), B.data(), n, &info);
          }
          if (info < 0) {
            std::ostringstream oss("Argument  has an illegal value at: ", std::ios_base::ate);
            oss << -info << ".";
            throw IllegalValue(oss.str());
          } else if (info > 0) {
            std::ostringstream oss("U(", std::ios_base::ate);
            oss << info << "," << info << ") is exactly zero.";
            throw SingularMatrix(oss.str());
          }
          auto res = Matrix<T>(n, nrhs);

          for(size_t i = 0; i < n; ++i){
            const size_t r = ptr_ipiv[i] - 1; ///< starting with 1
#ifndef NDEBUG 
            if (r > B.nrows()) {
              std::cout << "[linearSolveSq()] ptr_ipiv contains invalid numbers:\n";
              std::cout << "[linearSolveSq()] i = " << i << "; r = " << r << "\n";
              std::cout << "[linearSolveSq()] ptr_ipiv[i] = " << ptr_ipiv[i] << "\n";
            }
#endif
            for(size_t j = 0; j < nrhs; ++j){
              res(i, j) = B(r, j);
            }
          }
          return res;
        }

        inline bool isPowerOfTwo(size_t n)
        {
          return n != 0 && (n & (n-1)) == 0;
        }

        // summation of geometric sequence m + m^2 + m^3 + ... + m^n
        template<typename T>
        Matrix<T> geometricSum2Power(const Matrix<T> &sq_mat, size_t n)
        {
#ifndef NDEBUG
          if(sq_mat.nrows() != sq_mat.ncols()) {
            throw MatrixIsNotSquare("need square matrix to calculate geometric sequence.");
          }
          if(!isPowerOfTwo(n)){
            throw IllegalValue("n needs to be power of 2");
          }
#endif
          if(n == 1) {
            return sq_mat;
          }
          if(n == 0) {
            Matrix<T> t1 = zeros<T>(sq_mat.nrows(), sq_mat.ncols());
            return t1;
          }

          Matrix<T> ratio = sq_mat;
          Matrix<T> t1 = ratio;
          size_t i = 1;
          i <<= 1;
          while(i < n) {
            t1 = t1 + ratio * t1;
            ratio = ratio * ratio;
            i <<= 1;
          }
          t1 = t1 + ratio * t1; ///< calculate outside the loop to avoid one extra ratio * ratio calculation
#ifndef NDEBUG
          if(std::isnan(t1(0,0).real()) || std::isnan(t1(0,0).imag())) {
            std::cout << "[WARNING] geometricSum2Power saw NAN." << std::endl;
          }
#endif
          return t1;
        }

        template<typename T>
        Matrix<T> geometricSum(const Matrix<T> &mat, size_t n)
        {
#ifndef NDEBUG
          std::cout << "geometricSum for mat of order " << n << "\n";
#endif
          Matrix<T> result = zeros<T>(mat.nrows(), mat.ncols());
          size_t i = 1;
          while(i < n) {
            if(i & n) {
#ifndef NDEBUG
              std::cout << "\ti: " << i << "\n";
#endif
              if(i == 1) {
                result += mat;
              } else {
                if(n % i == 0) {
                  result += geometricSum2Power(mat, i);
                } else {
                  result += ::matrix::pow(mat, n%i) * geometricSum2Power(mat, i);
                }
              }
            }
            i <<= 1;
          }
          if(i == n) {
#ifndef NDEBUG
              std::cout << "\ti: " << i << "\n";
#endif
            return geometricSum2Power(mat, n);
          }
          return result;
        }

} // namespace v1
} // namespace matrix

