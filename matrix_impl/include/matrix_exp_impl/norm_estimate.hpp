namespace matrix { inline namespace v1 {
  namespace exponential {
    template<typename T>
    void populateRandomOnesAndNegOnes(Matrix<T> &mat, size_t c)
    {
#ifndef NDEBUG
      if(c >= mat.ncols()){
        throw IndexOutOfBound("column index out of bound.");
      }
#endif
      std::srand(std::time(nullptr));
      auto data_ptr = mat.data() + c * mat.nrows();
      for(size_t i = 0; i < mat.nrows(); ++i){
        int rand_var = std::rand() % 2;
        if(rand_var == 0){
          *data_ptr = static_cast<T>(-1);
        } else {
          *data_ptr = static_cast<T>(1);
        }
        ++data_ptr;
      }
    }

    template<typename T>
    bool areParallel(const Matrix<T> &A1, size_t c1, const Matrix<T> &A2, size_t c2)
    {
#ifndef NDEBUG
      if(c1 >= A1.ncols() && c2 >= A2.ncols()){
        throw IndexOutOfBound("column indices out of bound.");
      }
      if(A1.nrows() != A2.nrows()){
        throw MatrixSizeMismatchError("two columns need to have the same number of columns");
      }
#endif
      bool same_direction = true;
      bool opposite_direction = true;
      auto ptr_c1 = A1.data() + c1 * A1.nrows();
      auto ptr_c2 = A2.data() + c2 * A2.nrows();
      for(size_t i = 0; i < A1.nrows(); ++i){
        if (!approxEqual(*ptr_c1, *ptr_c2, 2*eps)) same_direction = false; 
        if (!approxEqual(*ptr_c1, -(*ptr_c2, 2*eps))) opposite_direction = false;
        if (!same_direction && !opposite_direction) return false;
      }
      return same_direction || opposite_direction;
    }

    template<typename T>
    double colNorm1(const Matrix<T> &mat, size_t c)
    {
#ifndef NDEBUG
      if(c >= mat.ncols()){
        throw IndexOutOfBound("column indices out of bound.");
      }
#endif
      double res = 0.0;
      const auto pos = mat.data() + c * mat.nrows();
      for(size_t i = 0; i < mat.nrows(); ++i){
        res += std::abs(*pos);
        ++pos;
      }
      return res;
    }

    template<typename T>
    double rowNorm1(const Matrix<T> &mat, size_t r)
    {
#ifndef NDEBUG
      if(r >= mat.nrows()){
        throw IndexOutOfBound("row indices out of bound.");
      }
#endif
      double res = 0.0;
      for(size_t i = 0; i < mat.ncols(); ++i){
        res += std::abs(mat(r,i));
      }
      return res;
    }

    /// Higham, Nicholas J., and Françoise Tisseur. "A block algorithm 
    /// for matrix 1-norm estimation, with an application to 1-norm 
    /// pseudospectra." SIAM Journal on Matrix Analysis and 
    /// Applications 21.4 (2000): 1185-1201.
    /// With generalization to A^m
    template<typename T>
    double normest(const Matrix<T> &A, size_t m)
    {
      constexpr size_t itmax = 2;
      constexpr size_t t = 2;

      auto x = Matrix<T>(A.nrows(), t);
      const double n_norm = static_cast<double>(x.nrows());
      for(size_t i = 0; i < x.nrows(); ++i){
        x(i, 0) = 1/n_norm;
      }
      for(size_t i = 0; i < x.nrows(); ++i){
        x(i, 1) = (i % 2 == 0) ? 1/n_norm : -1/n_norm;
      }

      double est_old = 0.0;
      double est = 0.0;
      auto S = zeros<T>(A.nrows(), t);
      auto S_old = S;
      std::unordered_set<size_t> ind_hist;
      for(size_t k = 1; k <= itmax; ++k){
        auto Y = A * x;   ///< Y is of dim n x 2
        for(size_t i = 1; i < m; ++i){
          Y = A * Y;
        }

        size_t ind_best = 0;
        for(size_t c = 0; c < t; ++c){
          double est_c = colNorm1(Y, c);
          if (est_c > est - eps){
            est = est_c;
            ind_best = c;
          }
        }

        if (k >= 2 && est <= est_old){
          est = est_old;
          return est;
        }

        est_old = est;
        S_old = S;

        if(k > itmax) return est;
        S = sign(Y);

        bool not_all_parallel = false;
        for(size_t i = 0; i < S.ncols(); ++i){
          for(size_t j = 0; j < S_old.ncols(); ++j){
            if(!areParallel(S, i, S_old, j)){
              not_all_parallel = true;
              break;
            }
          }
          if(not_all_parallel) break;
        }
        if(!not_all_parallel) return est;

        if constexpr(t > 1){
          bool has_parallel = true;
          while(has_parallel){
            for(size_t i = 0; i < S.ncols(); ++i){
              size_t j = i + 1;
              while(j < S.ncols()){
                if(areParallel(S, i, S, j)){
                  has_parallel = true;
                  populateRandomOnesAndNegOnes(S, i);
                  break;
                }
                ++j;
              }
              if(j == S.ncols()) has_parallel = false;
              if(has_parallel) break;
              j = 0;
              while(j < S_old.ncols()){
                if(areParallel(S, i, S_old, j)){
                  has_parallel = true;
                  populateRandomOnesAndNegOnes(S, i);
                  break;
                }
                ++j;
              }
              if(j == S_old.ncols()) has_parallel = false;
              if(has_parallel) break;
            }
          }
        }
        
        auto Z = S.adjoint() * A;
        for(size_t i = 1; i < m; ++i){
          Z *= A;
        }
        Z.adjointInplace();
        size_t h_best = 0;
        double max_h = 0.0;
        std::vector<size_t> ind_i;
        for(size_t r = 0; r < Z.ncols(); ++r){
          double est_r = rowNorm1(Z, r);
          if (est_r > max_h - eps){
            max_h = est_r;
            h_best = r;
          }
          ind_i.push_back(r);
        }
        if(k >= 2 && h_best == ind_best) return est;

        std::sort(ind_i.begin(), ind_i.end(), std::less<double>());

        if constexpr(t > 1){
          bool has_extra = false;
          for(size_t i = 0; i < t; ++i){
            if(ind_hist.find(ind_i[i]) == ind_hist.end()) has_extra = true;
          }
          if (!has_extra) return est;
        }

        for(size_t i = 0; i < t; ++i){
          ind_hist.insert(ind_i[i]);
        }

        for(size_t j = 0; j < t; ++j){
          for(size_t i = 0; i < x.nrows(); ++i){
            x(i, j) = (i == j) ? 1 : 0;
          }
        }
      }
      return est;
    }   /// normest()

    template<typename T>
    int ell(const Matrix<T> &A, size_t m)
    {
      const double m_d = static_cast<double>(m);
      const double alpha = std::abs(c(2*m+1)) * normest(abs(A), 2*m_d+1)/norm1(A);
      const double val = std::log2(alpha/eps) * 0.5 / m_d;
      int res = (int)std::max(std::ceil(val), 0.0);
      return res;
    }
  }   // namespace exponential
}  // namespace v1
} // namespace matrix
