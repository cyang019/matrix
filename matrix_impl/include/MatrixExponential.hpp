namespace matrix { inline namespace v1 {
  constexpr double b1 = 64764752532480000;
  constexpr double b2 = 32382376266240000;
  constexpr double b3 = 7771770303897600;
  constexpr double b4 = 1187353796428800;
  constexpr double b5 = 129060195264000;
  constexpr double b6 = 10559470521600;
  constexpr double b7 = 670442572800;
  constexpr double b8 = 33522128640;
  constexpr double b9 = 0;
  constexpr double b10 = 0;
  constexpr double b11 = 0;
  constexpr double b12 = 0;
  constexpr double b13 = 0;

  // squaring and scaling & Pade Approximation
  inline Matrix<cxdbl> exp(const Matrix<cxdbl> &mat)
  {
    Matrix<cxdbl> res(mat.nrows(), mat.ncols());
    return res;
  }

  inline Matrix<double> exp(const Matrix<double> &mat)
  {
    Matrix<double> res(mat.nrows(), mat.ncols());
    return res;
  }

}  // namespace v1
} // namespace matrix

