target_sources(matrix
  PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/blas_wrapper/cblas_level1.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/blas_wrapper/cblas_level2.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/blas_wrapper/cblas_level3.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/common.cpp"
  )
