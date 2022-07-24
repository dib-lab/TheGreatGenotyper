add_brwt_library(brwt
  "binary_relation.cpp"
  "bit_vector.cpp"
  "bitmap.cpp"
  "int_vector.cpp"
  "serialization.cpp"
  "wavelet_tree/algorithms.cpp"
  "wavelet_tree/wavelet_tree.cpp"
)

target_include_directories(brwt PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
)

install(TARGETS brwt DESTINATION lib)
