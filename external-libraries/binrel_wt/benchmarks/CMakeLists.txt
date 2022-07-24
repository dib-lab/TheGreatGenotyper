find_package(Benchmark) # google-benchmark

if(NOT Benchmark_FOUND)
  message(WARNING "Benchmark library not found. "
                  "Benchmark targets will not be included.")
  return()
endif()

# ---------------------------------------
# Helpers
# ---------------------------------------

function(disable_global_constructor_warning TARGET_ID)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${TARGET_ID} PRIVATE -Wno-global-constructors)
  endif()
endfunction()

# ---------------------------------------
# Targets configuration
# ---------------------------------------

add_custom_target(brwt_benchmarks)

if(BRWT_BUILD_BENCHMARKS)
  set_target_properties(brwt_benchmarks PROPERTIES EXCLUDE_FROM_ALL OFF)
endif()

function(add_benchmark_test BM_NAME)
  set(TARGET_ID "benchmark_${BM_NAME}")
  add_brwt_executable(${TARGET_ID} EXCLUDE_FROM_ALL "${BM_NAME}.bench.cpp")
  disable_global_constructor_warning(${TARGET_ID})

  target_link_libraries(${TARGET_ID} brwt Benchmark::benchmark)
  add_dependencies(brwt_benchmarks ${TARGET_ID})
endfunction()

# ---------------------------------------
# Registration of benchmark tests
# ---------------------------------------

add_benchmark_test("binary_relation")
add_benchmark_test("bit_hacks")
add_benchmark_test("bitmap")
add_benchmark_test("wavelet_tree")
