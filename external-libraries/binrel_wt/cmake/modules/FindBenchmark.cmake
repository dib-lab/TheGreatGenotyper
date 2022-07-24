find_path(Benchmark_INCLUDE_DIRS "benchmark/benchmark.h")
find_library(Benchmark_LIBRARIES NAMES "benchmark")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Benchmark
  FOUND_VAR Benchmark_FOUND
  REQUIRED_VARS Benchmark_INCLUDE_DIRS Benchmark_LIBRARIES
)

if(Benchmark_FOUND)
  add_library(Benchmark::benchmark INTERFACE IMPORTED)

  set_target_properties(Benchmark::benchmark PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${Benchmark_INCLUDE_DIRS}
    INTERFACE_LINK_LIBRARIES ${Benchmark_LIBRARIES}
  )
endif()
