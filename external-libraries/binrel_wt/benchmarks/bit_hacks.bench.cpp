#include <brwt/bit_hacks.h>
#include <benchmark/benchmark.h>

#include "utility.h" // get_random_engine
#include <algorithm> // generate
#include <array>     // array
#include <random>    // uniform_int_distribution
#include <string>    // to_string

using benchmark::DoNotOptimize;

static void bm_pop_counts(benchmark::State& state) {
  using int_t = unsigned long long;
  using brwt::pop_count;

  std::array<int_t, 8> inputs{};
  std::generate(begin(inputs), end(inputs), [] {
    std::uniform_int_distribution<int_t> dist;
    return dist(brwt::benchmark::get_random_engine());
  });

  while (state.KeepRunning()) {
    for (const auto elem : inputs) {
      DoNotOptimize(pop_count(elem));
    }
  }

  state.SetLabel("pop_count(s) per iteration: " +
                 std::to_string(inputs.size()));
}

BENCHMARK(bm_pop_counts);

BENCHMARK_MAIN();
