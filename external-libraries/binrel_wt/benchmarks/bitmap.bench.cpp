#include <brwt/bitmap.h>
#include <benchmark/benchmark.h>

#include "utility.h" // gen_integer
#include <cassert>   // assert
#include <random>    // bernoulli_distribution
#include <utility>   // move

using brwt::bitmap;
using brwt::index_type;
using brwt::size_type;

using brwt::benchmark::cyclic_input;
using brwt::benchmark::gen_integer;
using brwt::benchmark::pow_2;

using benchmark::DoNotOptimize;

// ==========================================
// Random data generation
// ==========================================

static bitmap gen_bitmap(const size_type size, const double density = 0.5) {
  brwt::bit_vector vec(size);
  std::bernoulli_distribution gen_bit(density);

  for (size_type i = 0; i < size; ++i) {
    const bool value = gen_bit(brwt::benchmark::get_random_engine());
    vec.set(i, value);
  }
  return bitmap(std::move(vec));
}

static index_type gen_index(const bitmap& bm) {
  assert(bm.size() > 0);
  return gen_integer<index_type>(0, bm.size() - 1);
}

static auto generate_random_indices(const bitmap& bm,
                                    const std::size_t count = 1024) {
  assert(count > 0);
  cyclic_input<index_type> indices;
  indices.generate(count, [&] { return gen_index(bm); });
  return indices;
}

// ==========================================
// Benchmark tests
// ==========================================

static void bm_access(benchmark::State& state) {
  const auto bm = gen_bitmap(state.range(0));
  auto indices = generate_random_indices(bm);

  while (state.KeepRunning()) {
    const auto idx = indices.next();
    DoNotOptimize(bm.access(idx));
  }
}
BENCHMARK(bm_access)->Range(pow_2(12), pow_2(20));

static void bm_rank_1(benchmark::State& state) {
  const auto bm = gen_bitmap(state.range(0));
  auto indices = generate_random_indices(bm, 1024);

  while (state.KeepRunning()) {
    const auto idx = indices.next();
    DoNotOptimize(bm.rank_1(idx));
  }
}
BENCHMARK(bm_rank_1)->Range(pow_2(12), pow_2(20));

static void bm_select_1(benchmark::State& state) {
  const auto bm = gen_bitmap(state.range(0));
  auto input = cyclic_input<size_type>();
  input.generate(1024,
                 [&] { return gen_integer<size_type>(1, bm.num_ones()); });

  while (state.KeepRunning()) {
    const auto nth = input.next();
    assert(nth >= 1 && nth <= bm.num_ones());

    DoNotOptimize(bm.select_1(nth));
  }
}
BENCHMARK(bm_select_1)->Range(pow_2(12), pow_2(20));

static void bm_select_0(benchmark::State& state) {
  const auto bm = gen_bitmap(state.range(0));
  auto input = cyclic_input<size_type>();
  input.generate(1024,
                 [&] { return gen_integer<size_type>(1, bm.num_zeros()); });

  while (state.KeepRunning()) {
    const auto nth = input.next();
    assert(nth >= 1 && nth <= bm.num_zeros());

    DoNotOptimize(bm.select_0(nth));
  }
}
BENCHMARK(bm_select_0)->Range(pow_2(12), pow_2(20));

BENCHMARK_MAIN();
