#include <brwt/wavelet_tree.h>
#include <benchmark/benchmark.h>

#include "utility.h"        // gen_integer
#include <brwt/bit_hacks.h> // used_bits
#include <cassert>          // assert
#include <type_traits>      // underlying_type_t, make_unsigned_t
#include <utility>          // pair, make_pair

using brwt::index_type;
using brwt::size_type;
using brwt::symbol_id;
using brwt::wavelet_tree;

using benchmark::DoNotOptimize;
using brwt::benchmark::cyclic_input;
using brwt::benchmark::gen_integer;
using brwt::benchmark::pow_2;

// ==========================================
// Random data generation
// ==========================================

static constexpr symbol_id max_symbol_id(const int bpe) {
  using U = std::underlying_type_t<symbol_id>;
  return symbol_id((U{1} << bpe) - 1);
}

static symbol_id gen_symbol(const symbol_id min, const symbol_id max) {
  return symbol_id(gen_integer<size_type>(min, max));
}

static wavelet_tree gen_wavelet_tree(const size_type count,
                                     const size_type sigma) {
  using U = std::make_unsigned_t<size_type>;
  const auto bpe = brwt::used_bits(static_cast<U>(sigma - 1));
  brwt::int_vector vec(count, bpe);

  const auto min_val = symbol_id(0);
  const auto max_val = max_symbol_id(bpe);

  assert(static_cast<size_type>(max_val + 1) == sigma &&
         "Sigma must be a power of two");

  std::generate(vec.begin(), vec.end(),
                [&] { return gen_symbol(min_val, max_val); });

  return wavelet_tree(vec);
}

static index_type gen_index(const wavelet_tree& wt) {
  assert(wt.size() > 0);
  return gen_integer<index_type>(0, wt.size() - 1);
}

static symbol_id gen_symbol(const wavelet_tree& wt) {
  return gen_symbol(symbol_id(0), wt.max_symbol_id());
}

static auto generate_random_indices(const wavelet_tree& wt,
                                    const std::size_t count) {
  assert(count > 0);

  cyclic_input<index_type> indices;
  indices.generate(count, [&wt] { return gen_index(wt); });
  return indices;
}

static auto generate_random_symbols(const wavelet_tree& wt,
                                    const std::size_t count) {
  assert(count > 0);
  assert(wt.size() > 0);

  cyclic_input<symbol_id> indices;
  indices.generate(count, [&] { return gen_symbol(wt); });
  return indices;
}

static auto generate_select_queries(const wavelet_tree& wt,
                                    const std::size_t count) {
  auto gen_query = [&wt] {
    const auto symbol = gen_symbol(wt);
    const auto total = exclusive_rank(wt, symbol, wt.size());
    const auto nth = (total == 0) ? 1 : gen_integer<size_type>(1, total);
    return std::make_pair(symbol, nth);
  };

  cyclic_input<std::pair<symbol_id, size_type>> queries;
  queries.generate(count, gen_query);
  return queries;
}

// ==========================================
// Benchmark tests
// ==========================================

static void bm_access(benchmark::State& state) {
  const auto wt = gen_wavelet_tree(pow_2(16), state.range(0));
  auto indices = generate_random_indices(wt, 1024);

  while (state.KeepRunning()) {
    const auto idx = indices.next();
    DoNotOptimize(wt.access(idx));
  }
}
BENCHMARK(bm_access)->Range(pow_2(1), pow_2(20));

static void bm_rank(benchmark::State& state) {
  const auto wt = gen_wavelet_tree(pow_2(16), state.range(0));
  auto indices = generate_random_indices(wt, 1024);
  auto symbols = generate_random_symbols(wt, 1019);

  while (state.KeepRunning()) {
    const auto idx = indices.next();
    const auto sym = symbols.next();
    DoNotOptimize(wt.rank(sym, idx));
  }
}
BENCHMARK(bm_rank)->Range(pow_2(1), pow_2(20));

static void bm_select(benchmark::State& state) {
  const auto wt = gen_wavelet_tree(pow_2(16), state.range(0));
  auto queries = generate_select_queries(wt, 1024);
  while (state.KeepRunning()) {
    const auto q = queries.next();
    DoNotOptimize(wt.select(q.first, q.second));
  }
}
BENCHMARK(bm_select)->Range(pow_2(1), pow_2(20));

BENCHMARK_MAIN();
