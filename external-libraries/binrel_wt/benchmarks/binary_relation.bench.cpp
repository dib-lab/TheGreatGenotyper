#include <brwt/binary_relation.h>
#include <benchmark/benchmark.h>

#include "utility.h" // gen_integer, pow_2
#include <cassert>   // assert
#include <cstddef>   // size_t
#include <vector>    // vector

using brwt::binary_relation;
using object_id = binary_relation::object_id;
using label_id = binary_relation::label_id;
using pair_type = binary_relation::pair_type;

using benchmark::DoNotOptimize;
using brwt::benchmark::gen_integer;
using brwt::benchmark::pow_2;

template <typename T>
static constexpr T clamp(const T& v, const T& lo, const T& hi) {
  assert(lo <= hi);
  return v < lo ? lo : hi < v ? hi : v;
}

// ==========================================
// Random data generation
// ==========================================

static object_id gen_object(const object_id min, const object_id max) {
  assert(min <= max);
  return static_cast<object_id>(gen_integer<brwt::size_type>(min, max));
}

static label_id gen_label(const label_id min, const label_id max) {
  assert(min <= max);
  return static_cast<label_id>(gen_integer<brwt::size_type>(min, max));
}

static binary_relation gen_binary_relation(const std::size_t max_size,
                                           const object_id max_object,
                                           const label_id max_label) {

  std::vector<pair_type> pairs;
  pairs.reserve(max_size);

  for (std::size_t i = 0; i < max_size; ++i) {
    const auto obj = gen_object(object_id(0), max_object);
    const auto lab = gen_label(label_id(0), max_label);
    pairs.push_back({obj, lab});
  }

  return binary_relation(pairs);
}

static object_id gen_object(const binary_relation& br) {
  assert(br.object_alphabet_size() > 0);
  const auto max_value = br.object_alphabet_size() - 1;
  return gen_object(object_id(0), object_id(max_value));
}

static label_id gen_label(const binary_relation& br) {
  assert(br.label_alphabet_size() > 0);
  const auto max_value = br.label_alphabet_size() - 1;
  return gen_label(label_id(0), label_id(max_value));
}
static pair_type gen_pair(const binary_relation& br) {
  return pair_type{gen_object(br), gen_label(br)};
}

/// Returns a pair [obj_min, obj_max] valid for \p br.
///
static std::pair<object_id, object_id>
gen_object_range(const binary_relation& br) {
  const auto max = gen_object(br);
  const auto min = gen_object(object_id(0), max);
  assert(min <= max);
  return {min, max};
}

/// Returns a pair [lab_min, lab_max] valid for \p br.
///
static std::pair<label_id, label_id>
gen_label_range(const binary_relation& br) {
  const auto max = gen_label(br);
  const auto min = gen_label(label_id(0), max);
  assert(min <= max);
  return {min, max};
}

// ==========================================
// Benchmark tests
// ==========================================

// All tests are discriminated by the size of the label alphabet as it is the
// most  influential parameter.

static void bm_rank(benchmark::State& state) {
  const auto br = gen_binary_relation(/*max_size=*/1'000'000,
                                      /*max_object=*/object_id(100'000),
                                      /*max_label=*/label_id(state.range(0)));

  while (state.KeepRunning()) {
    auto const max_object = gen_object(br);
    auto const max_label = gen_label(br);

    DoNotOptimize(br.rank(max_object, max_label));
  }
}
BENCHMARK(bm_rank)->Range(pow_2(1), pow_2(20));

static void bm_nth_element_lab_maj(benchmark::State& state) {
  const auto br = gen_binary_relation(/*max_size=*/1'000'000,
                                      /*max_object=*/object_id(100'000),
                                      /*max_label=*/label_id(state.range(0)));

  while (state.KeepRunning()) {
    auto obj_range = gen_object_range(br);
    auto lab_start = gen_label(br);

    DoNotOptimize(br.nth_element(obj_range.first, obj_range.second, lab_start,
                                 42, brwt::lab_major));
  }
}
BENCHMARK(bm_nth_element_lab_maj)->Range(pow_2(1), pow_2(20));

static void bm_nth_element_obj_maj(benchmark::State& state) {
  const auto br = gen_binary_relation(/*max_size=*/1'000'000,
                                      /*max_object=*/object_id(100'000),
                                      /*max_label=*/label_id(state.range(0)));

  while (state.KeepRunning()) {
    auto obj_start = gen_object(br);
    auto lab_range = gen_label_range(br);

    DoNotOptimize(br.nth_element(obj_start, lab_range.first, lab_range.second,
                                 42, brwt::obj_major));
  }
}
BENCHMARK(bm_nth_element_obj_maj)->Range(pow_2(1), pow_2(20));

static void bm_lower_bound(benchmark::State& state) {
  const auto br = gen_binary_relation(/*max_size=*/1'000'000,
                                      /*max_object=*/object_id(100'000),
                                      /*max_label=*/label_id(state.range(0)));

  while (state.KeepRunning()) {
    const auto lab_range = gen_label_range(br);
    auto start = gen_pair(br);
    start.label = clamp(start.label, lab_range.first, lab_range.second);

    DoNotOptimize(br.lower_bound(start, lab_range.first, lab_range.second,
                                 brwt::obj_major));
  }
}
BENCHMARK(bm_lower_bound)->Range(pow_2(1), pow_2(20));

static void bm_obj_select(benchmark::State& state) {
  const auto br = gen_binary_relation(/*max_size=*/1'000'000,
                                      /*max_object=*/object_id(100'000),
                                      /*max_label=*/label_id(state.range(0)));

  while (state.KeepRunning()) {
    auto start = gen_object(br);
    auto label = gen_label(br);

    DoNotOptimize(br.obj_select(start, label, 42));
  }
}
BENCHMARK(bm_obj_select)->Range(pow_2(1), pow_2(20));

static void bm_count_distinct_labels(benchmark::State& state) {
  const auto br = gen_binary_relation(/*max_size=*/1'000'000,
                                      /*max_object=*/object_id(100'000),
                                      /*max_label=*/label_id(state.range(0)));

  while (state.KeepRunning()) {
    const auto obj_range = gen_object_range(br);
    const auto lab_range = gen_label_range(br);

    DoNotOptimize(br.count_distinct_labels(obj_range.first, obj_range.second,
                                           lab_range.first, lab_range.second));
  }
}
BENCHMARK(bm_count_distinct_labels)->Range(pow_2(1), pow_2(20));

BENCHMARK_MAIN();
