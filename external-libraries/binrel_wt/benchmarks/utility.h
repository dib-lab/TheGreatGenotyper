#ifndef BINREL_WT_BENCHMARKS_UTILITY_H // NOLINT
#define BINREL_WT_BENCHMARKS_UTILITY_H

#include <algorithm>   // generate
#include <cassert>     // assert
#include <cstddef>     // size_t
#include <random>      // mt19937_64, random_device, uniform_int_distribution
#include <type_traits> // is_integral
#include <vector>      // vector

namespace brwt {
namespace benchmark {

/// Returns two nth power of two.
constexpr int pow_2(const int nth) {
  return 1 << nth;
}

inline auto& get_random_engine() {
  static std::mt19937_64 engine{std::random_device{}()};
  return engine;
}

template <typename T>
T gen_integer(const T min, const T max) {
  static_assert(std::is_integral<T>::value, "");
  assert(min <= max);

  std::uniform_int_distribution<T> dist{min, max};
  return dist(get_random_engine());
}

template <typename T>
class cyclic_input {
public:
  template <typename Generator>
  void generate(std::size_t count, Generator g) {
    input.resize(count);
    std::generate(input.begin(), input.end(), g);
    current = input.begin();
  }

  T next() {
    assert(!input.empty());
    if (current == input.end()) {
      current = input.begin();
    }
    return *(current++);
  }

private:
  std::vector<T> input;

  using iter_type = typename decltype(input)::const_iterator;
  iter_type current = input.end();
};

} // end namespace benchmark
} // end namespace brwt

#endif // BINREL_WT_BENCHMARKS_UTILITY_H
