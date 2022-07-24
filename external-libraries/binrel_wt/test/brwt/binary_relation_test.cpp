#include "doctest.h"
#include <brwt/binary_relation.h>

#include <algorithm>        // shuffle, for_each
#include <cassert>          // assert
#include <initializer_list> // initializer_list
#include <ostream>          // ostream
#include <random>           // default_random_engine
#include <stdexcept>        // logic_error
#include <type_traits>      // underlying_type_t, is_enum, decay_t
#include <utility>          // pair, move
#include <vector>           // vector

// TODO(Diego): Remove the use of make_test_binary_relation. Use
// make_test_binary_relation_2 instead.

using brwt::binary_relation;
using brwt::index_type;
using brwt::nullopt;
using brwt::size_type;
using object_id = binary_relation::object_id;
using label_id = binary_relation::label_id;
using pair_type = binary_relation::pair_type;

static object_id operator"" _obj(unsigned long long x) {
  return static_cast<object_id>(x);
}

static label_id operator"" _lab(unsigned long long x) {
  return static_cast<label_id>(x);
}

static pair_type pair(const object_id object, const label_id label) {
  return pair_type{object, label};
}

static binary_relation
make_test_binary_relation(const bool remove_labels_from_obj_6 = false) {
  std::vector<pair_type> pairs;
  pairs.reserve(40);

  auto add_pairs = [&](const object_id object,
                       const std::initializer_list<unsigned> labels) {
    std::for_each(begin(labels), end(labels), [&](const auto label) {
      const auto p = pair(object, static_cast<label_id>(label));

      // Add the pair multiple times to ensure that the constructor manages
      // correctly inputs with duplicate entries.
      pairs.push_back(p);
      pairs.push_back(p);
      pairs.push_back(p);
    });
  };

  // In the following matrix, objects are rows and labels are columns.
  // For now, it is required that each object_id has at least one associated
  // pair.
  // Conversely, labels with no associated pair are allowed.

  // The returned binary relation will have the following points.
  //    |0|1|2|3|4|5|6|7|8|9|
  //  0 |_|_|_|_|x|_|_|_|x|_|
  //  1 |_|_|x|_|x|_|_|_|_|_|
  //  2 |_|_|_|_|x|_|_|_|_|_|
  //  3 |_|_|x|_|x|_|x|_|_|_|
  //  4 |_|_|x|_|x|_|_|x|x|_|
  //  5 |_|x|_|x|_|_|_|_|x|x|
  //  6 |x|_|x|_|_|_|_|_|_|_|
  //  7 |_|x|_|x|_|_|x|_|x|_|
  //  8 |_|x|_|x|x|_|_|x|x|x|
  //  9 |x|_|x|_|_|_|x|x|_|_|
  // 10 |_|_|_|x|x|_|_|x|_|x|
  // 11 |_|x|x|_|x|_|_|_|x|_|
  //
  // Intentionally, no pair has label_id = 5
  add_pairs(0_obj, {4, 8});
  add_pairs(1_obj, {2, 4});
  add_pairs(2_obj, {4});
  add_pairs(3_obj, {2, 4, 6});
  add_pairs(4_obj, {2, 4, 7, 8});
  add_pairs(5_obj, {1, 3, 8, 9});
  if (!remove_labels_from_obj_6) {
    add_pairs(6_obj, {0, 2});
  }
  add_pairs(7_obj, {1, 3, 6, 8});
  add_pairs(8_obj, {1, 3, 4, 7, 8, 9});
  add_pairs(9_obj, {0, 2, 6, 7});
  add_pairs(10_obj, {3, 4, 7, 9});
  add_pairs(11_obj, {1, 2, 4, 8});

  assert(pairs.size() == 3 * (remove_labels_from_obj_6 ? 38 : 40));
  assert(std::none_of(begin(pairs), end(pairs),
                      [](const pair_type& p) { return p.label == 5_lab; }));

  std::shuffle(begin(pairs), end(pairs), std::default_random_engine{});
  return binary_relation(pairs);
}

static auto make_test_binary_relation_2() {
  //    |0|1|2|3|4|5|6|7|8|9|
  //  0 |_|_|_|_|x|_|_|_|x|_|
  //  1 |_|_|x|_|x|_|_|_|_|_|
  //  2 |_|_|_|_|x|_|_|_|_|_|
  //  3 |_|_|x|_|x|_|x|_|_|_|
  //  4 |_|_|x|_|x|_|_|x|x|_|
  //  5 |_|x|_|x|_|_|_|_|x|x|
  //  6 |_|_|_|_|_|_|_|_|_|_|
  //  7 |_|x|_|x|_|_|x|_|x|_|
  //  8 |_|x|_|x|x|_|_|x|x|x|
  //  9 |x|_|x|_|_|_|x|x|_|_|
  // 10 |_|_|_|x|x|_|_|x|_|x|
  // 11 |_|x|x|_|x|_|_|_|x|_|
  return make_test_binary_relation(true); // remove labels from obj 6
}

namespace brwt {
template <typename T>
static constexpr auto to_underlying_type(const T value) {
  static_assert(std::is_enum<T>::value, "");
  return static_cast<std::underlying_type_t<T>>(value);
}

static std::ostream& operator<<(std::ostream& os, const object_id value) {
  return os << to_underlying_type(value);
}

static std::ostream& operator<<(std::ostream& os, const label_id value) {
  return os << to_underlying_type(value);
}

static std::ostream& operator<<(std::ostream& os, const pair_type p) {
  return os << '(' << p.object << ", " << p.label << ')';
}

template <typename T>
static std::ostream& operator<<(std::ostream& os, const optional<T>& opt) {
  if (opt) {
    return os << *opt;
  }
  return os << "nullopt";
}

} // end namespace brwt

template <typename Select>
static auto make_select_list(const Select& select) {
  using value_type = std::decay_t<decltype(*select(1))>;
  std::vector<value_type> list;
  for (size_type nth = 1; nth < 1000; ++nth) {
    const auto opt = select(nth);
    if (!opt) {
      return list;
    }
    list.push_back(std::move(*opt));
  }
  throw std::logic_error("Something was wrong with the select function");
}

static auto as_objects(std::initializer_list<int> ilist) {
  std::vector<object_id> vec;
  for (const int elem : ilist) {
    vec.push_back(static_cast<object_id>(elem));
  }
  assert(vec.size() == ilist.size());
  return vec;
}

TEST_CASE("as_objects test") {
  {
    const auto v = as_objects({3, 4, 0, 4, 5});
    REQUIRE(v.size() == 5);
    CHECK(v[0] == 3_obj);
    CHECK(v[1] == 4_obj);
    CHECK(v[2] == 0_obj);
    CHECK(v[3] == 4_obj);
    CHECK(v[4] == 5_obj);
  }
  {
    const auto v = as_objects({42});
    REQUIRE(v.size() == 1);
    CHECK(v[0] == 42_obj);
  }
  {
    const auto v = as_objects({});
    CHECK(v.empty());
  }
}

TEST_SUITE("binary_relation");

TEST_CASE("Vector of pairs constructor") {
  using vec_t = std::vector<pair_type>;
  SUBCASE("Empty vector") {
    const vec_t pairs{};
    const binary_relation br(pairs);
    CHECK(br.size() == 0);
    CHECK(br.object_alphabet_size() == 0);
    CHECK(br.label_alphabet_size() <= 1);
  }
  SUBCASE("Vector with unordered unique entries") {
    const vec_t pairs = {{0_obj, 1_lab}, {1_obj, 2_lab}, {0_obj, 2_lab},
                         {0_obj, 4_lab}, {3_obj, 4_lab}, {3_obj, 2_lab},
                         {2_obj, 1_lab}, {0_obj, 5_lab}, {5_obj, 0_lab}};
    const binary_relation br(pairs);
    CHECK(br.size() == 9);
    CHECK(br.object_alphabet_size() == 6);
    CHECK(br.label_alphabet_size() <= 8);
  }
  SUBCASE("Vector with unordered duplicate entries") {
    const vec_t pairs = {
        {0_obj, 1_lab}, {1_obj, 2_lab}, {0_obj, 2_lab}, {0_obj, 4_lab},
        {3_obj, 4_lab}, {3_obj, 2_lab}, {2_obj, 1_lab}, {0_obj, 5_lab},
        {5_obj, 0_lab}, {0_obj, 4_lab}, {3_obj, 4_lab}, {0_obj, 5_lab},
        {0_obj, 1_lab}, {3_obj, 2_lab}, {0_obj, 4_lab}, {1_obj, 2_lab},
        {3_obj, 2_lab}, {3_obj, 2_lab}, {3_obj, 2_lab}, {0_obj, 4_lab}};
    const binary_relation br(pairs);
    CHECK(br.size() == 9);
    CHECK(br.object_alphabet_size() == 6);
    CHECK(br.label_alphabet_size() <= 8);
  }
  SUBCASE("Vector where the max-label is a power of two") {
    const vec_t pairs = {{0_obj, 4_lab}, {10_obj, 8_lab}, {3_obj, 8_lab}};
    const binary_relation br(pairs);
    CHECK(br.size() == 3);
    CHECK(br.object_alphabet_size() == 11);
    CHECK(br.label_alphabet_size() > 8);
    CHECK(br.label_alphabet_size() <= 16);
  }
  SUBCASE("Vector where the max-label is a power of two minus one") {
    const vec_t pairs = {{54_obj, 5_lab}, {10_obj, 31_lab}, {42_obj, 7_lab}};
    const binary_relation br(pairs);
    CHECK(br.size() == 3);
    CHECK(br.object_alphabet_size() == 55);
    CHECK(br.label_alphabet_size() == 32);
  }
  SUBCASE("Main test vector") {
    const auto br = make_test_binary_relation_2();
    CHECK(br.size() == 38);
    CHECK(br.object_alphabet_size() == 12); // max object is 11
    CHECK(br.label_alphabet_size() <= 16);  // max label is 9
  }
}

TEST_CASE("Size and alphabets size") {
  const auto binrel = make_test_binary_relation_2();
  CHECK(binrel.size() == 38);
  CHECK(binrel.object_alphabet_size() == 12);
  CHECK(binrel.label_alphabet_size() <= 16);
}

TEST_CASE("Rank with max object, max label") {
  const auto binrel = make_test_binary_relation();
  CHECK(binrel.rank(0_obj, 0_lab) == 0);
  CHECK(binrel.rank(0_obj, 9_lab) == 2);
  CHECK(binrel.rank(11_obj, 0_lab) == 2);

  CHECK(binrel.rank(1_obj, 2_lab) == 1);
  CHECK(binrel.rank(2_obj, 3_lab) == 1);
  CHECK(binrel.rank(3_obj, 4_lab) == 6);
  CHECK(binrel.rank(7_obj, 6_lab) == 16);
  CHECK(binrel.rank(8_obj, 2_lab) == 8);
  CHECK(binrel.rank(8_obj, 3_lab) == 11);

  CHECK(binrel.rank(7_obj, 4_lab) == 14);
  CHECK(binrel.rank(7_obj, 5_lab) == 14);
  CHECK(binrel.rank(7_obj, 6_lab) == 16);
  CHECK(binrel.rank(11_obj, 4_lab) == 24);
  CHECK(binrel.rank(11_obj, 5_lab) == 24);
  CHECK(binrel.rank(11_obj, 6_lab) == 27);

  CHECK(binrel.rank(10_obj, 8_lab) == 33);
  CHECK(binrel.rank(10_obj, 9_lab) == 36);
  CHECK(binrel.rank(11_obj, 8_lab) == 37);
  CHECK(binrel.rank(11_obj, 9_lab) == 40);
}

TEST_CASE("Rank with min object, max object, max label") {
  const auto br = make_test_binary_relation();
  CHECK(br.rank(0_obj, 0_obj, 0_lab) == 0);
  CHECK(br.rank(0_obj, 8_obj, 0_lab) == 1);
  CHECK(br.rank(0_obj, 9_obj, 0_lab) == 2);

  CHECK(br.rank(0_obj, 0_obj, 5_lab) == 1);
  CHECK(br.rank(4_obj, 4_obj, 7_lab) == 3);
  CHECK(br.rank(8_obj, 8_obj, 4_lab) == 3);
  CHECK(br.rank(11_obj, 11_obj, 8_lab) == 4);

  CHECK(br.rank(0_obj, 0_obj, 9_lab) == 2);
  CHECK(br.rank(4_obj, 4_obj, 9_lab) == 4);
  CHECK(br.rank(8_obj, 8_obj, 9_lab) == 6);
  CHECK(br.rank(11_obj, 11_obj, 9_lab) == 4);

  CHECK(br.rank(4_obj, 5_obj, 4_lab) == 4);
  CHECK(br.rank(4_obj, 6_obj, 4_lab) == 6);
  CHECK(br.rank(2_obj, 7_obj, 6_lab) == 13);
  CHECK(br.rank(2_obj, 7_obj, 7_lab) == 14);

  CHECK(br.rank(7_obj, 10_obj, 6_lab) == 11);
  CHECK(br.rank(7_obj, 10_obj, 8_lab) == 16);
  CHECK(br.rank(7_obj, 10_obj, 9_lab) == 18);
  CHECK(br.rank(7_obj, 11_obj, 6_lab) == 14);
  CHECK(br.rank(7_obj, 11_obj, 8_lab) == 20);
  CHECK(br.rank(7_obj, 11_obj, 9_lab) == 22);

  CHECK(br.rank(3_obj, 11_obj, 9_lab) == 35);
  CHECK(br.rank(2_obj, 11_obj, 9_lab) == 36);
  CHECK(br.rank(1_obj, 11_obj, 9_lab) == 38);
  CHECK(br.rank(0_obj, 11_obj, 9_lab) == 40);
}

static auto make_label_major_nth_element() {
  return [br = make_test_binary_relation()](const auto... args) {
    using brwt::lab_major;
    return br.nth_element(args..., lab_major);
  };
}

// ==========================================
// n-th element in label major order
// ==========================================

TEST_CASE("[nth_element,lab_major]: Line queries") {
  const auto nth_element = make_label_major_nth_element();

  // first line.
  CHECK(nth_element(0_obj, 0_obj, 0_lab, 1) == pair(0_obj, 4_lab));
  CHECK(nth_element(0_obj, 0_obj, 0_lab, 2) == pair(0_obj, 8_lab));
  CHECK(nth_element(0_obj, 0_obj, 0_lab, 3) == nullopt);
  CHECK(nth_element(0_obj, 0_obj, 4_lab, 1) == pair(0_obj, 4_lab));
  CHECK(nth_element(0_obj, 0_obj, 5_lab, 1) == pair(0_obj, 8_lab));
  CHECK(nth_element(0_obj, 0_obj, 8_lab, 1) == pair(0_obj, 8_lab));

  // random lines
  CHECK(nth_element(2_obj, 2_obj, 2_lab, 1) == pair(2_obj, 4_lab));
  CHECK(nth_element(5_obj, 5_obj, 0_lab, 2) == pair(5_obj, 3_lab));
  CHECK(nth_element(8_obj, 8_obj, 3_lab, 3) == pair(8_obj, 7_lab));
  CHECK(nth_element(10_obj, 10_obj, 6_lab, 2) == pair(10_obj, 9_lab));
}

TEST_CASE("[nth_element,lab_major]: Last pairs and out of range") {
  const auto nth_element = make_label_major_nth_element();

  CHECK(nth_element(2_obj, 6_obj, 0_lab, 14) == pair(5_obj, 9_lab));
  CHECK(nth_element(2_obj, 6_obj, 0_lab, 15) == nullopt);
  CHECK(nth_element(2_obj, 6_obj, 5_lab, 130) == nullopt);
  CHECK(nth_element(2_obj, 6_obj, 5_lab, 5) == pair(5_obj, 9_lab));
  CHECK(nth_element(2_obj, 6_obj, 5_lab, 6) == nullopt);
  CHECK(nth_element(2_obj, 6_obj, 5_lab, 192) == nullopt);

  CHECK(nth_element(4_obj, 9_obj, 0_lab, 24) == pair(8_obj, 9_lab));
  CHECK(nth_element(4_obj, 9_obj, 0_lab, 25) == nullopt);
  CHECK(nth_element(4_obj, 9_obj, 7_lab, 9) == pair(8_obj, 9_lab));
  CHECK(nth_element(4_obj, 9_obj, 7_lab, 10) == nullopt);

  CHECK(nth_element(0_obj, 3_obj, 0_lab, 8) == pair(0_obj, 8_lab));
  CHECK(nth_element(0_obj, 3_obj, 2_lab, 8) == pair(0_obj, 8_lab));
  CHECK(nth_element(0_obj, 3_obj, 2_lab, 9) == nullopt);

  CHECK(nth_element(6_obj, 11_obj, 0_lab, 24) == pair(10_obj, 9_lab));
  CHECK(nth_element(6_obj, 11_obj, 3_lab, 16) == pair(10_obj, 9_lab));
  CHECK(nth_element(6_obj, 11_obj, 3_lab, 17) == nullopt);

  CHECK(nth_element(0_obj, 11_obj, 0_lab, 40) == pair(10_obj, 9_lab));
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 41) == nullopt);
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 2323) == nullopt);
  CHECK(nth_element(0_obj, 11_obj, 4_lab, 24) == pair(10_obj, 9_lab));
  CHECK(nth_element(0_obj, 11_obj, 4_lab, 25) == nullopt);
  CHECK(nth_element(0_obj, 11_obj, 4_lab, 5343) == nullopt);
}

TEST_CASE("[nth_element,lab_major]: Objects ranges") {
  const auto nth_element = make_label_major_nth_element();

  // range of objects [6, 8]
  CHECK(nth_element(6_obj, 8_obj, 0_lab, 1) == pair(6_obj, 0_lab));
  CHECK(nth_element(6_obj, 8_obj, 0_lab, 2) == pair(7_obj, 1_lab));
  CHECK(nth_element(6_obj, 8_obj, 0_lab, 3) == pair(8_obj, 1_lab));
  CHECK(nth_element(6_obj, 8_obj, 0_lab, 7) == pair(8_obj, 4_lab));
  CHECK(nth_element(6_obj, 8_obj, 0_lab, 8) == pair(7_obj, 6_lab));
  CHECK(nth_element(6_obj, 8_obj, 0_lab, 12) == pair(8_obj, 9_lab));
  CHECK(nth_element(6_obj, 8_obj, 4_lab, 1) == pair(8_obj, 4_lab));
  CHECK(nth_element(6_obj, 8_obj, 4_lab, 2) == pair(7_obj, 6_lab));
  CHECK(nth_element(6_obj, 8_obj, 4_lab, 6) == pair(8_obj, 9_lab));

  // Random ranges
  CHECK(nth_element(2_obj, 8_obj, 3_lab, 7) == pair(8_obj, 4_lab));
  CHECK(nth_element(4_obj, 10_obj, 5_lab, 13) == pair(10_obj, 9_lab));
  CHECK(nth_element(5_obj, 11_obj, 3_lab, 15) == pair(8_obj, 8_lab));
  CHECK(nth_element(0_obj, 5_obj, 2_lab, 4) == pair(5_obj, 3_lab));
  CHECK(nth_element(1_obj, 9_obj, 3_lab, 17) == pair(7_obj, 8_lab));
  CHECK(nth_element(0_obj, 11_obj, 4_lab, 21) == pair(11_obj, 8_lab));

  // Full range
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 1) == pair(6_obj, 0_lab));
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 2) == pair(9_obj, 0_lab));
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 3) == pair(5_obj, 1_lab));
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 10) == pair(6_obj, 2_lab));
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 20) == pair(3_obj, 4_lab));
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 30) == pair(9_obj, 7_lab));
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 38) == pair(5_obj, 9_lab));
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 39) == pair(8_obj, 9_lab));
  CHECK(nth_element(0_obj, 11_obj, 0_lab, 40) == pair(10_obj, 9_lab));
}

// ==========================================
// n-th element in object major order
// ==========================================

static auto make_object_major_nth_element() {
  return [br = make_test_binary_relation()](const auto... args) {
    using brwt::obj_major;
    return br.nth_element(args..., obj_major);
  };
}

TEST_CASE("[nth_element,obj_major]: One column test") {
  const auto nth_element = make_object_major_nth_element();

  CHECK(nth_element(0_obj, 0_lab, 0_lab, 1) == pair(6_obj, 0_lab));
  CHECK(nth_element(0_obj, 0_lab, 0_lab, 2) == pair(9_obj, 0_lab));
  CHECK(nth_element(0_obj, 0_lab, 0_lab, 3) == nullopt);
  CHECK(nth_element(0_obj, 2_lab, 2_lab, 1) == pair(1_obj, 2_lab));
  CHECK(nth_element(0_obj, 2_lab, 2_lab, 6) == pair(11_obj, 2_lab));
  CHECK(nth_element(0_obj, 2_lab, 2_lab, 7) == nullopt);
  CHECK(nth_element(0_obj, 5_lab, 5_lab, 1) == nullopt);
  CHECK(nth_element(0_obj, 8_lab, 8_lab, 1) == pair(0_obj, 8_lab));
  CHECK(nth_element(0_obj, 8_lab, 8_lab, 6) == pair(11_obj, 8_lab));
  CHECK(nth_element(0_obj, 8_lab, 8_lab, 7) == nullopt);
  CHECK(nth_element(0_obj, 9_lab, 9_lab, 1) == pair(5_obj, 9_lab));
  CHECK(nth_element(0_obj, 9_lab, 9_lab, 2) == pair(8_obj, 9_lab));
  CHECK(nth_element(0_obj, 9_lab, 9_lab, 3) == pair(10_obj, 9_lab));
  CHECK(nth_element(0_obj, 9_lab, 9_lab, 4) == nullopt);

  CHECK(nth_element(5_obj, 0_lab, 0_lab, 1) == pair(6_obj, 0_lab));
  CHECK(nth_element(6_obj, 0_lab, 0_lab, 1) == pair(6_obj, 0_lab));
  CHECK(nth_element(7_obj, 0_lab, 0_lab, 1) == pair(9_obj, 0_lab));
  CHECK(nth_element(10_obj, 0_lab, 0_lab, 1) == nullopt);
  CHECK(nth_element(5_obj, 5_lab, 5_lab, 1) == nullopt);
  CHECK(nth_element(11_obj, 5_lab, 5_lab, 1) == nullopt);
  CHECK(nth_element(5_obj, 9_lab, 9_lab, 1) == pair(5_obj, 9_lab));
  CHECK(nth_element(6_obj, 9_lab, 9_lab, 1) == pair(8_obj, 9_lab));

  CHECK(nth_element(5_obj, 8_lab, 8_lab, 1) == pair(5_obj, 8_lab));
  CHECK(nth_element(5_obj, 8_lab, 8_lab, 2) == pair(7_obj, 8_lab));
  CHECK(nth_element(5_obj, 8_lab, 8_lab, 3) == pair(8_obj, 8_lab));
  CHECK(nth_element(5_obj, 8_lab, 8_lab, 4) == pair(11_obj, 8_lab));
  CHECK(nth_element(5_obj, 8_lab, 8_lab, 5) == nullopt);
  CHECK(nth_element(5_obj, 8_lab, 8_lab, 42) == nullopt);
}

TEST_CASE("[nth_element,obj_major]: Left range test") {
  const auto nth_element = make_object_major_nth_element();

  CHECK(nth_element(0_obj, 0_lab, 4_lab, 1) == pair(0_obj, 4_lab));
  CHECK(nth_element(0_obj, 0_lab, 4_lab, 2) == pair(1_obj, 2_lab));
  CHECK(nth_element(0_obj, 0_lab, 4_lab, 12) == pair(6_obj, 2_lab));
  CHECK(nth_element(0_obj, 0_lab, 4_lab, 23) == pair(11_obj, 2_lab));
  CHECK(nth_element(0_obj, 0_lab, 4_lab, 24) == pair(11_obj, 4_lab));
  CHECK(nth_element(0_obj, 0_lab, 4_lab, 25) == nullopt);
  CHECK(nth_element(0_obj, 0_lab, 4_lab, 42) == nullopt);

  CHECK(nth_element(6_obj, 0_lab, 4_lab, 1) == pair(6_obj, 0_lab));
  CHECK(nth_element(6_obj, 0_lab, 4_lab, 2) == pair(6_obj, 2_lab));
  CHECK(nth_element(6_obj, 0_lab, 4_lab, 9) == pair(9_obj, 2_lab));
  CHECK(nth_element(6_obj, 0_lab, 4_lab, 13) == pair(11_obj, 2_lab));
  CHECK(nth_element(6_obj, 0_lab, 4_lab, 14) == pair(11_obj, 4_lab));
  CHECK(nth_element(6_obj, 0_lab, 4_lab, 15) == nullopt);
  CHECK(nth_element(6_obj, 0_lab, 4_lab, 42) == nullopt);

  CHECK(nth_element(11_obj, 0_lab, 4_lab, 1) == pair(11_obj, 1_lab));
  CHECK(nth_element(11_obj, 0_lab, 4_lab, 2) == pair(11_obj, 2_lab));
  CHECK(nth_element(11_obj, 0_lab, 4_lab, 3) == pair(11_obj, 4_lab));
  CHECK(nth_element(11_obj, 0_lab, 4_lab, 4) == nullopt);
  CHECK(nth_element(11_obj, 0_lab, 4_lab, 42) == nullopt);
}

TEST_CASE("[nth_element,obj_major]: Right range test") {
  const auto nth_element = make_object_major_nth_element();

  CHECK(nth_element(0_obj, 5_lab, 9_lab, 1) == pair(0_obj, 8_lab));
  CHECK(nth_element(0_obj, 5_lab, 9_lab, 2) == pair(3_obj, 6_lab));
  CHECK(nth_element(0_obj, 5_lab, 9_lab, 10) == pair(8_obj, 8_lab));
  CHECK(nth_element(0_obj, 5_lab, 9_lab, 15) == pair(10_obj, 9_lab));
  CHECK(nth_element(0_obj, 5_lab, 9_lab, 16) == pair(11_obj, 8_lab));
  CHECK(nth_element(0_obj, 5_lab, 9_lab, 17) == nullopt);
  CHECK(nth_element(0_obj, 5_lab, 9_lab, 42) == nullopt);

  CHECK(nth_element(6_obj, 5_lab, 9_lab, 1) == pair(7_obj, 6_lab));
  CHECK(nth_element(6_obj, 5_lab, 9_lab, 2) == pair(7_obj, 8_lab));
  CHECK(nth_element(6_obj, 5_lab, 9_lab, 6) == pair(9_obj, 6_lab));
  CHECK(nth_element(6_obj, 5_lab, 9_lab, 9) == pair(10_obj, 9_lab));
  CHECK(nth_element(6_obj, 5_lab, 9_lab, 10) == pair(11_obj, 8_lab));
  CHECK(nth_element(6_obj, 5_lab, 9_lab, 11) == nullopt);
  CHECK(nth_element(6_obj, 5_lab, 9_lab, 42) == nullopt);

  CHECK(nth_element(11_obj, 5_lab, 9_lab, 1) == pair(11_obj, 8_lab));
  CHECK(nth_element(11_obj, 5_lab, 9_lab, 2) == nullopt);
  CHECK(nth_element(11_obj, 5_lab, 9_lab, 42) == nullopt);
}

TEST_CASE("[nth_element,obj_major]: Center range test") {
  const auto nth_element = make_object_major_nth_element();

  CHECK(nth_element(0_obj, 2_lab, 7_lab, 1) == pair(0_obj, 4_lab));
  CHECK(nth_element(0_obj, 2_lab, 7_lab, 2) == pair(1_obj, 2_lab));
  CHECK(nth_element(0_obj, 2_lab, 7_lab, 12) == pair(6_obj, 2_lab));
  CHECK(nth_element(0_obj, 2_lab, 7_lab, 24) == pair(11_obj, 2_lab));
  CHECK(nth_element(0_obj, 2_lab, 7_lab, 25) == pair(11_obj, 4_lab));
  CHECK(nth_element(0_obj, 2_lab, 7_lab, 26) == nullopt);
  CHECK(nth_element(0_obj, 2_lab, 7_lab, 42) == nullopt);

  CHECK(nth_element(6_obj, 2_lab, 7_lab, 1) == pair(6_obj, 2_lab));
  CHECK(nth_element(6_obj, 2_lab, 7_lab, 2) == pair(7_obj, 3_lab));
  CHECK(nth_element(6_obj, 2_lab, 7_lab, 6) == pair(8_obj, 7_lab));
  CHECK(nth_element(6_obj, 2_lab, 7_lab, 13) == pair(11_obj, 2_lab));
  CHECK(nth_element(6_obj, 2_lab, 7_lab, 14) == pair(11_obj, 4_lab));
  CHECK(nth_element(6_obj, 2_lab, 7_lab, 15) == nullopt);
  CHECK(nth_element(6_obj, 2_lab, 7_lab, 42) == nullopt);

  CHECK(nth_element(11_obj, 2_lab, 7_lab, 1) == pair(11_obj, 2_lab));
  CHECK(nth_element(11_obj, 2_lab, 7_lab, 2) == pair(11_obj, 4_lab));
  CHECK(nth_element(11_obj, 2_lab, 7_lab, 3) == nullopt);
  CHECK(nth_element(11_obj, 2_lab, 7_lab, 42) == nullopt);
}

TEST_CASE("[nth_element,obj_major]: Tiny range test") {
  const auto nth_element = make_object_major_nth_element();

  CHECK(nth_element(0_obj, 5_lab, 7_lab, 1) == pair(3_obj, 6_lab));
  CHECK(nth_element(0_obj, 5_lab, 7_lab, 2) == pair(4_obj, 7_lab));
  CHECK(nth_element(0_obj, 5_lab, 7_lab, 3) == pair(7_obj, 6_lab));
  CHECK(nth_element(0_obj, 5_lab, 7_lab, 4) == pair(8_obj, 7_lab));
  CHECK(nth_element(0_obj, 5_lab, 7_lab, 5) == pair(9_obj, 6_lab));
  CHECK(nth_element(0_obj, 5_lab, 7_lab, 6) == pair(9_obj, 7_lab));
  CHECK(nth_element(0_obj, 5_lab, 7_lab, 7) == pair(10_obj, 7_lab));
  CHECK(nth_element(0_obj, 5_lab, 7_lab, 8) == nullopt);
  CHECK(nth_element(0_obj, 5_lab, 7_lab, 42) == nullopt);

  CHECK(nth_element(6_obj, 5_lab, 7_lab, 1) == pair(7_obj, 6_lab));
  CHECK(nth_element(6_obj, 5_lab, 7_lab, 2) == pair(8_obj, 7_lab));
  CHECK(nth_element(6_obj, 5_lab, 7_lab, 3) == pair(9_obj, 6_lab));
  CHECK(nth_element(6_obj, 5_lab, 7_lab, 4) == pair(9_obj, 7_lab));
  CHECK(nth_element(6_obj, 5_lab, 7_lab, 5) == pair(10_obj, 7_lab));
  CHECK(nth_element(6_obj, 5_lab, 7_lab, 6) == nullopt);
  CHECK(nth_element(6_obj, 5_lab, 7_lab, 42) == nullopt);

  CHECK(nth_element(11_obj, 5_lab, 7_lab, 1) == nullopt);
  CHECK(nth_element(11_obj, 5_lab, 7_lab, 2) == nullopt);
  CHECK(nth_element(11_obj, 5_lab, 7_lab, 3) == nullopt);
  CHECK(nth_element(11_obj, 5_lab, 7_lab, 42) == nullopt);
}

TEST_CASE("[nth_element,obj_major]: Full range test") {
  const auto nth_element = make_object_major_nth_element();

  CHECK(nth_element(0_obj, 0_lab, 9_lab, 1) == pair(0_obj, 4_lab));
  CHECK(nth_element(0_obj, 0_lab, 9_lab, 2) == pair(0_obj, 8_lab));
  CHECK(nth_element(0_obj, 0_lab, 9_lab, 8) == pair(3_obj, 6_lab));
  CHECK(nth_element(0_obj, 0_lab, 9_lab, 24) == pair(8_obj, 3_lab));
  CHECK(nth_element(0_obj, 0_lab, 9_lab, 31) == pair(9_obj, 6_lab));
  CHECK(nth_element(0_obj, 0_lab, 9_lab, 39) == pair(11_obj, 4_lab));
  CHECK(nth_element(0_obj, 0_lab, 9_lab, 40) == pair(11_obj, 8_lab));
  CHECK(nth_element(0_obj, 0_lab, 9_lab, 41) == nullopt);
  CHECK(nth_element(0_obj, 0_lab, 9_lab, 42) == nullopt);
  CHECK(nth_element(0_obj, 0_lab, 9_lab, 3141) == nullopt);

  CHECK(nth_element(6_obj, 0_lab, 9_lab, 1) == pair(6_obj, 0_lab));
  CHECK(nth_element(6_obj, 0_lab, 9_lab, 2) == pair(6_obj, 2_lab));
  CHECK(nth_element(6_obj, 0_lab, 9_lab, 5) == pair(7_obj, 6_lab));
  CHECK(nth_element(6_obj, 0_lab, 9_lab, 10) == pair(8_obj, 7_lab));
  CHECK(nth_element(6_obj, 0_lab, 9_lab, 17) == pair(10_obj, 3_lab));
  CHECK(nth_element(6_obj, 0_lab, 9_lab, 23) == pair(11_obj, 4_lab));
  CHECK(nth_element(6_obj, 0_lab, 9_lab, 24) == pair(11_obj, 8_lab));
  CHECK(nth_element(6_obj, 0_lab, 9_lab, 25) == nullopt);
  CHECK(nth_element(6_obj, 0_lab, 9_lab, 42) == nullopt);
  CHECK(nth_element(6_obj, 0_lab, 9_lab, 3141) == nullopt);

  CHECK(nth_element(11_obj, 0_lab, 9_lab, 1) == pair(11_obj, 1_lab));
  CHECK(nth_element(11_obj, 0_lab, 9_lab, 2) == pair(11_obj, 2_lab));
  CHECK(nth_element(11_obj, 0_lab, 9_lab, 3) == pair(11_obj, 4_lab));
  CHECK(nth_element(11_obj, 0_lab, 9_lab, 4) == pair(11_obj, 8_lab));
  CHECK(nth_element(11_obj, 0_lab, 9_lab, 5) == nullopt);
  CHECK(nth_element(11_obj, 0_lab, 9_lab, 42) == nullopt);
  CHECK(nth_element(11_obj, 0_lab, 9_lab, 3141) == nullopt);
}

// ==========================================
// lower_bound in object major order
// ==========================================

TEST_CASE("[lower_bound,obj_major]: Left, right and empty columns") {
  auto lower_bound = [br = make_test_binary_relation_2()](const pair_type p) {
    using brwt::obj_major;
    return br.lower_bound(p, p.label, p.label, obj_major);
  };

  // left column
  CHECK(lower_bound({0_obj, 0_lab}) == pair(9_obj, 0_lab));
  CHECK(lower_bound({5_obj, 0_lab}) == pair(9_obj, 0_lab));
  CHECK(lower_bound({8_obj, 0_lab}) == pair(9_obj, 0_lab));
  CHECK(lower_bound({9_obj, 0_lab}) == pair(9_obj, 0_lab));
  CHECK(lower_bound({10_obj, 0_lab}) == nullopt);
  CHECK(lower_bound({11_obj, 0_lab}) == nullopt);

  // right column
  CHECK(lower_bound({0_obj, 9_lab}) == pair(5_obj, 9_lab));
  CHECK(lower_bound({4_obj, 9_lab}) == pair(5_obj, 9_lab));
  CHECK(lower_bound({5_obj, 9_lab}) == pair(5_obj, 9_lab));
  CHECK(lower_bound({6_obj, 9_lab}) == pair(8_obj, 9_lab));
  CHECK(lower_bound({7_obj, 9_lab}) == pair(8_obj, 9_lab));
  CHECK(lower_bound({8_obj, 9_lab}) == pair(8_obj, 9_lab));
  CHECK(lower_bound({9_obj, 9_lab}) == pair(10_obj, 9_lab));
  CHECK(lower_bound({10_obj, 9_lab}) == pair(10_obj, 9_lab));
  CHECK(lower_bound({11_obj, 9_lab}) == nullopt);

  // empty column
  CHECK(lower_bound({0_obj, 5_lab}) == nullopt);
  CHECK(lower_bound({5_obj, 5_lab}) == nullopt);
  CHECK(lower_bound({11_obj, 5_lab}) == nullopt);
}

TEST_CASE("[lower_bound,obj_major]: Center columns") {
  auto lower_bound = [br = make_test_binary_relation_2()](const pair_type p) {
    using brwt::obj_major;
    return br.lower_bound(p, p.label, p.label, obj_major);
  };

  // column-id = 3
  CHECK(lower_bound({0_obj, 3_lab}) == pair(5_obj, 3_lab));
  CHECK(lower_bound({6_obj, 3_lab}) == pair(7_obj, 3_lab));
  CHECK(lower_bound({7_obj, 3_lab}) == pair(7_obj, 3_lab));
  CHECK(lower_bound({8_obj, 3_lab}) == pair(8_obj, 3_lab));
  CHECK(lower_bound({9_obj, 3_lab}) == pair(10_obj, 3_lab));
  CHECK(lower_bound({10_obj, 3_lab}) == pair(10_obj, 3_lab));
  CHECK(lower_bound({11_obj, 3_lab}) == nullopt);

  // column-id = 4
  CHECK(lower_bound({0_obj, 4_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({1_obj, 4_lab}) == pair(1_obj, 4_lab));
  CHECK(lower_bound({2_obj, 4_lab}) == pair(2_obj, 4_lab));
  CHECK(lower_bound({3_obj, 4_lab}) == pair(3_obj, 4_lab));
  CHECK(lower_bound({4_obj, 4_lab}) == pair(4_obj, 4_lab));
  CHECK(lower_bound({5_obj, 4_lab}) == pair(8_obj, 4_lab));
  CHECK(lower_bound({6_obj, 4_lab}) == pair(8_obj, 4_lab));
  CHECK(lower_bound({7_obj, 4_lab}) == pair(8_obj, 4_lab));
  CHECK(lower_bound({8_obj, 4_lab}) == pair(8_obj, 4_lab));
  CHECK(lower_bound({9_obj, 4_lab}) == pair(10_obj, 4_lab));
  CHECK(lower_bound({10_obj, 4_lab}) == pair(10_obj, 4_lab));
  CHECK(lower_bound({11_obj, 4_lab}) == pair(11_obj, 4_lab));
}

TEST_CASE("[lower_bound,obj_major]: Left range") {
  auto lower_bound = [br = make_test_binary_relation_2()](const pair_type p) {
    using brwt::obj_major;
    return br.lower_bound(p, 0_lab, 4_lab, obj_major);
  };

  // first row
  CHECK(lower_bound({0_obj, 0_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 1_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 2_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 3_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 4_lab}) == pair(0_obj, 4_lab));

  // middle rows
  CHECK(lower_bound({1_obj, 0_lab}) == pair(1_obj, 2_lab));
  CHECK(lower_bound({1_obj, 1_lab}) == pair(1_obj, 2_lab));
  CHECK(lower_bound({1_obj, 2_lab}) == pair(1_obj, 2_lab));
  CHECK(lower_bound({1_obj, 3_lab}) == pair(1_obj, 4_lab));
  CHECK(lower_bound({1_obj, 4_lab}) == pair(1_obj, 4_lab));

  CHECK(lower_bound({2_obj, 0_lab}) == pair(2_obj, 4_lab));
  CHECK(lower_bound({2_obj, 2_lab}) == pair(2_obj, 4_lab));
  CHECK(lower_bound({2_obj, 4_lab}) == pair(2_obj, 4_lab));

  CHECK(lower_bound({5_obj, 0_lab}) == pair(5_obj, 1_lab));
  CHECK(lower_bound({5_obj, 1_lab}) == pair(5_obj, 1_lab));
  CHECK(lower_bound({5_obj, 2_lab}) == pair(5_obj, 3_lab));
  CHECK(lower_bound({5_obj, 3_lab}) == pair(5_obj, 3_lab));
  CHECK(lower_bound({5_obj, 4_lab}) == pair(7_obj, 1_lab));

  CHECK(lower_bound({9_obj, 0_lab}) == pair(9_obj, 0_lab));
  CHECK(lower_bound({9_obj, 2_lab}) == pair(9_obj, 2_lab));
  CHECK(lower_bound({9_obj, 4_lab}) == pair(10_obj, 3_lab));

  // empty row
  CHECK(lower_bound({6_obj, 0_lab}) == pair(7_obj, 1_lab));
  CHECK(lower_bound({6_obj, 2_lab}) == pair(7_obj, 1_lab));
  CHECK(lower_bound({6_obj, 4_lab}) == pair(7_obj, 1_lab));

  // last row
  CHECK(lower_bound({11_obj, 0_lab}) == pair(11_obj, 1_lab));
  CHECK(lower_bound({11_obj, 1_lab}) == pair(11_obj, 1_lab));
  CHECK(lower_bound({11_obj, 2_lab}) == pair(11_obj, 2_lab));
  CHECK(lower_bound({11_obj, 3_lab}) == pair(11_obj, 4_lab));
  CHECK(lower_bound({11_obj, 4_lab}) == pair(11_obj, 4_lab));
}

TEST_CASE("[lower_bound,obj_major]: Right range") {
  auto lower_bound = [br = make_test_binary_relation_2()](const pair_type p) {
    using brwt::obj_major;
    return br.lower_bound(p, 5_lab, 9_lab, obj_major);
  };

  // first row
  CHECK(lower_bound({0_obj, 5_lab}) == pair(0_obj, 8_lab));
  CHECK(lower_bound({0_obj, 6_lab}) == pair(0_obj, 8_lab));
  CHECK(lower_bound({0_obj, 7_lab}) == pair(0_obj, 8_lab));
  CHECK(lower_bound({0_obj, 8_lab}) == pair(0_obj, 8_lab));
  CHECK(lower_bound({0_obj, 9_lab}) == pair(3_obj, 6_lab));

  // middle rows
  CHECK(lower_bound({3_obj, 5_lab}) == pair(3_obj, 6_lab));
  CHECK(lower_bound({3_obj, 6_lab}) == pair(3_obj, 6_lab));
  CHECK(lower_bound({3_obj, 7_lab}) == pair(4_obj, 7_lab));
  CHECK(lower_bound({3_obj, 8_lab}) == pair(4_obj, 7_lab));
  CHECK(lower_bound({3_obj, 9_lab}) == pair(4_obj, 7_lab));

  CHECK(lower_bound({5_obj, 5_lab}) == pair(5_obj, 8_lab));
  CHECK(lower_bound({5_obj, 7_lab}) == pair(5_obj, 8_lab));
  CHECK(lower_bound({5_obj, 9_lab}) == pair(5_obj, 9_lab));

  CHECK(lower_bound({7_obj, 5_lab}) == pair(7_obj, 6_lab));
  CHECK(lower_bound({7_obj, 6_lab}) == pair(7_obj, 6_lab));
  CHECK(lower_bound({7_obj, 7_lab}) == pair(7_obj, 8_lab));
  CHECK(lower_bound({7_obj, 8_lab}) == pair(7_obj, 8_lab));
  CHECK(lower_bound({7_obj, 9_lab}) == pair(8_obj, 7_lab));

  CHECK(lower_bound({9_obj, 5_lab}) == pair(9_obj, 6_lab));
  CHECK(lower_bound({9_obj, 7_lab}) == pair(9_obj, 7_lab));
  CHECK(lower_bound({9_obj, 9_lab}) == pair(10_obj, 7_lab));

  // empty rows
  CHECK(lower_bound({1_obj, 5_lab}) == pair(3_obj, 6_lab));
  CHECK(lower_bound({1_obj, 7_lab}) == pair(3_obj, 6_lab));
  CHECK(lower_bound({1_obj, 9_lab}) == pair(3_obj, 6_lab));
  CHECK(lower_bound({2_obj, 6_lab}) == pair(3_obj, 6_lab));
  CHECK(lower_bound({2_obj, 8_lab}) == pair(3_obj, 6_lab));

  CHECK(lower_bound({6_obj, 5_lab}) == pair(7_obj, 6_lab));
  CHECK(lower_bound({6_obj, 8_lab}) == pair(7_obj, 6_lab));
  CHECK(lower_bound({6_obj, 9_lab}) == pair(7_obj, 6_lab));

  // last row
  CHECK(lower_bound({11_obj, 5_lab}) == pair(11_obj, 8_lab));
  CHECK(lower_bound({11_obj, 6_lab}) == pair(11_obj, 8_lab));
  CHECK(lower_bound({11_obj, 7_lab}) == pair(11_obj, 8_lab));
  CHECK(lower_bound({11_obj, 8_lab}) == pair(11_obj, 8_lab));
  CHECK(lower_bound({11_obj, 9_lab}) == nullopt);
}

TEST_CASE("[lower_bound,obj_major]: Center range") {
  auto lower_bound = [br = make_test_binary_relation_2()](const pair_type p) {
    using brwt::obj_major;
    return br.lower_bound(p, 2_lab, 7_lab, obj_major);
  };

  // first row
  CHECK(lower_bound({0_obj, 2_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 3_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 4_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 5_lab}) == pair(1_obj, 2_lab));
  CHECK(lower_bound({0_obj, 6_lab}) == pair(1_obj, 2_lab));
  CHECK(lower_bound({0_obj, 7_lab}) == pair(1_obj, 2_lab));

  // middle rows
  CHECK(lower_bound({2_obj, 2_lab}) == pair(2_obj, 4_lab));
  CHECK(lower_bound({2_obj, 3_lab}) == pair(2_obj, 4_lab));
  CHECK(lower_bound({2_obj, 4_lab}) == pair(2_obj, 4_lab));
  CHECK(lower_bound({2_obj, 5_lab}) == pair(3_obj, 2_lab));
  CHECK(lower_bound({2_obj, 6_lab}) == pair(3_obj, 2_lab));
  CHECK(lower_bound({2_obj, 7_lab}) == pair(3_obj, 2_lab));

  CHECK(lower_bound({5_obj, 2_lab}) == pair(5_obj, 3_lab));
  CHECK(lower_bound({5_obj, 3_lab}) == pair(5_obj, 3_lab));
  CHECK(lower_bound({5_obj, 4_lab}) == pair(7_obj, 3_lab));

  CHECK(lower_bound({7_obj, 2_lab}) == pair(7_obj, 3_lab));
  CHECK(lower_bound({7_obj, 3_lab}) == pair(7_obj, 3_lab));
  CHECK(lower_bound({7_obj, 4_lab}) == pair(7_obj, 6_lab));

  CHECK(lower_bound({9_obj, 5_lab}) == pair(9_obj, 6_lab));
  CHECK(lower_bound({9_obj, 6_lab}) == pair(9_obj, 6_lab));
  CHECK(lower_bound({9_obj, 7_lab}) == pair(9_obj, 7_lab));

  // empty row
  CHECK(lower_bound({6_obj, 2_lab}) == pair(7_obj, 3_lab));
  CHECK(lower_bound({6_obj, 5_lab}) == pair(7_obj, 3_lab));
  CHECK(lower_bound({6_obj, 7_lab}) == pair(7_obj, 3_lab));

  // last row
  CHECK(lower_bound({11_obj, 2_lab}) == pair(11_obj, 2_lab));
  CHECK(lower_bound({11_obj, 3_lab}) == pair(11_obj, 4_lab));
  CHECK(lower_bound({11_obj, 4_lab}) == pair(11_obj, 4_lab));
  CHECK(lower_bound({11_obj, 5_lab}) == nullopt);
  CHECK(lower_bound({11_obj, 6_lab}) == nullopt);
  CHECK(lower_bound({11_obj, 7_lab}) == nullopt);
}

TEST_CASE("[lower_bound,obj_major]: Full range") {
  auto lower_bound = [br = make_test_binary_relation_2()](const pair_type p) {
    using brwt::obj_major;
    return br.lower_bound(p, 0_lab, 9_lab, obj_major);
  };

  // first row
  CHECK(lower_bound({0_obj, 0_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 3_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 4_lab}) == pair(0_obj, 4_lab));
  CHECK(lower_bound({0_obj, 5_lab}) == pair(0_obj, 8_lab));

  // middle rows
  CHECK(lower_bound({1_obj, 1_lab}) == pair(1_obj, 2_lab));
  CHECK(lower_bound({1_obj, 2_lab}) == pair(1_obj, 2_lab));
  CHECK(lower_bound({1_obj, 3_lab}) == pair(1_obj, 4_lab));
  CHECK(lower_bound({1_obj, 4_lab}) == pair(1_obj, 4_lab));
  CHECK(lower_bound({1_obj, 5_lab}) == pair(2_obj, 4_lab));
  CHECK(lower_bound({1_obj, 9_lab}) == pair(2_obj, 4_lab));

  CHECK(lower_bound({2_obj, 4_lab}) == pair(2_obj, 4_lab));
  CHECK(lower_bound({2_obj, 5_lab}) == pair(3_obj, 2_lab));
  CHECK(lower_bound({4_obj, 3_lab}) == pair(4_obj, 4_lab));
  CHECK(lower_bound({5_obj, 3_lab}) == pair(5_obj, 3_lab));
  CHECK(lower_bound({7_obj, 4_lab}) == pair(7_obj, 6_lab));

  // empty row
  CHECK(lower_bound({6_obj, 0_lab}) == pair(7_obj, 1_lab));
  CHECK(lower_bound({6_obj, 5_lab}) == pair(7_obj, 1_lab));
  CHECK(lower_bound({6_obj, 9_lab}) == pair(7_obj, 1_lab));

  // last row
  CHECK(lower_bound({11_obj, 0_lab}) == pair(11_obj, 1_lab));
  CHECK(lower_bound({11_obj, 1_lab}) == pair(11_obj, 1_lab));
  CHECK(lower_bound({11_obj, 2_lab}) == pair(11_obj, 2_lab));
  CHECK(lower_bound({11_obj, 3_lab}) == pair(11_obj, 4_lab));
  CHECK(lower_bound({11_obj, 7_lab}) == pair(11_obj, 8_lab));
  CHECK(lower_bound({11_obj, 8_lab}) == pair(11_obj, 8_lab));
  CHECK(lower_bound({11_obj, 9_lab}) == nullopt);
}

TEST_CASE("[lower_bound,obj_major]: Thin range") {
  auto lower_bound = [br = make_test_binary_relation_2()](const pair_type p) {
    using brwt::obj_major;
    return br.lower_bound(p, 5_lab, 6_lab, obj_major);
  };

  CHECK(lower_bound({0_obj, 5_lab}) == pair(3_obj, 6_lab));
  CHECK(lower_bound({3_obj, 5_lab}) == pair(3_obj, 6_lab));
  CHECK(lower_bound({3_obj, 6_lab}) == pair(3_obj, 6_lab));

  CHECK(lower_bound({4_obj, 5_lab}) == pair(7_obj, 6_lab));
  CHECK(lower_bound({6_obj, 6_lab}) == pair(7_obj, 6_lab));
  CHECK(lower_bound({7_obj, 5_lab}) == pair(7_obj, 6_lab));
  CHECK(lower_bound({7_obj, 6_lab}) == pair(7_obj, 6_lab));

  CHECK(lower_bound({8_obj, 5_lab}) == pair(9_obj, 6_lab));
  CHECK(lower_bound({8_obj, 6_lab}) == pair(9_obj, 6_lab));
  CHECK(lower_bound({9_obj, 5_lab}) == pair(9_obj, 6_lab));
  CHECK(lower_bound({9_obj, 6_lab}) == pair(9_obj, 6_lab));

  CHECK(lower_bound({10_obj, 5_lab}) == nullopt);
  CHECK(lower_bound({10_obj, 6_lab}) == nullopt);

  CHECK(lower_bound({11_obj, 5_lab}) == nullopt);
  CHECK(lower_bound({11_obj, 6_lab}) == nullopt);
}

// ==========================================
// obj[_exclusive]_rank with fixed_label
// ==========================================

TEST_CASE("[obj_rank, fixed_label]: Left, center, empty and right") {
  using p = std::pair<size_type, size_type>;
  auto rank_pair = [br = make_test_binary_relation_2()](auto... args) {
    return p(br.obj_exclusive_rank(args...), br.obj_rank(args...));
  };

  CHECK(rank_pair(0_obj, 0_lab) == p(0, 0));
  CHECK(rank_pair(8_obj, 0_lab) == p(0, 0));
  CHECK(rank_pair(9_obj, 0_lab) == p(0, 1));
  CHECK(rank_pair(10_obj, 0_lab) == p(1, 1));
  CHECK(rank_pair(11_obj, 0_lab) == p(1, 1));

  CHECK(rank_pair(0_obj, 3_lab) == p(0, 0));
  CHECK(rank_pair(4_obj, 3_lab) == p(0, 0));
  CHECK(rank_pair(5_obj, 3_lab) == p(0, 1));
  CHECK(rank_pair(6_obj, 3_lab) == p(1, 1));
  CHECK(rank_pair(7_obj, 3_lab) == p(1, 2));
  CHECK(rank_pair(8_obj, 3_lab) == p(2, 3));
  CHECK(rank_pair(9_obj, 3_lab) == p(3, 3));
  CHECK(rank_pair(10_obj, 3_lab) == p(3, 4));
  CHECK(rank_pair(11_obj, 3_lab) == p(4, 4));

  CHECK(rank_pair(0_obj, 5_lab) == p(0, 0));
  CHECK(rank_pair(5_obj, 5_lab) == p(0, 0));
  CHECK(rank_pair(11_obj, 5_lab) == p(0, 0));

  CHECK(rank_pair(0_obj, 9_lab) == p(0, 0));
  CHECK(rank_pair(4_obj, 9_lab) == p(0, 0));
  CHECK(rank_pair(5_obj, 9_lab) == p(0, 1));
  CHECK(rank_pair(6_obj, 9_lab) == p(1, 1));
  CHECK(rank_pair(7_obj, 9_lab) == p(1, 1));
  CHECK(rank_pair(8_obj, 9_lab) == p(1, 2));
  CHECK(rank_pair(9_obj, 9_lab) == p(2, 2));
  CHECK(rank_pair(10_obj, 9_lab) == p(2, 3));
  CHECK(rank_pair(11_obj, 9_lab) == p(3, 3));
}

TEST_CASE("[obj_rank, fixed_label]: First and last object associated") {
  using p = std::pair<size_type, size_type>;
  auto rank_pair = [br = make_test_binary_relation_2()](auto... args) {
    return p(br.obj_exclusive_rank(args...), br.obj_rank(args...));
  };

  CHECK(rank_pair(0_obj, 4_lab) == p(0, 1));
  CHECK(rank_pair(1_obj, 4_lab) == p(1, 2));
  CHECK(rank_pair(2_obj, 4_lab) == p(2, 3));
  CHECK(rank_pair(3_obj, 4_lab) == p(3, 4));
  CHECK(rank_pair(4_obj, 4_lab) == p(4, 5));
  CHECK(rank_pair(5_obj, 4_lab) == p(5, 5));
  CHECK(rank_pair(6_obj, 4_lab) == p(5, 5));
  CHECK(rank_pair(7_obj, 4_lab) == p(5, 5));
  CHECK(rank_pair(8_obj, 4_lab) == p(5, 6));
  CHECK(rank_pair(9_obj, 4_lab) == p(6, 6));
  CHECK(rank_pair(10_obj, 4_lab) == p(6, 7));
  CHECK(rank_pair(11_obj, 4_lab) == p(7, 8));
}

// ==========================================
// obj[_exclusive]_rank with label range
// ==========================================

TEST_CASE("[obj_rank,label_range]: Single label range") {
  using p = std::pair<size_type, size_type>;
  auto rank_pair = [br = make_test_binary_relation_2()](auto... args) {
    return p(br.obj_exclusive_rank(args...), br.obj_rank(args...));
  };

  CHECK(rank_pair(0_obj, 0_lab, 0_lab) == p(0, 0));
  CHECK(rank_pair(8_obj, 0_lab, 0_lab) == p(0, 0));
  CHECK(rank_pair(9_obj, 0_lab, 0_lab) == p(0, 1));
  CHECK(rank_pair(10_obj, 0_lab, 0_lab) == p(1, 1));
  CHECK(rank_pair(11_obj, 0_lab, 0_lab) == p(1, 1));

  CHECK(rank_pair(0_obj, 4_lab, 4_lab) == p(0, 1));
  CHECK(rank_pair(4_obj, 4_lab, 4_lab) == p(4, 5));
  CHECK(rank_pair(5_obj, 4_lab, 4_lab) == p(5, 5));
  CHECK(rank_pair(9_obj, 4_lab, 4_lab) == p(6, 6));
  CHECK(rank_pair(10_obj, 4_lab, 4_lab) == p(6, 7));
  CHECK(rank_pair(11_obj, 4_lab, 4_lab) == p(7, 8));

  CHECK(rank_pair(0_obj, 5_lab, 5_lab) == p(0, 0));
  CHECK(rank_pair(5_obj, 5_lab, 5_lab) == p(0, 0));
  CHECK(rank_pair(11_obj, 5_lab, 5_lab) == p(0, 0));

  CHECK(rank_pair(0_obj, 9_lab, 9_lab) == p(0, 0));
  CHECK(rank_pair(7_obj, 9_lab, 9_lab) == p(1, 1));
  CHECK(rank_pair(8_obj, 9_lab, 9_lab) == p(1, 2));
  CHECK(rank_pair(9_obj, 9_lab, 9_lab) == p(2, 2));
  CHECK(rank_pair(10_obj, 9_lab, 9_lab) == p(2, 3));
  CHECK(rank_pair(11_obj, 9_lab, 9_lab) == p(3, 3));
}

TEST_CASE("[obj_rank,label_range]: Left range") {
  using p = std::pair<size_type, size_type>;
  auto rank_pair = [br = make_test_binary_relation_2()](auto... args) {
    return p(br.obj_exclusive_rank(args...), br.obj_rank(args...));
  };

  CHECK(rank_pair(0_obj, 0_lab, 3_lab) == p(0, 0));
  CHECK(rank_pair(1_obj, 0_lab, 3_lab) == p(0, 1));
  CHECK(rank_pair(2_obj, 0_lab, 3_lab) == p(1, 1));
  CHECK(rank_pair(3_obj, 0_lab, 3_lab) == p(1, 2));
  CHECK(rank_pair(4_obj, 0_lab, 3_lab) == p(2, 3));
  CHECK(rank_pair(5_obj, 0_lab, 3_lab) == p(3, 5));
  CHECK(rank_pair(6_obj, 0_lab, 3_lab) == p(5, 5));
  CHECK(rank_pair(7_obj, 0_lab, 3_lab) == p(5, 7));
  CHECK(rank_pair(8_obj, 0_lab, 3_lab) == p(7, 9));
  CHECK(rank_pair(9_obj, 0_lab, 3_lab) == p(9, 11));
  CHECK(rank_pair(10_obj, 0_lab, 3_lab) == p(11, 12));
  CHECK(rank_pair(11_obj, 0_lab, 3_lab) == p(12, 14));
}

TEST_CASE("[obj_rank,label_range]: Right range") {
  using p = std::pair<size_type, size_type>;
  auto rank_pair = [br = make_test_binary_relation_2()](auto... args) {
    return p(br.obj_exclusive_rank(args...), br.obj_rank(args...));
  };

  CHECK(rank_pair(0_obj, 6_lab, 9_lab) == p(0, 1));
  CHECK(rank_pair(1_obj, 6_lab, 9_lab) == p(1, 1));
  CHECK(rank_pair(2_obj, 6_lab, 9_lab) == p(1, 1));
  CHECK(rank_pair(3_obj, 6_lab, 9_lab) == p(1, 2));
  CHECK(rank_pair(4_obj, 6_lab, 9_lab) == p(2, 4));
  CHECK(rank_pair(5_obj, 6_lab, 9_lab) == p(4, 6));
  CHECK(rank_pair(6_obj, 6_lab, 9_lab) == p(6, 6));
  CHECK(rank_pair(7_obj, 6_lab, 9_lab) == p(6, 8));
  CHECK(rank_pair(8_obj, 6_lab, 9_lab) == p(8, 11));
  CHECK(rank_pair(9_obj, 6_lab, 9_lab) == p(11, 13));
  CHECK(rank_pair(10_obj, 6_lab, 9_lab) == p(13, 15));
  CHECK(rank_pair(11_obj, 6_lab, 9_lab) == p(15, 16));
}

TEST_CASE("[obj_rank,label_range]: Center range") {
  using p = std::pair<size_type, size_type>;
  auto rank_pair = [br = make_test_binary_relation_2()](auto... args) {
    return p(br.obj_exclusive_rank(args...), br.obj_rank(args...));
  };

  CHECK(rank_pair(0_obj, 2_lab, 7_lab) == p(0, 1));
  CHECK(rank_pair(1_obj, 2_lab, 7_lab) == p(1, 3));
  CHECK(rank_pair(2_obj, 2_lab, 7_lab) == p(3, 4));
  CHECK(rank_pair(3_obj, 2_lab, 7_lab) == p(4, 7));
  CHECK(rank_pair(4_obj, 2_lab, 7_lab) == p(7, 10));
  CHECK(rank_pair(5_obj, 2_lab, 7_lab) == p(10, 11));
  CHECK(rank_pair(6_obj, 2_lab, 7_lab) == p(11, 11));
  CHECK(rank_pair(7_obj, 2_lab, 7_lab) == p(11, 13));
  CHECK(rank_pair(8_obj, 2_lab, 7_lab) == p(13, 16));
  CHECK(rank_pair(9_obj, 2_lab, 7_lab) == p(16, 19));
  CHECK(rank_pair(10_obj, 2_lab, 7_lab) == p(19, 22));
  CHECK(rank_pair(11_obj, 2_lab, 7_lab) == p(22, 24));
}

TEST_CASE("[obj_rank,label_range]: Full range") {
  using p = std::pair<size_type, size_type>;
  auto rank_pair = [br = make_test_binary_relation_2()](auto... args) {
    return p(br.obj_exclusive_rank(args...), br.obj_rank(args...));
  };

  CHECK(rank_pair(0_obj, 0_lab, 9_lab) == p(0, 2));
  CHECK(rank_pair(1_obj, 0_lab, 9_lab) == p(2, 4));
  CHECK(rank_pair(2_obj, 0_lab, 9_lab) == p(4, 5));
  CHECK(rank_pair(3_obj, 0_lab, 9_lab) == p(5, 8));
  CHECK(rank_pair(4_obj, 0_lab, 9_lab) == p(8, 12));
  CHECK(rank_pair(5_obj, 0_lab, 9_lab) == p(12, 16));
  CHECK(rank_pair(6_obj, 0_lab, 9_lab) == p(16, 16));
  CHECK(rank_pair(7_obj, 0_lab, 9_lab) == p(16, 20));
  CHECK(rank_pair(8_obj, 0_lab, 9_lab) == p(20, 26));
  CHECK(rank_pair(9_obj, 0_lab, 9_lab) == p(26, 30));
  CHECK(rank_pair(10_obj, 0_lab, 9_lab) == p(30, 34));
  CHECK(rank_pair(11_obj, 0_lab, 9_lab) == p(34, 38));
}

// ==========================================
// obj_select with fixed label
// ==========================================

TEST_CASE("[obj_select,fixed_label]") {
  auto select_list = [br = make_test_binary_relation_2()](
                         const object_id start, const label_id label) {
    auto select = [&](const size_type nth) {
      return br.obj_select(start, label, nth);
    };
    return make_select_list(select);
  };

  // Left column
  CHECK(select_list(0_obj, 0_lab) == as_objects({9}));
  CHECK(select_list(5_obj, 0_lab) == as_objects({9}));
  CHECK(select_list(9_obj, 0_lab) == as_objects({9}));
  CHECK(select_list(10_obj, 0_lab).empty());
  CHECK(select_list(11_obj, 0_lab).empty());

  // Right column.
  CHECK(select_list(0_obj, 9_lab) == as_objects({5, 8, 10}));
  CHECK(select_list(5_obj, 9_lab) == as_objects({5, 8, 10}));
  CHECK(select_list(8_obj, 9_lab) == as_objects({8, 10}));
  CHECK(select_list(9_obj, 9_lab) == as_objects({10}));
  CHECK(select_list(10_obj, 9_lab) == as_objects({10}));
  CHECK(select_list(11_obj, 9_lab).empty());

  // Center column with few objects.
  CHECK(select_list(0_obj, 3_lab) == as_objects({5, 7, 8, 10}));
  CHECK(select_list(5_obj, 3_lab) == as_objects({5, 7, 8, 10}));
  CHECK(select_list(7_obj, 3_lab) == as_objects({7, 8, 10}));
  CHECK(select_list(8_obj, 3_lab) == as_objects({8, 10}));
  CHECK(select_list(10_obj, 3_lab) == as_objects({10}));
  CHECK(select_list(11_obj, 3_lab).empty());

  // Center column with many objects.
  CHECK(select_list(0_obj, 4_lab) == as_objects({0, 1, 2, 3, 4, 8, 10, 11}));
  CHECK(select_list(3_obj, 4_lab) == as_objects({3, 4, 8, 10, 11}));
  CHECK(select_list(4_obj, 4_lab) == as_objects({4, 8, 10, 11}));
  CHECK(select_list(9_obj, 4_lab) == as_objects({10, 11}));
  CHECK(select_list(10_obj, 4_lab) == as_objects({10, 11}));
  CHECK(select_list(11_obj, 4_lab) == as_objects({11}));

  // Empty column.
  CHECK(select_list(0_obj, 5_lab).empty());
  CHECK(select_list(5_obj, 5_lab).empty());
  CHECK(select_list(11_obj, 5_lab).empty());
}

// ==========================================
// Count distinct labels tests
// ==========================================

TEST_CASE("count_distinct_labels, basic") {
  const auto br = make_test_binary_relation();
  // Rectangle starting from (0, 0)
  CHECK(br.count_distinct_labels(0_obj, 0_obj, 0_lab, 0_lab) == 0);
  CHECK(br.count_distinct_labels(0_obj, 4_obj, 0_lab, 4_lab) == 2);
  CHECK(br.count_distinct_labels(0_obj, 8_obj, 0_lab, 8_lab) == 8);
  CHECK(br.count_distinct_labels(0_obj, 11_obj, 0_lab, 9_lab) == 9);

  // Rows
  CHECK(br.count_distinct_labels(0_obj, 0_obj, 0_lab, 9_lab) == 2);
  CHECK(br.count_distinct_labels(5_obj, 5_obj, 0_lab, 9_lab) == 4);
  CHECK(br.count_distinct_labels(8_obj, 8_obj, 0_lab, 9_lab) == 6);
  CHECK(br.count_distinct_labels(11_obj, 11_obj, 0_lab, 9_lab) == 4);

  // Columns
  CHECK(br.count_distinct_labels(0_obj, 11_obj, 0_lab, 0_lab) == 1);
  CHECK(br.count_distinct_labels(0_obj, 11_obj, 2_lab, 2_lab) == 1);
  CHECK(br.count_distinct_labels(0_obj, 11_obj, 4_lab, 4_lab) == 1);
  CHECK(br.count_distinct_labels(0_obj, 11_obj, 5_lab, 5_lab) == 0);
  CHECK(br.count_distinct_labels(0_obj, 11_obj, 6_lab, 6_lab) == 1);
  CHECK(br.count_distinct_labels(0_obj, 11_obj, 9_lab, 9_lab) == 1);

  // Spread (all labels)
  CHECK(br.count_distinct_labels(0_obj, 3_obj, 0_lab, 9_lab) == 4);
  CHECK(br.count_distinct_labels(3_obj, 7_obj, 0_lab, 9_lab) == 9);
  CHECK(br.count_distinct_labels(9_obj, 11_obj, 0_lab, 9_lab) == 9);
  CHECK(br.count_distinct_labels(10_obj, 11_obj, 0_lab, 9_lab) == 7);
}

TEST_CASE("count_distinct_labels, complex") {
  const auto br = make_test_binary_relation();
  // Different wavelet tree nodes:
  CHECK(br.count_distinct_labels(1_obj, 3_obj, 0_lab, 7_lab) == 3);
  CHECK(br.count_distinct_labels(1_obj, 3_obj, 0_lab, 3_lab) == 1);
  CHECK(br.count_distinct_labels(1_obj, 3_obj, 4_lab, 7_lab) == 2);
  CHECK(br.count_distinct_labels(1_obj, 3_obj, 0_lab, 1_lab) == 0);
  CHECK(br.count_distinct_labels(1_obj, 3_obj, 2_lab, 3_lab) == 1);
  CHECK(br.count_distinct_labels(1_obj, 3_obj, 4_lab, 5_lab) == 1);
  CHECK(br.count_distinct_labels(1_obj, 3_obj, 6_lab, 7_lab) == 1);

  // Random label pairs
  auto count_labels = [&](const label_id alpha, const label_id beta) {
    return br.count_distinct_labels(1_obj, 4_obj, alpha, beta);
  };
  CHECK(count_labels(1_lab, 6_lab) == 3);
  CHECK(count_labels(2_lab, 9_lab) == 5);
  CHECK(count_labels(3_lab, 8_lab) == 4);
  CHECK(count_labels(0_lab, 5_lab) == 2);
  CHECK(count_labels(3_lab, 5_lab) == 1);
  CHECK(count_labels(7_lab, 9_lab) == 2);
  CHECK(count_labels(1_lab, 8_lab) == 5);
  CHECK(count_labels(0_lab, 5_lab) == 2);
  CHECK(count_labels(3_lab, 9_lab) == 4);
  CHECK(count_labels(4_lab, 8_lab) == 4);
}

TEST_CASE("serialize_and_load, basic") {
  const auto br = make_test_binary_relation();
  const std::string test_dump_basename = "../test/test_dump";
  std::ofstream outstream(test_dump_basename);
  br.serialize();
  outstream.close();

  binary_relation br_loaded();
  std::ifstream instream(test_dump_basename);
  br_loaded.load();
  instream.close();
  CHECK(br == br_loaded);
}

TEST_SUITE_END();
