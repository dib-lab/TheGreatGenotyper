#include "doctest.h"
#include <brwt/int_vector.h>

#include <algorithm>   // all_of, sort, reverse, equal, find
#include <iterator>    // distance
#include <limits>      // numeric_limits
#include <stdexcept>   // domain_error
#include <type_traits> // is_signed, is_unsigned, ...
#include <type_traits> // is_same

// TODO(Diego): Consider adding more static asserts.

using brwt::int_vector;
using value_t = int_vector::value_type;

template <typename T, typename U>
constexpr bool same = std::is_same<T, U>::value;

static auto std_vec(const int_vector& seq) {
  return std::vector<value_t>(seq.begin(), seq.end());
}

static auto std_vec(std::initializer_list<value_t> ilist) {
  return std::vector<value_t>(ilist);
}

// From here on, 'bpe' stands for 'bits per element'

static_assert(std::is_unsigned<int_vector::value_type>::value, "");
static_assert(std::is_signed<int_vector::size_type>::value, "");

static_assert(std::is_nothrow_default_constructible<int_vector>::value, "");
static_assert(std::is_copy_constructible<int_vector>::value, "");
static_assert(std::is_copy_assignable<int_vector>::value, "");
static_assert(std::is_nothrow_move_constructible<int_vector>::value, "");
static_assert(std::is_nothrow_move_assignable<int_vector>::value, "");
static_assert(std::is_nothrow_destructible<int_vector>::value, "");

TEST_SUITE("int_vector");

TEST_CASE("int_vector::int_vector()") {
  const int_vector v{};
  CHECK(v.size() == 0);
  CHECK(v.get_bpe() == 0);
  CHECK(v.allocated_bytes() == 0);
}

TEST_CASE("int_vector::int_vector(size_type, int)") {
  const int_vector v(/*count=*/10, /*bpe=*/3);
  CHECK(v.size() == 10);
  CHECK(v.get_bpe() == 3);
  CHECK(v.allocated_bytes() >= 8);

  auto is_zero = [](const auto value) { return value == 0; };
  CHECK(std::all_of(v.begin(), v.end(), is_zero));

  // bits per block
  constexpr auto bpb = std::numeric_limits<int_vector::value_type>::digits;

  CHECK_NOTHROW(int_vector(10, bpb - 1));
  CHECK_THROWS_AS(int_vector(10, bpb), std::domain_error);     // NOLINT
  CHECK_THROWS_AS(int_vector(10, bpb + 1), std::domain_error); // NOLINT
  CHECK_THROWS_AS(int_vector(10, bpb + 2), std::domain_error); // NOLINT
}

TEST_CASE("int_vector::int_vector(initializer_list<value_type>)") {
  SUBCASE("Empty initializer_list") {
    std::initializer_list<value_t> empty{};
    const int_vector seq = empty;
    CHECK(seq.size() == 0);
    CHECK(seq.get_bpe() == 0);
  }
  SUBCASE("Basic test 1") {
    const int_vector seq = {10, 20, 30, 40}; // 40 needs 6 bits.
    CHECK(seq.size() == 4);
    CHECK(seq.get_bpe() == 6);
    CHECK(std_vec(seq) == std_vec({10, 20, 30, 40}));
  }
  SUBCASE("Basic test 2") {
    const int_vector seq = {2, 4, 256, 3, 100, 255, 30}; // 256 needs 9 bits.
    CHECK(seq.size() == 7);
    CHECK(seq.get_bpe() == 9);
    CHECK(std_vec(seq) == std_vec({2, 4, 256, 3, 100, 255, 30}));
  }
  SUBCASE("All zero except one") {
    const int_vector seq = {0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
    CHECK(seq.size() == 10);
    CHECK(seq.get_bpe() == 1);
    CHECK(std_vec(seq) == std_vec({0, 0, 0, 0, 1, 0, 0, 0, 0, 0}));
  }
  SUBCASE("All zero") {
    const int_vector seq = {0, 0, 0, 0, 0, 0};
    CHECK(seq.size() == 6);
    CHECK(seq.get_bpe() == 1);
    CHECK(std_vec(seq) == std_vec({0, 0, 0, 0, 0, 0}));
  }
}

TEST_CASE("int_vector::set_value") {
  SUBCASE("No invalid values") {
    int_vector vec(/*count=*/100, /*bpe=*/10);
    // The maximum storable value is 1023.

    for (int i = 0; i < vec.size(); ++i) {
      vec[i] = static_cast<unsigned>(i * 10);
    }

    CHECK(vec[0] == 0);
    CHECK(vec[1] == 10);
    CHECK(vec[2] == 20);
    CHECK(vec[99] == 990);
    CHECK(vec[98] == 980);
    CHECK(vec[31] == 310);

    CHECK(vec[44] == 440);
    CHECK(vec[45] == 450);
    CHECK(vec[46] == 460);

    vec[45] = 1023;

    CHECK(vec[44] == 440);
    CHECK(vec[45] == 1023);
    CHECK(vec[46] == 460);

    static_assert(noexcept(vec[0] = 140), "");
  }
  SUBCASE("There are invalid values") {
    int_vector vec(/*count=*/50, /*bpe=*/9);
    // The maximum storable value is 511

    vec[30] = 511;

    CHECK(vec[29] == 0);
    CHECK(vec[30] == 511);
    CHECK(vec[31] == 0);

    CHECK(vec[29] == 0);
    CHECK(vec[30] == 511);
    CHECK(vec[31] == 0);

    CHECK(vec[40] == 0);

    vec[40] = 500;

    CHECK(vec[39] == 0);
    CHECK(vec[40] == 500);
    CHECK(vec[41] == 0);
  }
}

TEST_CASE("int_vector::front") {
  {
    int_vector seq = {10, 20, 30, 40};
    CHECK(seq.front() == 10);
    CHECK(seq[0] == 10);

    seq.front() = 42;
    CHECK(seq.front() == 42);
    CHECK(seq[0] == 42);
    CHECK(std_vec(seq) == std_vec({42, 20, 30, 40}));

    static_assert(same<decltype(seq.front()), int_vector::reference>, "");
  }
  {
    const int_vector seq = {10, 20, 30, 40};
    CHECK(seq.front() == 10);

    static_assert(same<decltype(seq.front()), int_vector::const_reference>, "");
  }
}

TEST_CASE("int_vector::back") {
  {
    int_vector seq = {10, 20, 30, 40};
    CHECK(seq.back() == 40);
    CHECK(seq[3] == 40);

    seq.back() = 13;
    CHECK(seq.back() == 13);
    CHECK(seq[3] == 13);
    CHECK(std_vec(seq) == std_vec({10, 20, 30, 13}));

    static_assert(same<decltype(seq.back()), int_vector::reference>, "");
  }
  {
    const int_vector seq = {10, 20, 30, 40};
    CHECK(seq.back() == 40);

    static_assert(same<decltype(seq.front()), int_vector::const_reference>, "");
  }
}

TEST_CASE("int_vector::size") {
  // int_vector(count, bpe)
  CHECK(int_vector(10, 14).size() == 10);
  CHECK(int_vector(10, 41).size() == 10);
  CHECK(int_vector(0, 14).size() == 0);
  CHECK(int_vector(42, 0).size() == 42);
  CHECK(int_vector(42, 25).size() == 42);
}

TEST_CASE("int_vector::empty()") {
  // int_vector(count, bpe)
  CHECK(int_vector(0, 0).empty());
  CHECK(int_vector(0, 10).empty());
  CHECK(!int_vector(10, 0).empty());
  CHECK(!int_vector(10, 10).empty());
}

TEST_CASE("int_vector::get_bpe()") {
  // int_vector(count, bpe)
  CHECK(int_vector(10, 14).get_bpe() == 14);
  CHECK(int_vector(30, 14).get_bpe() == 14);
  CHECK(int_vector(0, 12).get_bpe() == 12);
  CHECK(int_vector(21, 0).get_bpe() == 0);
  CHECK(int_vector(1, 32).get_bpe() == 32);
}

TEST_CASE("int_vector::allocated_bytes()") {
  CHECK(int_vector(10, 14).allocated_bytes() >= 24);
  CHECK(int_vector(0, 15).allocated_bytes() == 0);
  CHECK(int_vector(19, 0).allocated_bytes() == 0);
  CHECK(int_vector(0, 0).allocated_bytes() == 0);
  CHECK(int_vector(41, 37).allocated_bytes() >= 192);
}

// iterators

TEST_CASE("int_vector:: non-const iterators") {
  int_vector seq = {10, 20, 30, 40, 50};
  REQUIRE(seq.size() == 5);
  REQUIRE(seq.get_bpe() == 6); // 50 = 32 + 16 + 2 (max value is 63)

  CHECK(std::distance(seq.begin(), seq.end()) == 5);
  CHECK(seq.begin() + 5 == seq.end());
  CHECK(*seq.begin() == 10);
  CHECK(seq.end()[-1] == 50);

  std::reverse(seq.begin(), seq.end());
  CHECK(std_vec(seq) == std_vec({50, 40, 30, 20, 10}));

  *(seq.begin() + 2) = 63;
  *(seq.begin() + 4) = 42;

  CHECK(std_vec(seq) == std_vec({50, 40, 63, 20, 42}));

  std::sort(seq.begin(), seq.end());
  CHECK(std_vec(seq) == std_vec({20, 40, 42, 50, 63}));

  std::reverse(seq.begin(), seq.begin() + 3);
  CHECK(std_vec(seq) == std_vec({42, 40, 20, 50, 63}));
}

TEST_CASE("int_vector:: const iterators") {
  {
    int_vector seq = {10, 20, 30, 40};
    CHECK(std::distance(seq.cbegin(), seq.cend()) == 4);
    CHECK(std::equal(seq.begin(), seq.end(), seq.cbegin(), seq.cend()));
    CHECK(std::find(seq.cbegin(), seq.cend(), 30) == seq.begin() + 2);
  }

  const int_vector seq = {11, 22, 33, 44, 55};
  CHECK(std::distance(seq.begin(), seq.end()) == 5);
  CHECK(std::find(seq.begin(), seq.end(), 55) == seq.end() - 1);
  CHECK(std::find(seq.begin(), seq.end(), 42) == seq.end());
}

// modifiers

TEST_CASE("int_vector::clear") {
  int_vector seq = {10, 20, 30, 40};

  REQUIRE(seq.size() == 4);
  REQUIRE(seq.get_bpe() == 6); // 40 needs 6 bits

  seq.clear();
  CHECK(seq.size() == 0);
  CHECK(seq.get_bpe() == 0);
}

TEST_CASE("int_vector::erase(iterator)") {
  int_vector seq = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  REQUIRE(seq.size() == 10);

  // Erase the first element.
  auto it = seq.erase(seq.cbegin());
  CHECK(seq.size() == 9);
  CHECK(it == seq.begin());
  CHECK(*it == 2);
  CHECK(std_vec(seq) == std_vec({2, 3, 4, 5, 6, 7, 8, 9, 10}));

  // Erase the last element.
  it = seq.erase(seq.cend() - 1);
  CHECK(seq.size() == 8);
  CHECK(it == seq.end());
  CHECK(std_vec(seq) == std_vec({2, 3, 4, 5, 6, 7, 8, 9}));

  // Erase an element in the middle.
  it = seq.erase(seq.cbegin() + 4); // The element with value 6
  CHECK(seq.size() == 7);
  CHECK(it == seq.begin() + 4);
  CHECK(*it == 7);
  CHECK(std_vec(seq) == std_vec({2, 3, 4, 5, 7, 8, 9}));

  // Erase all.
  while (!seq.empty()) {
    it = seq.erase(seq.cend() - 1);
  }

  CHECK(seq.size() == 0);
  CHECK(it == seq.end());
  CHECK(std_vec(seq) == std_vec({}));

  // Finally, checks whether the bpe remains unchanged.
  CHECK(seq.get_bpe() == 4); // 9 uses 4 bits
}

namespace {
class erase_fixture {
protected:
  erase_fixture() {
    assert(seq.size() == 10);
  }
  ~erase_fixture() {
    REQUIRE(seq.get_bpe() == 4); // The bpe must not change.
  }

  int_vector seq = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int_vector::iterator it{};
};
} // end namespace

TEST_CASE_FIXTURE(erase_fixture,
                  "int_vector::erase(iter, iter): Remove from left") {

  it = seq.erase(seq.cbegin(), seq.cbegin() + 1); // Just one
  CHECK(seq.size() == 9);
  CHECK(it == seq.begin());
  CHECK(*it == 2);
  CHECK(std_vec(seq) == std_vec({2, 3, 4, 5, 6, 7, 8, 9, 10}));

  it = seq.erase(seq.cbegin(), seq.cbegin()); // Empty
  CHECK(seq.size() == 9);
  CHECK(it == seq.begin());
  CHECK(*it == 2);
  CHECK(std_vec(seq) == std_vec({2, 3, 4, 5, 6, 7, 8, 9, 10}));

  it = seq.erase(seq.cbegin(), seq.cbegin() + 3); // Three elements.
  CHECK(seq.size() == 6);
  CHECK(it == seq.begin());
  CHECK(*it == 5);
  CHECK(std_vec(seq) == std_vec({5, 6, 7, 8, 9, 10}));
}

TEST_CASE_FIXTURE(erase_fixture,
                  "int_vector::erase(iter, iter): Remove from right") {

  it = seq.erase(seq.cend() - 1, seq.cend()); // Just one
  CHECK(seq.size() == 9);
  CHECK(it == seq.end());
  CHECK(std_vec(seq) == std_vec({1, 2, 3, 4, 5, 6, 7, 8, 9}));

  it = seq.erase(seq.cend(), seq.cend()); // Empty
  CHECK(seq.size() == 9);
  CHECK(it == seq.end());
  CHECK(std_vec(seq) == std_vec({1, 2, 3, 4, 5, 6, 7, 8, 9}));

  it = seq.erase(seq.cend() - 3, seq.cend()); // Three elements
  CHECK(seq.size() == 6);
  CHECK(it == seq.end());
  CHECK(std_vec(seq) == std_vec({1, 2, 3, 4, 5, 6}));
}

TEST_CASE_FIXTURE(erase_fixture,
                  "int_vector::erase(iter, iter): Remove from center") {
  it = seq.erase(seq.cbegin() + 2, seq.cbegin() + 3);
  CHECK(seq.size() == 9);
  CHECK(it == seq.cbegin() + 2);
  CHECK(*it == 4);
  CHECK(std_vec(seq) == std_vec({1, 2, 4, 5, 6, 7, 8, 9, 10}));

  it = seq.erase(seq.cbegin() + 4, seq.cbegin() + 4); // Empty
  CHECK(seq.size() == 9);
  CHECK(it == seq.cbegin() + 4);
  CHECK(*it == 6);
  CHECK(std_vec(seq) == std_vec({1, 2, 4, 5, 6, 7, 8, 9, 10}));

  it = seq.erase(seq.cbegin() + 3, seq.cbegin() + 6); // Three elements
  CHECK(seq.size() == 6);
  CHECK(it == seq.cbegin() + 3);
  CHECK(*it == 8);
  CHECK(std_vec(seq) == std_vec({1, 2, 4, 8, 9, 10}));
}

TEST_CASE_FIXTURE(erase_fixture, "int_vector::erase(iter, iter): Remove all") {
  it = seq.erase(seq.cbegin(), seq.cend());
  CHECK(seq.size() == 0);
  CHECK(it == seq.begin());
  CHECK(it == seq.end());
  CHECK(seq.begin() == seq.end());
  CHECK(std_vec(seq) == std_vec({}));
}

TEST_CASE("int_vector:: member swap") {
  int_vector a = {10, 20, 30, 40};
  int_vector b = {1, 2, 3, 4, 5, 6, 7};

  REQUIRE(std_vec(a) == std_vec({10, 20, 30, 40}));
  REQUIRE(a.size() == 4);
  REQUIRE(a.get_bpe() == 6);
  REQUIRE(std_vec(b) == std_vec({1, 2, 3, 4, 5, 6, 7}));
  REQUIRE(b.size() == 7);
  REQUIRE(b.get_bpe() == 3);

  a.swap(b);

  CHECK(std_vec(a) == std_vec({1, 2, 3, 4, 5, 6, 7}));
  CHECK(a.size() == 7);
  CHECK(a.get_bpe() == 3);
  CHECK(std_vec(b) == std_vec({10, 20, 30, 40}));
  CHECK(b.size() == 4);
  CHECK(b.get_bpe() == 6);
}

TEST_CASE("int_vector:: non-member swap") {
  int_vector a = {10, 20, 30, 40};
  int_vector b = {1, 2, 3, 4, 5, 6, 7};

  REQUIRE(std_vec(a) == std_vec({10, 20, 30, 40}));
  REQUIRE(a.size() == 4);
  REQUIRE(a.get_bpe() == 6);
  REQUIRE(std_vec(b) == std_vec({1, 2, 3, 4, 5, 6, 7}));
  REQUIRE(b.size() == 7);
  REQUIRE(b.get_bpe() == 3);

  swap(a, b);

  CHECK(std_vec(a) == std_vec({1, 2, 3, 4, 5, 6, 7}));
  CHECK(a.size() == 7);
  CHECK(a.get_bpe() == 3);
  CHECK(std_vec(b) == std_vec({10, 20, 30, 40}));
  CHECK(b.size() == 4);
  CHECK(b.get_bpe() == 6);
}

namespace {
class int_vector_fixture {
protected:
  int_vector_fixture() {
    auto x = {10, 20, 30, 40};
    auto y = {10, 20, 30, 40, 50};

    std::copy(x.begin(), x.end(), a.begin());
    std::copy(x.begin(), x.end(), b.begin());
    std::copy(y.begin(), y.end(), c.begin());
    std::copy(y.begin(), y.end(), d.begin());
  }

  int_vector a = int_vector(/*count=*/4, /*bpe=*/10);
  int_vector b = int_vector(/*count=*/4, /*bpe=*/20);
  int_vector c = int_vector(/*count=*/5, /*bpe=*/10);
  int_vector d = int_vector(/*count=*/5, /*bpe=*/20);
};
} // end namespace

TEST_CASE_FIXTURE(int_vector_fixture, "operator==(int_vector,int_vector)") {
  CHECK(a == a);
  CHECK(a == b);
  CHECK_FALSE(a == c);
  CHECK_FALSE(a == d);

  CHECK(b == a);
  CHECK(b == b);
  CHECK_FALSE(b == c);
  CHECK_FALSE(b == d);

  CHECK_FALSE(c == a);
  CHECK_FALSE(c == b);
  CHECK(c == c);
  CHECK(c == d);

  CHECK_FALSE(d == a);
  CHECK_FALSE(d == b);
  CHECK(d == c);
  CHECK(d == d);
}

TEST_CASE_FIXTURE(int_vector_fixture, "operator!=(int_vector,int_vector)") {
  CHECK_FALSE(a != a);
  CHECK_FALSE(a != b);
  CHECK(a != c);
  CHECK(a != d);

  CHECK_FALSE(b != a);
  CHECK_FALSE(b != b);
  CHECK(b != c);
  CHECK(b != d);

  CHECK(c != a);
  CHECK(c != b);
  CHECK_FALSE(c != c);
  CHECK_FALSE(c != d);

  CHECK(d != a);
  CHECK(d != b);
  CHECK_FALSE(d != c);
  CHECK_FALSE(d != d);
}

TEST_CASE("swap(reference, reference)") {
  int_vector v = {10, 20, 30, 40};
  swap(v.front(), v.back());

  CHECK(v.front() == 40);
  CHECK(v.back() == 10);
}

TEST_CASE("swap(reference, value_type)") {
  int_vector v = {10, 20, 30, 40};
  value_t x = 42;
  swap(v.front(), x);

  CHECK(v.front() == 42);
  CHECK(x == 10);
}

TEST_CASE("swap(value_type, reference)") {
  int_vector v = {10, 20, 30, 40};
  value_t x = 42;
  swap(x, v.front());

  CHECK(v.front() == 42);
  CHECK(x == 10);
}

TEST_SUITE_END();
