#include "doctest.h"
#include <brwt/detail/iterator.h>

#include <algorithm>   // find, reverse, sort
#include <cstddef>     // size_t, ptrdiff_t
#include <iterator>    // iterator_traits
#include <memory>      // unique_ptr
#include <type_traits> // is_same, is_convertible
#include <utility>     // pair, swap

// TODO(Diego): Review random_access_iterator in more detail when possible.

using brwt::detail::random_access_iterator;

namespace {
template <typename T>
class dyn_array {
public:
  // public types

  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T&;
  using const_reference = const T&;
  using iterator =
      random_access_iterator<dyn_array, T, reference, difference_type>;
  using const_iterator =
      random_access_iterator<const dyn_array, T, const_reference,
                             difference_type>;

  // constructors

  dyn_array() = default;
  explicit dyn_array(size_type count)
      : m_data{std::make_unique<T[]>(count)}, m_size{count} {}

  dyn_array(std::initializer_list<T> ilist) : dyn_array(ilist.size()) {
    std::copy(ilist.begin(), ilist.end(), begin());
  }

  // element access

  reference operator[](size_type pos) {
    return m_data[pos];
  }
  const_reference operator[](size_type pos) const {
    return m_data[pos];
  }

  // capacity

  size_type size() const {
    return m_size;
  }

  // iterators

  iterator begin() {
    return iterator(*this, 0);
  }
  const_iterator begin() const {
    return cbegin();
  }
  const_iterator cbegin() const {
    return const_iterator(*this, 0);
  }
  iterator end() {
    return iterator(*this, static_cast<difference_type>(size()));
  }
  const_iterator end() const {
    return cend();
  }
  const_iterator cend() const {
    return const_iterator(*this, static_cast<difference_type>(size()));
  }

private:
  std::unique_ptr<T[]> m_data{};
  size_type m_size{};
};
} // end namespace

// ===------------------------------------------------------===
//                  Static tests
// ===------------------------------------------------------===

namespace {
template <typename T, typename U>
constexpr bool same = std::is_same<T, U>::value;

template <typename T, typename... Args>
constexpr bool constructible = std::is_constructible<T, Args...>::value;

template <typename T, typename U>
constexpr bool assignable = std::is_assignable<T, U>::value;

template <typename From, typename To>
constexpr bool convertible = std::is_convertible<From, To>::value;
} // end namespace

template <typename Container>
static void test_iterators_types_and_conversions() {
  using iterator = typename Container::iterator;
  using const_iterator = typename Container::const_iterator;

  static_assert(!same<iterator, const_iterator>, "");

  static_assert(!constructible<iterator, const_iterator>, "");
  static_assert(constructible<const_iterator, iterator>, "");
  static_assert(!assignable<iterator, const_iterator>, "");
  static_assert(assignable<const_iterator, iterator>, "");

  static_assert(convertible<iterator, const_iterator>, "");
  static_assert(!convertible<const_iterator, iterator>, "");

  {
    Container c{};
    static_assert(same<decltype(c.begin()), iterator>, "");
    static_assert(same<decltype(c.cbegin()), const_iterator>, "");
    static_assert(same<decltype(c.end()), iterator>, "");
    static_assert(same<decltype(c.cend()), const_iterator>, "");
  }
  {
    const Container c{};
    static_assert(same<decltype(c.begin()), const_iterator>, "");
    static_assert(same<decltype(c.cbegin()), const_iterator>, "");
    static_assert(same<decltype(c.end()), const_iterator>, "");
    static_assert(same<decltype(c.cend()), const_iterator>, "");
  }
}

template <typename I>
void test_iterator_properties() {
  static_assert(std::is_nothrow_default_constructible<I>::value, "");
  static_assert(std::is_trivially_copy_constructible<I>::value, "");
  static_assert(std::is_trivially_move_constructible<I>::value, "");
  static_assert(std::is_trivially_copy_assignable<I>::value, "");
  static_assert(std::is_trivially_move_assignable<I>::value, "");
  static_assert(std::is_trivially_destructible<I>::value, "");

  static_assert(same<typename std::iterator_traits<I>::iterator_category,
                     std::random_access_iterator_tag>,
                "");

  using D = typename I::difference_type;

  I mut_it{}; // mutable iterator
  const I it{};
  const I it2{};
  const D d{};

  // element access
  static_assert(same<decltype(*it), typename I::reference>, "");
  static_assert(same<decltype(it[d]), typename I::reference>, "");
  static_assert(same<decltype(it.operator->()), typename I::pointer>, "");

  // modifiers
  static_assert(same<decltype(++mut_it), I&>, "");
  static_assert(same<decltype(mut_it++), I>, "");
  static_assert(same<decltype(--mut_it), I&>, "");
  static_assert(same<decltype(mut_it--), I>, "");
  static_assert(same<decltype(mut_it += d), I&>, "");
  static_assert(same<decltype(mut_it -= d), I&>, "");

  // arithmetic operators
  static_assert(same<decltype(it + d), I>, "");
  static_assert(same<decltype(d + it), I>, "");
  static_assert(same<decltype(it - d), I>, "");
  static_assert(same<decltype(it - it2), D>, "");

  // relational operators
  static_assert(same<decltype(it == it2), bool>, "");
  static_assert(same<decltype(it != it2), bool>, "");
  static_assert(same<decltype(it < it2), bool>, "");
  static_assert(same<decltype(it > it2), bool>, "");
  static_assert(same<decltype(it <= it2), bool>, "");
  static_assert(same<decltype(it >= it2), bool>, "");

  // Check swappable property.
  using std::swap;
  static_assert(noexcept(swap(mut_it, mut_it)), "");
}

// ===------------------------------------------------------===
//                  Test cases
// ===------------------------------------------------------===

TEST_SUITE("detail::random_access_iterator");

TEST_CASE("iterator: Static tests") {
  using container = dyn_array<int>;
  test_iterators_types_and_conversions<container>();
  test_iterator_properties<container::iterator>();
  test_iterator_properties<container::const_iterator>();
}

TEST_CASE("Element access: read only") {
  const dyn_array<int> seq = {10, 20, 30, 40};
  const auto b = seq.begin();
  const auto m = seq.begin() + 2;
  const auto e = seq.end();

  CHECK(*b == 10);
  CHECK(*m == 30);
  CHECK(*(e - 1) == 40);

  CHECK(b[0] == 10);
  CHECK(b[1] == 20);
  CHECK(b[2] == 30);
  CHECK(b[3] == 40);

  CHECK(m[-2] == 10);
  CHECK(m[-1] == 20);
  CHECK(m[0] == 30);
  CHECK(m[1] == 40);

  CHECK(e[-4] == 10);
  CHECK(e[-3] == 20);
  CHECK(e[-2] == 30);
  CHECK(e[-1] == 40);

  CHECK(m == std::find(b, e, 30));
}

TEST_CASE("Element access: Modify") {
  dyn_array<int> seq = {10, 20, 30, 40};
  const auto b = seq.begin();
  const auto e = seq.end();

  std::reverse(b, e); // Uses operator*()

  CHECK(b[0] == 40);
  CHECK(b[1] == 30);
  CHECK(b[2] == 20);
  CHECK(b[3] == 10);

  b[0] = 15;
  b[3] = 5;

  CHECK(b[0] == 15);
  CHECK(b[1] == 30);
  CHECK(b[2] == 20);
  CHECK(b[3] == 5);

  std::sort(b, e);

  CHECK(b[0] == 5);
  CHECK(b[1] == 15);
  CHECK(b[2] == 20);
  CHECK(b[3] == 30);
}

TEST_CASE("operator->()") {
  using p = std::pair<int, char>;
  dyn_array<p> seq = {{1, 'a'}, {2, 'b'}, {3, 'c'}, {4, 'd'}};

  auto i = seq.begin();
  CHECK(i->first == 1);
  CHECK(i->second == 'a');

  i->first = 10;
  CHECK(*i == p(10, 'a'));
  CHECK(i->first == 10);
  CHECK(i->second == 'a');

  i->second = 'A';
  CHECK(*i == p(10, 'A'));
  CHECK(i->first == 10);
  CHECK(i->second == 'A');

  for (; i != seq.end(); ++i) {
    i->first *= 11;
    i->second += 2;
  }

  CHECK(seq[0] == p(110, 'C'));
  CHECK(seq[1] == p(22, 'd'));
  CHECK(seq[2] == p(33, 'e'));
  CHECK(seq[3] == p(44, 'f'));
}

namespace {
class iterator_fixture {
  dyn_array<int> seq = {10, 20, 30};

protected:
  using iterator = decltype(seq)::iterator;
  using const_iterator = decltype(seq)::const_iterator;

  iterator_fixture() {
    assert(&*i0 == &seq[0]);
    assert(&*i1 == &seq[1]);
    assert(&*i2 == &seq[2]);

    assert(i0 == seq.begin());
    assert(i3 == seq.end());

    assert(i0 != i1);
    assert(i0 != i2);
    assert(i0 != i3);
    assert(i1 != i2);
    assert(i1 != i3);
    assert(i2 != i3);
  }

  // member data
  const_iterator i{};

  const const_iterator i0 = seq.begin();
  const const_iterator i1 = seq.begin() + 1;
  const const_iterator i2 = seq.begin() + 2;
  const const_iterator i3 = seq.begin() + 3;

  const iterator j0 = seq.begin();
};
} // end namespace

// modifiers

TEST_CASE_FIXTURE(iterator_fixture, "pre-increment") {
  i = i0;
  CHECK(++i == i1);
  CHECK(++i == i2);
  CHECK(++i == i3);
  i = i0;
  CHECK(++(++i) == i2);
  (++i) = i0;
  CHECK(i == i0);
}

TEST_CASE_FIXTURE(iterator_fixture, "post-increment") {
  i = i0;
  CHECK(i++ == i0);
  CHECK(i == i1);
  CHECK(i++ == i1);
  CHECK(i == i2);
  CHECK(i++ == i2);
  CHECK(i == i3);
  i = i0;
  CHECK((i++)++ == i0);
  CHECK(i == i1);
  CHECK(++(i++) == i2);
  CHECK(i == i2);
}

TEST_CASE_FIXTURE(iterator_fixture, "pre-decrement") {
  i = i3;
  CHECK(--i == i2);
  CHECK(--i == i1);
  CHECK(--i == i0);
  i = i3;
  CHECK(--(--i) == i1);
  (--i) = i3;
  CHECK(i == i3);
}

TEST_CASE_FIXTURE(iterator_fixture, "post-decrement") {
  i = i3;
  CHECK(i-- == i3);
  CHECK(i == i2);
  CHECK(i-- == i2);
  CHECK(i == i1);
  CHECK(i-- == i1);
  CHECK(i == i0);
  i = i3;
  CHECK((i--)-- == i3);
  CHECK(i == i2);
  CHECK(--(i--) == i1);
  CHECK(i == i1);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator +=") {
  i = i0;
  CHECK((i += 2) == i2);
  CHECK((i += 1) == i3);
  i = i0;
  CHECK(((i += 1) += 1) == i2);
  CHECK((i += 0) == i2);
  (i += 1) = i0;
  CHECK(i == i0);
  i = i0;
  i += 3;
  CHECK(i == i3);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator -=") {
  i = i3;
  CHECK((i -= 2) == i1);
  CHECK((i -= 1) == i0);
  i = i3;
  CHECK(((i -= 1) -= 1) == i1);
  CHECK((i -= 0) == i1);
  (i -= 1) = i3;
  CHECK(i == i3);
  i = i3;
  i -= 3;
  CHECK(i == i0);
}

// arithmetic operators

TEST_CASE_FIXTURE(iterator_fixture, "operator+(iter,diff)") {
  CHECK(i0 + 0 == i0);
  CHECK(i0 + 1 == i1);
  CHECK(i0 + 2 == i2);
  CHECK(i0 + 3 == i3);

  CHECK(i1 + (-1) == i0);
  CHECK(i1 + 0 == i1);
  CHECK(i1 + 1 == i2);
  CHECK(i1 + 2 == i3);

  CHECK(i2 + (-2) == i0);
  CHECK(i2 + (-1) == i1);
  CHECK(i2 + 0 == i2);
  CHECK(i2 + 1 == i3);

  CHECK(i3 + (-3) == i0);
  CHECK(i3 + (-2) == i1);
  CHECK(i3 + (-1) == i2);
  CHECK(i3 + 0 == i3);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator+(diff,iter)") {
  CHECK(0 + i0 == i0);
  CHECK(1 + i0 == i1);
  CHECK(2 + i0 == i2);
  CHECK(3 + i0 == i3);

  CHECK((-1) + i1 == i0);
  CHECK(0 + i1 == i1);
  CHECK(1 + i1 == i2);
  CHECK(2 + i1 == i3);

  CHECK((-2) + i2 == i0);
  CHECK((-1) + i2 == i1);
  CHECK(0 + i2 == i2);
  CHECK(1 + i2 == i3);

  CHECK((-3) + i3 == i0);
  CHECK((-2) + i3 == i1);
  CHECK((-1) + i3 == i2);
  CHECK(0 + i3 == i3);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator-(iter,diff)") {
  CHECK(i0 - 0 == i0);
  CHECK(i0 - (-1) == i1);
  CHECK(i0 - (-2) == i2);
  CHECK(i0 - (-3) == i3);

  CHECK(i1 - 1 == i0);
  CHECK(i1 - 0 == i1);
  CHECK(i1 - (-1) == i2);
  CHECK(i1 - (-2) == i3);

  CHECK(i2 - 2 == i0);
  CHECK(i2 - 1 == i1);
  CHECK(i2 - 0 == i2);
  CHECK(i2 - (-1) == i3);

  CHECK(i3 - 3 == i0);
  CHECK(i3 - 2 == i1);
  CHECK(i3 - 1 == i2);
  CHECK(i3 - 0 == i3);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator-(iter,iter)") {
  CHECK(i0 - i0 == 0);
  CHECK(i0 - i1 == -1);
  CHECK(i0 - i2 == -2);
  CHECK(i0 - i3 == -3);

  CHECK(i1 - i0 == 1);
  CHECK(i1 - i1 == 0);
  CHECK(i1 - i2 == -1);
  CHECK(i1 - i3 == -2);

  CHECK(i2 - i0 == 2);
  CHECK(i2 - i1 == 1);
  CHECK(i2 - i2 == 0);
  CHECK(i2 - i3 == -1);

  CHECK(i3 - i0 == 3);
  CHECK(i3 - i1 == 2);
  CHECK(i3 - i2 == 1);
  CHECK(i3 - i3 == 0);
}

// relational operators

TEST_CASE_FIXTURE(iterator_fixture, "operator==") {
  CHECK(i0 == i0);
  CHECK(!(i0 == i1));
  CHECK(i0 == j0);

  CHECK(!(i1 == i0));
  CHECK(i1 == i1);
  CHECK(!(i1 == j0));

  CHECK(j0 == i0);
  CHECK(!(j0 == i1));
  CHECK(j0 == j0);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator!=") {
  CHECK_FALSE(i0 != i0);
  CHECK(i0 != i1);
  CHECK_FALSE(i0 != j0);

  CHECK(i1 != i0);
  CHECK_FALSE(i1 != i1);
  CHECK(i1 != j0);

  CHECK_FALSE(j0 != i0);
  CHECK(j0 != i1);
  CHECK_FALSE(j0 != j0);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator<") {
  CHECK_FALSE(i0 < i0);
  CHECK(i0 < i1);
  CHECK_FALSE(i0 < j0);

  CHECK_FALSE(i1 < i0);
  CHECK_FALSE(i1 < i1);
  CHECK_FALSE(i1 < j0);

  CHECK_FALSE(j0 < i0);
  CHECK(j0 < i1);
  CHECK_FALSE(j0 < j0);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator>") {
  CHECK_FALSE(i0 > i0);
  CHECK_FALSE(i0 > i1);
  CHECK_FALSE(i0 > j0);

  CHECK(i1 > i0);
  CHECK_FALSE(i1 > i1);
  CHECK(i1 > j0);

  CHECK_FALSE(j0 > i0);
  CHECK_FALSE(j0 > i1);
  CHECK_FALSE(j0 > j0);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator<=") {
  CHECK(i0 <= i0);
  CHECK(i0 <= i1);
  CHECK(i0 <= j0);

  CHECK_FALSE(i1 <= i0);
  CHECK(i1 <= i1);
  CHECK_FALSE(i1 <= j0);

  CHECK(j0 <= i0);
  CHECK(j0 <= i1);
  CHECK(j0 <= j0);
}

TEST_CASE_FIXTURE(iterator_fixture, "operator>=") {
  CHECK(i0 >= i0);
  CHECK_FALSE(i0 >= i1);
  CHECK(i0 >= j0);

  CHECK(i1 >= i0);
  CHECK(i1 >= i1);
  CHECK(i1 >= j0);

  CHECK(j0 >= i0);
  CHECK_FALSE(j0 >= i1);
  CHECK(j0 >= j0);
}

TEST_SUITE_END();
