#include "doctest.h"
#include <brwt/detail/utility.h>

#include <brwt/detail/iterator.h> // random_access_iterator
#include <algorithm>              // max
#include <cassert>                // assert
#include <cstddef>                // size_t, ptrdiff_t
#include <iterator>               // distance, begin, end
#include <type_traits>            // conditional_t

// TODO(Diego): Review reference_traits in more detail when possible.

using brwt::detail::reference_traits;

namespace {
struct point {
  point(int x_, int y_) : x{x_}, y{y_} {}

  int x;
  int y;
};

bool operator==(const point& lhs, const point& rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y;
}

class positive_array;

template <typename Reference, bool IsConst>
class basic_iterator {
  using base_type =
      std::conditional_t<IsConst, const positive_array, positive_array>;

public:
  using difference_type = typename base_type::difference_type;
  using reference = Reference;
  using pointer = typename reference_traits<reference>::pointer;

  reference operator*() const {
    return (*m_array)[m_pos];
  }
  pointer operator->() const {
    return reference_traits<reference>::as_pointer(*(*this));
  }
  basic_iterator& operator++() {
    ++m_pos;
    return *this;
  }

  friend bool operator==(const basic_iterator& lhs, const basic_iterator& rhs) {
    assert(lhs.m_array == rhs.m_array);
    return lhs.m_pos == rhs.m_pos;
  }

private:
  basic_iterator(base_type* array, difference_type pos)
      : m_array{array}, m_pos{pos} {}
  friend positive_array;

  base_type* m_array;
  difference_type m_pos;
};

class positive_array {
public:
  using value_type = point;
  using size_type = std::ptrdiff_t;
  using difference_type = std::ptrdiff_t;
  class reference;
  using const_reference = const point&;
  using iterator = basic_iterator<reference, false>;
  using const_iterator = basic_iterator<const_reference, true>;

  reference operator[](size_type pos);

  const_reference operator[](size_type pos) const {
    return m_data[pos];
  }

  size_type size() const {
    return std::distance(std::begin(m_data), std::end(m_data));
  }

  iterator begin() {
    return iterator(this, 0);
  }
  const_iterator begin() const {
    return const_iterator(this, 0);
  }
  iterator end() {
    return iterator(this, size());
  }
  const_iterator end() const {
    return const_iterator(this, size());
  }

private:
  point m_data[5] = {{10, 1}, {20, 2}, {30, 3}, {40, 4}, {50, 5}};
};

class positive_array::reference {
public:
  reference& operator=(const point& rhs) {
    auto& lhs = m_point;
    lhs.x = std::max(lhs.x, rhs.x);
    lhs.y = std::max(lhs.y, rhs.y);
    return *this;
  }

  int get_x() const {
    return m_point.x;
  }

  int get_y() const {
    return m_point.y;
  }

  operator point() const {
    return m_point;
  }

private:
  reference(positive_array* array, size_type pos)
      : m_point{array->m_data[pos]} {}

  friend positive_array;

  point& m_point;
};

auto positive_array::operator[](size_type pos) -> reference {
  return reference(this, pos);
}

} // end namespace

// ===------------------------------------------------------===
//                      Tests
// ===------------------------------------------------------===

template <typename T, typename U>
constexpr bool same = std::is_same<T, U>::value;

TEST_SUITE("detail::reference_traits");

TEST_CASE("static checks") {
  using brwt::detail::pointed_reference;

  using reference = positive_array::reference;
  using const_reference = positive_array::const_reference;

  using pointer = typename reference_traits<reference>::pointer;
  using const_pointer = typename reference_traits<const_reference>::pointer;

  static_assert(same<pointer, pointed_reference<reference>>, "");
  static_assert(same<const_pointer, const point*>, "");
}

TEST_CASE("const_iterator test") {
  const positive_array array{};
  auto it = array.begin();

  // point m_data[5] = {{10, 1}, {20, 2}, {30, 3}, {40, 4}, {50, 5}};
  CHECK(*it == point(10, 1));
  CHECK(it->x == 10);
  CHECK(it->y == 1);
  ++it;
  CHECK(*it == point(20, 2));
  CHECK(it->x == 20);
  CHECK(it->y == 2);
  ++it;
  CHECK(&*it == &array[2]);
  CHECK(&(it->x) == &array[2].x);
  CHECK(&(it->y) == &array[2].y);
  ++it;
  CHECK(&*it == &array[3]);
  CHECK(&(it->x) == &array[3].x);
  CHECK(&(it->y) == &array[3].y);
  ++it;
  CHECK(*it == point(50, 5));
  CHECK(it->x == 50);
  CHECK(it->y == 5);
  ++it;
  CHECK(it == array.end());
}

TEST_CASE("iterator test") {
  positive_array array{};
  auto it = array.begin();

  CHECK(*it == point(10, 1));
  CHECK(it->get_x() == 10);
  CHECK(it->get_y() == 1);
  ++it;
  CHECK(*it == point(20, 2));
  CHECK(it->get_x() == 20);
  CHECK(it->get_y() == 2);

  *it = point(5, 50);
  CHECK(*it == point(20, 50));
  CHECK(it->get_x() == 20);
  CHECK(it->get_y() == 50);
  ++it;
  CHECK(*it == point(30, 3));
  *it = point(10, 10);
  CHECK(*it == point(30, 10));
  CHECK(it->get_x() == 30);
  CHECK(it->get_y() == 10);
  *++it = point(42, 0);
  CHECK(*it == point(42, 4));
  CHECK(it->get_x() == 42);
  CHECK(it->get_y() == 4);
  *++it = point(100, 150);
  CHECK(*it == point(100, 150));
  CHECK(it->get_x() == 100);
  CHECK(it->get_y() == 150);
  ++it;
  CHECK(it == array.end());
}

TEST_SUITE_END();
