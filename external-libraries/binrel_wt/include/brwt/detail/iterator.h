#ifndef BRWT_DETAIL_ITERATOR_H
#define BRWT_DETAIL_ITERATOR_H

#include <brwt/detail/utility.h> // reference_traits
#include <cassert>               // assert
#include <iterator>              // random_access_iterator_tag
#include <memory>                // addressof
#include <type_traits>           // remove_const_t

namespace brwt {
namespace detail {

/// \brief This class provide \e random-access iterators to containers that have
/// overloaded the subscript operator (\c operator[]).
///
/// For the sake of performance, two iterators can be compared if and only if
/// they refer to the same container, otherwise it is undefined behaviour. This
/// makes possible to compare them by using only its stored position.
///
/// \note This class fully satisfies the \c RandomAccessIterator concept.
///
template <typename Container, typename ValueType, typename Reference,
          typename DifferenceType>
class random_access_iterator {
public:
  // public types

  using value_type = ValueType;
  using reference = Reference;
  using pointer = typename reference_traits<reference>::pointer;
  using difference_type = DifferenceType;
  using iterator_category = std::random_access_iterator_tag;

  template <typename Ref>
  using non_const_iterator =
      random_access_iterator<std::remove_const_t<Container>, value_type, Ref,
                             difference_type>;

  // constructors

  constexpr random_access_iterator() = default;

  template <typename Ref>
  constexpr /*implicit*/ random_access_iterator( // NOLINT
      const non_const_iterator<Ref>& other) noexcept
      : m_cont{other.m_cont}, m_pos{other.m_pos} {}

  // element access

  reference operator*() const {
    using sz = typename Container::size_type;
    return (*m_cont)[static_cast<sz>(m_pos)];
  }

  pointer operator->() const {
    return reference_traits<reference>::as_pointer(*(*this));
  }

  reference operator[](difference_type n) const {
    using sz = typename Container::size_type;
    return (*m_cont)[static_cast<sz>(m_pos + n)];
  }

  // modifiers

  random_access_iterator& operator++() noexcept {
    ++m_pos;
    return *this;
  }
  random_access_iterator operator++(int) noexcept {
    auto old = *this;
    ++(*this);
    return old;
  }
  random_access_iterator& operator--() noexcept {
    --m_pos;
    return *this;
  }
  random_access_iterator operator--(int) noexcept {
    auto old = *this;
    --(*this);
    return old;
  }
  random_access_iterator& operator+=(difference_type d) noexcept {
    m_pos += d;
    return *this;
  }
  random_access_iterator& operator-=(difference_type d) noexcept {
    m_pos -= d;
    return *this;
  }

  // arithmetic operators

  friend random_access_iterator operator+(random_access_iterator it,
                                          difference_type d) noexcept {
    return it += d;
  }

  friend random_access_iterator operator+(difference_type d,
                                          random_access_iterator it) noexcept {
    return it += d;
  }

  friend random_access_iterator operator-(random_access_iterator it,
                                          difference_type d) noexcept {
    return it -= d;
  }

  friend difference_type operator-(const random_access_iterator& lhs,
                                   const random_access_iterator& rhs) noexcept {
    return lhs.m_pos - rhs.m_pos;
  }

  // relational operators

  friend bool operator==(const random_access_iterator& lhs,
                         const random_access_iterator& rhs) noexcept {
    assert(lhs.m_cont == rhs.m_cont);
    return lhs.m_pos == rhs.m_pos;
  }
  friend bool operator!=(const random_access_iterator& lhs,
                         const random_access_iterator& rhs) noexcept {
    return !(lhs == rhs);
  }
  friend bool operator<(const random_access_iterator& lhs,
                        const random_access_iterator& rhs) noexcept {
    assert(lhs.m_cont == rhs.m_cont);
    return lhs.m_pos < rhs.m_pos;
  }
  friend bool operator>(const random_access_iterator& lhs,
                        const random_access_iterator& rhs) noexcept {
    return rhs < lhs;
  }
  friend bool operator<=(const random_access_iterator& lhs,
                         const random_access_iterator& rhs) noexcept {
    return !(rhs < lhs);
  }
  friend bool operator>=(const random_access_iterator& lhs,
                         const random_access_iterator& rhs) noexcept {
    return !(lhs < rhs);
  }

private:
  constexpr random_access_iterator(Container& c, difference_type pos) noexcept
      : m_cont{std::addressof(c)}, m_pos{pos} {}

  // frienship to allow the container using the private constructor.
  friend Container;

  // frienship to allow non-const to const conversion.
  template <typename C, typename V, typename R, typename D>
  friend class random_access_iterator;

  // member data
  Container* m_cont{};
  difference_type m_pos{};
};

} // namespace detail
} // namespace brwt

#endif // BRWT_DETAIL_ITERATOR_H
