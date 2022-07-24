#ifndef BINREL_WT_LIB_UTILS_ARRAY_VIEW_H // NOLINT
#define BINREL_WT_LIB_UTILS_ARRAY_VIEW_H

#include <algorithm> // min
#include <cassert>   // assert
#include <cstddef>   // ptrdiff_t

namespace brwt {

template <typename T>
class array_view {
public:
  using const_pointer = const T*;
  using const_reference = const T&;
  using const_iterator = const_pointer;
  using size_type = std::ptrdiff_t;

  constexpr array_view() noexcept : m_data{nullptr}, m_size{0} {}

  constexpr array_view(const_pointer data, size_type count) noexcept
      : m_data{data}, m_size{count} {
    assert(count >= 0);
  }

  constexpr const_iterator begin() const noexcept {
    return m_data;
  }
  constexpr const_iterator end() const noexcept {
    return m_data + m_size;
  }
  constexpr size_type size() const noexcept {
    return m_size;
  }

  constexpr const_reference operator[](size_type pos) const noexcept {
    assert(pos >= 0 && pos < size());
    return m_data[pos];
  }

  constexpr array_view subarray(size_type pos, size_type count) const noexcept {
    assert(pos >= 0 && pos < size());
    assert(count >= 0);
    return array_view{m_data + pos, std::min(count, size() - pos)};
  }

private:
  const_pointer m_data;
  size_type m_size;
};

} // end namespace brwt

#endif // BINREL_WT_LIB_UTILS_ARRAY_VIEW_H
