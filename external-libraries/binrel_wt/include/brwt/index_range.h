#ifndef BRWT_INDEX_RANGE_H
#define BRWT_INDEX_RANGE_H

#include <brwt/common_types.h> // index_type, size_type
#include <cassert>             // assert

namespace brwt {

class index_range {
public:
  using index_type = types::index_type;
  using size_type = types::size_type;

  index_range() = default;

  index_range(index_type first, index_type last) noexcept
      : m_first{first}, m_size{last - first} {
    // preconditions
    assert(first <= last);
    // postconditions
    assert(begin() == first);
    assert(end() == last);
  }

  index_type begin() const noexcept {
    return m_first;
  }
  index_type end() const noexcept {
    return m_first + m_size;
  }
  size_type size() const noexcept {
    return m_size;
  }
  bool empty() const noexcept {
    return m_size == 0;
  }

private:
  index_type m_first{};
  size_type m_size{};
};

// Non-member helpers

inline auto begin(const index_range& range) {
  return range.begin();
}
inline auto end(const index_range& range) {
  return range.end();
}
inline auto size(const index_range& range) {
  return range.size();
}
inline auto empty(const index_range& range) {
  return range.empty();
}

} // namespace brwt

#endif // BRWT_INDEX_RANGE_H
