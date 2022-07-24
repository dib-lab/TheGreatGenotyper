#ifndef BRWT_DETAIL_UTILITY_H
#define BRWT_DETAIL_UTILITY_H

#include <memory>      // addressof
#include <type_traits> // is_reference
#include <utility>     // move

namespace brwt {
namespace detail {

/// \brief Class used to treat a smart reference (or proxy) as a pointer.
///
template <typename SmartRef>
class pointed_reference {
  static_assert(!std::is_reference<SmartRef>::value,
                "This class is intended for smart references only");

public:
  using reference = const SmartRef&;
  using pointer = const SmartRef*;

  explicit pointed_reference(SmartRef ref) : m_ref{std::move(ref)} {}

  reference operator*() const noexcept {
    return m_ref;
  }

  pointer operator->() const noexcept {
    return std::addressof(m_ref);
  }

public:
  SmartRef m_ref;
};

template <typename SmartRef>
struct reference_traits {
  using reference = SmartRef;
  using pointer = pointed_reference<SmartRef>;

  static pointer as_pointer(SmartRef ref) {
    return pointer{std::move(ref)};
  }
};

template <typename T>
struct reference_traits<T&> {
  using reference = T&;
  using pointer = T*;

  static pointer as_pointer(T& ref) noexcept {
    return std::addressof(ref);
  }
};

} // namespace detail
} // namespace brwt

#endif // BRWT_DETAIL_UTILITY_H
