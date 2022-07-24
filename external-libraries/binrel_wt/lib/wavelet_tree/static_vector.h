#ifndef BINREL_WT_SRC_STATIC_VECTOR_H // NOLINT
#define BINREL_WT_SRC_STATIC_VECTOR_H

#include <cassert>     // assert
#include <cstddef>     // size_t
#include <type_traits> // aligned_storage
#include <utility>     // forward

namespace brwt {

template <typename T, std::size_t N>
class static_vector {
public:
  using size_type = std::size_t;
  using const_pointer = const T*;
  using const_reference = const T&;

public:
  static_vector() : m_size{0} {
    // m_data is uninitialized.
  }

  // These functions are deleted for now.
  static_vector(const static_vector&) = delete;
  static_vector(static_vector&&) = delete;
  static_vector& operator=(const static_vector&) = delete;
  static_vector& operator=(static_vector&&) = delete;

  ~static_vector() {
    // If T is trivially destructible the compiler should optimize this away.
    while (!empty()) {
      pop_back();
    }
  }

  template <typename... Args>
  void emplace_back(Args&&... args) {
    assert(m_size < N);
    new (m_data + m_size) T(std::forward<Args>(args)...); // NOLINT
    ++m_size;
  }

  void pop_back() {
    assert(!empty());
    (data()[--m_size]).~T();
  }

  const_reference operator[](size_type pos) const noexcept {
    assert(pos < m_size);
    return data()[pos];
  }

  const_reference back() const noexcept {
    assert(!empty());
    return data()[m_size - 1];
  }

  bool empty() const noexcept {
    return m_size == 0;
  }

  const_pointer data() const noexcept {
    return reinterpret_cast<const T*>(m_data); // NOLINT
  }

private:
  using aligned_storage_t =
      typename std::aligned_storage<sizeof(T), alignof(T)>::type;

  aligned_storage_t m_data[N];
  size_type m_size;
};

} // end namespace brwt

#endif // BINREL_WT_SRC_STATIC_VECTOR_H
