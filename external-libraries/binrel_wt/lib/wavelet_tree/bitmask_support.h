#ifndef BINREL_WT_LIB_WAVELET_TREE_BITMASK_SUPPORT_H // NOLINT
#define BINREL_WT_LIB_WAVELET_TREE_BITMASK_SUPPORT_H

#include <brwt/common_types.h> // symbol_id
#include <cassert>             // assert
#include <type_traits>         // underlying_type_t

namespace brwt {

// ==========================================
// symbol_id extensions
// ==========================================

constexpr auto extract_value(const symbol_id symbol) noexcept {
  return static_cast<std::underlying_type_t<symbol_id>>(symbol);
}

template <typename Integer>
constexpr auto operator<<(const symbol_id symbol, const Integer d) noexcept {
  assert(d >= 0);
  const auto res = extract_value(symbol) << d;
  // TODO(Diego): Use brace initialization when possible.
  return symbol_id(res);
}

template <typename Integer>
constexpr auto operator>>(const symbol_id symbol, const Integer d) noexcept {
  assert(d >= 0);
  return static_cast<symbol_id>(extract_value(symbol) >> d);
}

template <typename Integer>
constexpr auto& operator<<=(symbol_id& symbol, const Integer d) {
  assert(d >= 0);
  return symbol = symbol << d;
}

template <typename Integer>
constexpr auto& operator|=(symbol_id& symbol, const Integer value) {
  return symbol = static_cast<symbol_id>(extract_value(symbol) | value);
}

constexpr auto prev(const symbol_id symbol) noexcept {
  assert(extract_value(symbol) > 0);
  return static_cast<symbol_id>(extract_value(symbol) - 1);
}

} // end namespace brwt

#endif // header guad
