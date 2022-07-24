#ifndef BRWT_COMMON_TYPES_H
#define BRWT_COMMON_TYPES_H

#include <optional>              // optional
#include <cstddef>               // ptrdiff_t
#include <cstdint>               // uint_fast64_t

namespace brwt {
namespace types {
// TODO(Diego): Deprecate the use of 'types' namespace, remove the use of it
// along the library, then remove it. Types will be in brwt namespace now.

/// Signed vocabulary type used to indicate quantities.
using size_type = std::ptrdiff_t;

/// Signed vocabulary type used to indicate positions.
using index_type = size_type;

using word_type = std::uint_fast64_t;

} // end namespace types

// generic types

using types::index_type;
using types::size_type;
using types::word_type;

constexpr index_type index_npos = -1;

using std::nullopt;
using std::optional;

// wavelet tree types

/// Unsigned type used to represent symbol codes.
enum symbol_id : word_type {};

// condition types (used as function parameters)

template <typename T>
struct less_equal {
  T max_value;
};

template <typename T>
struct greater_equal {
  T min_value;
};

template <typename T>
struct between {
  T min_value;
  T max_value;
};

} // namespace brwt

#endif // BRWT_COMMON_TYPES_H
