#ifndef BRWT_BITMAP_H
#define BRWT_BITMAP_H

#include <brwt/bit_vector.h>   // bit_vector
#include <brwt/common_types.h> // index_type, size_type
#include <brwt/int_vector.h>   // int_vector
#include <cassert>             // assert
#include <iostream>            // iostream

namespace brwt {

class bitmap {
public:
  using index_type = brwt::index_type;
  using size_type = brwt::size_type;

public:
  bitmap() noexcept = default;
  explicit bitmap(bit_vector vec);

  bool access(index_type pos) const noexcept;

  size_type rank_0(index_type pos) const noexcept;
  size_type rank_1(index_type pos) const noexcept;

  index_type select_0(size_type nth) const noexcept;
  index_type select_1(size_type nth) const noexcept;

  size_type length() const noexcept; // TODO(Diego): Remove this.
  size_type size() const noexcept;

  size_type num_ones() const noexcept;
  size_type num_zeros() const noexcept;

  bool load(std::istream &in) noexcept;
  void serialize(std::ostream &out) const;

private:
  auto blocks_of_super_block(index_type sb_idx) const noexcept;
  size_type num_super_blocks() const noexcept;

  template <bool B>
  size_type num_of() const noexcept;

  template <bool B>
  size_type sb_rank(index_type sb_idx) const noexcept;

  template <bool B>
  size_type sb_exclusive_rank(index_type sb_idx) const noexcept;

  template <bool B>
  index_type sb_select(size_type nth) const noexcept;

  template <bool B>
  index_type select(size_type nth) const noexcept;

  /// Original bit sequence.
  bit_vector bit_seq;

  /// Number of set bits until the i-th super block (included).
  int_vector sb_rank_1;
};

// ==========================================
// Inline definitions
// ==========================================

inline auto bitmap::access(const index_type pos) const noexcept -> bool {
  return bit_seq.get(pos);
}

inline auto bitmap::length() const noexcept -> size_type {
  return bit_seq.length();
}

inline auto bitmap::size() const noexcept -> size_type {
  return bit_seq.length();
}

inline auto bitmap::num_ones() const noexcept -> size_type {
  if (size() == 0) {
    return 0;
  }
  assert(!sb_rank_1.empty());
  return static_cast<size_type>(sb_rank_1.back());
}

inline auto bitmap::num_zeros() const noexcept -> size_type {
  return size() - num_ones();
}

} // namespace brwt

#endif // BRWT_BITMAP_H
