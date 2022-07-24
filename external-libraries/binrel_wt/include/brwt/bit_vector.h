#ifndef BRWT_BIT_VECTOR_H
#define BRWT_BIT_VECTOR_H

#include <cstddef> // ptrdiff_t
#include <cstdint> // uint_fast64_t
#include <limits>  // numeric_limits
#include <string>  // string
#include <vector>  // vector

namespace brwt {

class bit_vector {
public:
  using size_type = std::ptrdiff_t;
  using block_type = std::uint_fast64_t;

public:
  static constexpr size_type bits_per_block =
      std::numeric_limits<block_type>::digits;

public:
  bit_vector() = default;
  explicit bit_vector(size_type count);
  explicit bit_vector(size_type count, block_type value);
  explicit bit_vector(const std::string& s);

  size_type length() const noexcept;
  size_type size() const noexcept;
  size_type num_blocks() const noexcept;
  size_type allocated_bytes() const noexcept;

  bool get(size_type pos) const noexcept;
  void set(size_type pos, bool value) noexcept;

  block_type get_chunk(size_type pos, size_type count) const noexcept;
  void set_chunk(size_type pos, size_type count, block_type value) noexcept;

  block_type get_block(size_type num_block) const noexcept;
  void set_block(size_type num_block, block_type value) noexcept;

  /// \brief Returns a pointer to the underlying array of blocks.
  ///
  /// The number of blocks can be obtained by calling \c num_blocks(). Note that
  /// the unused bits from the last block are set to zero.
  ///
  const block_type* data() const noexcept;

  bool load(std::istream &in) noexcept;
  void serialize(std::ostream &out) const;

private:
  size_type m_len{};
  std::vector<block_type> blocks;
};

// ==========================================
// Inline definitions
// ==========================================

inline auto bit_vector::length() const noexcept -> size_type {
  return m_len;
}

inline auto bit_vector::size() const noexcept -> size_type {
  return m_len;
}

inline auto bit_vector::num_blocks() const noexcept -> size_type {
  return static_cast<size_type>(blocks.size());
}

inline auto bit_vector::get_block(const size_type num_block) const noexcept
    -> block_type {
  return blocks[static_cast<std::size_t>(num_block)];
}

inline void bit_vector::set_block(const size_type num_block,
                                  const block_type value) noexcept {
  blocks[static_cast<std::size_t>(num_block)] = value;
}

inline auto bit_vector::data() const noexcept -> const block_type* {
  return blocks.data();
}

} // namespace brwt

#endif // BRWT_BIT_VECTOR_H
