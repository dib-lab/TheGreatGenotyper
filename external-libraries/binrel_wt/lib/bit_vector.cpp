#include <brwt/bit_vector.h>

#include <brwt/serialization.h>
#include <brwt/bit_hacks.h>       // lsb_mask
#include <brwt/utility.h>         // ceil_div
#include <algorithm>              // min
#include <cassert>                // assert

using brwt::bit_vector;
using size_type = bit_vector::size_type;
using block_type = bit_vector::block_type;

// ==========================================
// helper functions
// ==========================================

template <typename Container, typename Integer>
static constexpr decltype(auto) at(Container& c, Integer pos) noexcept {
  using index_t = typename Container::size_type;
  assert(static_cast<index_t>(pos) < c.size());
  return c[static_cast<index_t>(pos)];
}

static constexpr block_type make_mask(const size_type count) noexcept {
  using brwt::lsb_mask;
  return (count == bit_vector::bits_per_block)
             ? (block_type{0} - 1)
             : lsb_mask<block_type>(static_cast<int>(count));
}

// ==========================================
// bit_vector implementation
// ==========================================

bit_vector::bit_vector(const size_type count) : m_len{count} {
  using sz = decltype(blocks)::size_type;
  assert(count >= 0);

  const auto num_blocks = ceil_div(count, bits_per_block);
  blocks.resize(static_cast<sz>(num_blocks));
}

bit_vector::bit_vector(const size_type count, const block_type value)
    : bit_vector(count) {
  if (blocks.empty()) {
    return;
  }
  const auto mask = make_mask(std::min(bits_per_block, count));
  blocks[0] = value & mask;
}

bit_vector::bit_vector(const std::string& s)
    : bit_vector(static_cast<size_type>(s.length())) {
  auto rev_s = s.rbegin();
  for (size_type i = 0; i < m_len; ++i) {
    assert(rev_s[i] == '0' || rev_s[i] == '1');
    set(i, rev_s[i] == '1');
  }
}

size_type bit_vector::allocated_bytes() const noexcept {
  return static_cast<size_type>(blocks.capacity() * sizeof(block_type));
}

bool bit_vector::get(const size_type pos) const noexcept {
  const auto block = at(blocks, pos / bits_per_block);
  const auto mask = (block_type{1} << (pos % bits_per_block));
  return (block & mask) != 0;
}

void bit_vector::set(const size_type pos, const bool value) noexcept {
  auto& block = at(blocks, pos / bits_per_block);
  const auto mask = block_type{1} << (pos % bits_per_block);

  if (value) {
    block |= mask;
  } else {
    block &= ~mask;
  }
}

block_type bit_vector::get_chunk(const size_type pos,
                                 const size_type count) const noexcept {
  assert(count >= 0 && count <= bits_per_block);
  assert(pos >= 0 && pos + count <= length());

  // l prefix refers to the left (or front) bit.
  const auto lblock = pos / bits_per_block;
  const auto loffset = pos % bits_per_block;

  if (loffset + count <= bits_per_block) {
    const auto mask = make_mask(count);
    return (at(blocks, lblock) >> loffset) & mask;
  }

  const auto lcount = bits_per_block - loffset;
  assert(lcount < count);
  const auto rcount = count - lcount;

  const auto rmask = make_mask(rcount);
  auto result = at(blocks, lblock) >> loffset;
  result |= (at(blocks, lblock + 1) & rmask) << lcount;
  return result;
}

void bit_vector::set_chunk(const size_type pos, const size_type count,
                           const block_type value) noexcept {
  assert(count >= 0 && count <= bits_per_block);
  assert(pos >= 0 && pos + count <= length());

  const auto lblock = pos / bits_per_block;
  const auto loffset = pos % bits_per_block;

  if (loffset + count <= bits_per_block) {
    const auto mask = make_mask(count);
    at(blocks, lblock) &= ~(mask << loffset);
    at(blocks, lblock) |= (value & mask) << loffset;
    return;
  }

  const auto lcount = bits_per_block - loffset;
  assert(lcount < count);
  const auto rcount = count - lcount;

  const auto lmask = make_mask(lcount);
  const auto rmask = make_mask(rcount);

  at(blocks, lblock) &= ~(lmask << loffset);
  at(blocks, lblock) |= (value & lmask) << loffset;

  at(blocks, lblock + 1) &= ~rmask;
  at(blocks, lblock + 1) |= (value >> lcount) & rmask;
}

bool bit_vector::load(std::istream &in) noexcept {
  if (!in.good())
    return false;
  try {
    m_len = static_cast<size_type>(load_number(in));
    blocks = load_number_vector(in);
    return true;
  } catch (...) {
    return false;
  }
}

void bit_vector::serialize(std::ostream &out) const {
  if (!out.good())
    throw std::ostream::failure("Bad stream");
  serialize_number(out, static_cast<uint64_t>(m_len));
  serialize_number_vector(out, blocks);
}
