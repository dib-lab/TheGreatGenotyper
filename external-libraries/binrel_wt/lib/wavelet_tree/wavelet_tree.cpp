#include <brwt/wavelet_tree.h>
#include <brwt/serialization.h> // serialize, load

#include "bitmask_support.h" // symbol_id-stuff
#include "static_vector.h"   // static_vector
#include <cassert>           // assert
#include <cstddef>           // size_t
#include <limits>            // numeric_limits
#include <utility>           // move, exchange
#include <vector>            // vector

using brwt::wavelet_tree;
using node_proxy = wavelet_tree::node_proxy;

// ==========================================
// Auxiliary algorithms
// ==========================================

template <typename InputIt, typename OutputIt, typename T>
static void exclusive_scan(InputIt first, const InputIt last, OutputIt d_first,
                           T init) {
  for (; first != last; ++first) {
    *d_first++ = std::exchange(init, init + (*first));
  }
}

// ==========================================
// wavelet_tree implementation
// ==========================================

wavelet_tree::wavelet_tree(const int_vector& sequence)
    : seq_len{sequence.size()}, bits_per_symbol{sequence.get_bpe()} {
  assert(bits_per_symbol >= 1);
  using std::size_t;
  using value_type = int_vector::value_type;

  // TODO(diegoramirez): Improve the constructor implementation. Consider
  // representing the tree with a Wavelet matrix.

  const auto alphabet_size = max_symbol_id() + 1;
  std::vector<size_type> next_pos(2 * alphabet_size);

  for (const auto symbol : sequence) {
    ++next_pos[alphabet_size + symbol];
  }
  {
    const auto first = next_pos.begin() + static_cast<size_type>(alphabet_size);
    exclusive_scan(first, next_pos.end(), first, size_type{0});
  }
  for (auto j = alphabet_size - 1; j > 0; --j) {
    const auto lhs = 2 * j; // rhs would be (2 * j + 1)
    next_pos[j] = next_pos[lhs];
  }

  // Finally we can fill the table.
  bit_vector bit_seq(bits_per_symbol * seq_len);

  auto push_symbol = [&](const value_type symbol) {
    value_type j = 1;
    value_type base_symbol = 0;
    auto num_symbols = alphabet_size;
    size_type level_pos = 0;

    while (num_symbols > 1) {
      const auto lhs_symbols = num_symbols / 2;
      const bool is_lhs = (symbol < base_symbol + lhs_symbols);

      bit_seq.set(level_pos + (next_pos[j]++), !is_lhs);

      if (is_lhs) {
        j = 2 * j;
        // base_symbol = base_symbol; // does not change
        num_symbols = lhs_symbols;
      } else {
        j = 2 * j + 1;
        base_symbol = (base_symbol + lhs_symbols);
        num_symbols = (num_symbols - lhs_symbols);
      }
      level_pos += seq_len;
    }
    assert(j == alphabet_size + symbol);
    assert(base_symbol == symbol);
    assert(num_symbols == 1);
    assert(level_pos == bit_seq.size());
  };
  for (const auto symbol : sequence) {
    push_symbol(symbol);
  }

  table = bitmap(std::move(bit_seq)); // The final magic.
}

auto wavelet_tree::access(index_type pos) const noexcept -> symbol_id {
  assert(pos >= 0 && pos < size());

  node_proxy node = make_root();
  symbol_id res{};
  while (!node.is_leaf()) {
    // each iteration invokes rank three times.
    if (!node.access(pos)) {
      // res |= 0;
      pos = node.rank_0(pos) - 1; // 1 rank
      node = node.make_lhs();     // 2 ranks
    } else {
      res |= static_cast<symbol_id>(1u);
      pos = node.rank_1(pos) - 1; // 1 rank
      node = node.make_rhs();     // 2 ranks
    }
    res <<= static_cast<symbol_id>(1);
  }
  res |= (node.access(pos)) ? static_cast<symbol_id>(1u) : static_cast<symbol_id>(0u);
  return res;
}

auto wavelet_tree::rank(const symbol_id symbol, index_type pos) const noexcept
    -> size_type {
  assert(symbol <= max_symbol_id());
  assert(pos >= 0 && pos < size());

  node_proxy node = make_root();
  while (!node.is_leaf()) {
    // each iteration invokes rank three times.
    if (node.is_lhs_symbol(symbol)) {
      pos = node.rank_0(pos) - 1; // 1 rank
      if (pos == -1) {
        return 0;
      }
      node = node.make_lhs(); // 2 ranks
    } else {
      pos = node.rank_1(pos) - 1; // 1 rank
      if (pos == -1) {
        return 0;
      }
      node = node.make_rhs(); // 2 ranks
    }
  }
  return node.is_lhs_symbol(symbol) ? node.rank_0(pos) : node.rank_1(pos);
}

auto wavelet_tree::select(const symbol_id symbol, const size_type nth) const
    noexcept -> index_type {
  assert(symbol <= max_symbol_id());
  assert(nth > 0);
  // Time complexity: Exactly 2 bitmap ranks and 1 bitmap select for level.

  constexpr size_t max_bits_per_symbol =
      std::numeric_limits<std::underlying_type_t<symbol_id>>::digits;
  static_vector<node_proxy, max_bits_per_symbol> stack;

  stack.emplace_back(make_root());
  while (true) {
    const auto& node = stack.back();
    if (node.size() < nth) {
      return -1;
    }
    if (node.is_leaf()) {
      break;
    }
    stack.emplace_back(node.is_lhs_symbol(symbol) ? node.make_lhs()
                                                  : node.make_rhs());
  }

  auto pos = [&] {
    const auto& node = stack.back();
    assert(node.size() >= nth);
    return node.is_lhs_symbol(symbol) ? node.select_0(nth) : node.select_1(nth);
  }();
  if (pos == -1) {
    return -1; // such element does not exists.
  }

  stack.pop_back();
  while (!stack.empty()) {
    const auto& node = stack.back();
    pos = node.is_lhs_symbol(symbol) ? node.select_0(pos + 1)
                                     : node.select_1(pos + 1);
    assert(pos >= 0 && pos < node.size()); // The element is supposed to exist.
    stack.pop_back();
  }
  return pos;
};

auto wavelet_tree::get_bits_per_symbol() const noexcept -> int {
  return static_cast<int>(bits_per_symbol);
}

auto wavelet_tree::max_symbol_id() const noexcept -> symbol_id {
  using limits = std::numeric_limits<symbol_id>;
  const auto res = bits_per_symbol == limits::digits
                       ? limits::max()
                       : (symbol_id(1) << bits_per_symbol) - 1;
  return static_cast<symbol_id>(res);
}

bool wavelet_tree::load(std::istream &in) noexcept {
  if (!in.good())
    return false;
  try {
    seq_len = static_cast<size_type>(load_number(in));
    bits_per_symbol = static_cast<size_type>(load_number(in));
    return table.load(in);
  } catch (...) {
    return false;
  }
}

void wavelet_tree::serialize(std::ostream &out) const {
  if (!out.good())
    throw std::ostream::failure("Bad stream");
  serialize_number(out, static_cast<uint64_t>(seq_len));
  serialize_number(out, static_cast<uint64_t>(bits_per_symbol));
  table.serialize(out);
}

// ==========================================
// node_proxy implementation
// ==========================================

node_proxy::node_proxy(const wavelet_tree& wt) noexcept
    : wt_ptr{&wt},
      range_begin{0},
      range_size{wt.seq_len},
      num_ones_before{0},
      level_mask{symbol_id(1) << (wt.bits_per_symbol - 1)} {
  assert(wt.bits_per_symbol >= 1);
}

auto node_proxy::access(const index_type pos) const noexcept -> bool {
  assert(pos >= 0 && pos < size());
  return get_table().access(begin() + pos);
}

// This function invokes table.rank_0 once.
auto node_proxy::rank_0(const index_type pos) const noexcept -> size_type {
  assert(pos >= 0 && pos < size());
  return get_table().rank_0(begin() + pos) - zeros_before();
}

// This function invokes table.rank_1 once.
auto node_proxy::rank_1(const index_type pos) const noexcept -> size_type {
  assert(pos >= 0 && pos < size());
  return get_table().rank_1(begin() + pos) - ones_before();
}

auto node_proxy::select_0(const size_type nth) const noexcept -> index_type {
  assert(nth > 0);
  const auto abs_pos = get_table().select_0(zeros_before() + nth);
  if (abs_pos == -1 || abs_pos >= end()) {
    return -1;
  }
  return abs_pos - begin();
}

auto node_proxy::select_1(const size_type nth) const noexcept -> index_type {
  assert(nth > 0);
  const auto abs_pos = get_table().select_1(ones_before() + nth);
  if (abs_pos == -1 || abs_pos >= end()) {
    return -1;
  }
  return abs_pos - begin();
}

// This function invokes table rank twice.
auto node_proxy::make_lhs() const noexcept -> node_proxy {
  assert(!is_leaf());
  const auto first = begin() + wt_ptr->seq_len;
  return node_proxy(/*wt_=*/*wt_ptr,
                    /*begin_=*/first,
                    /*size_=*/count_zeros(),
                    /*ones_before_=*/get_table().rank_1(first - 1),
                    /*level_mask_=*/(level_mask >> 1));
}

// This function invokes table rank twice.
auto node_proxy::make_rhs() const noexcept -> node_proxy {
  assert(!is_leaf());
  const auto num_zeros = count_zeros();
  const auto first = (begin() + wt_ptr->seq_len) + num_zeros;
  return node_proxy(/*wt_=*/*wt_ptr,
                    /*begin_=*/first,
                    /*size_=*/(size() - num_zeros),
                    /*ones_before_=*/get_table().rank_1(first - 1),
                    /*level_mask_=*/(level_mask >> 1));
}

// This function invokes table rank thrice.
auto node_proxy::make_lhs_and_rhs() const noexcept
    -> std::pair<node_proxy, node_proxy> {
  const auto num_zeros = count_zeros();
  const auto lhs_first = begin() + wt_ptr->seq_len;
  const auto rhs_first = lhs_first + num_zeros;
  return {node_proxy(
              /*wt_=*/*wt_ptr,
              /*begin_=*/lhs_first,
              /*size_=*/num_zeros,
              /*ones_before_=*/get_table().rank_1(lhs_first - 1),
              /*level_mask_=*/level_mask >> 1),
          node_proxy(
              /*wt_=*/*wt_ptr,
              /*begin_=*/rhs_first,
              /*size_=*/(size() - num_zeros),
              /*ones_before_=*/get_table().rank_1(rhs_first - 1),
              /*level_mask_=*/(level_mask >> 1))};
}

node_proxy::node_proxy(const wavelet_tree& wt_, const index_type begin_,
                       const size_type size_, const size_type ones_before_,
                       const symbol_id level_mask_) noexcept
    : wt_ptr{&wt_},
      range_begin{begin_},
      range_size{size_},
      num_ones_before{ones_before_},
      level_mask{level_mask_} {}

auto node_proxy::begin() const noexcept -> index_type {
  return range_begin;
}

auto node_proxy::end() const noexcept -> index_type {
  return begin() + size();
}

auto node_proxy::zeros_before() const noexcept -> size_type {
  return begin() - ones_before();
}

auto node_proxy::ones_before() const noexcept -> size_type {
  return num_ones_before;
}

auto node_proxy::count_zeros() const noexcept -> size_type {
  return rank_0(size() - 1);
}

auto node_proxy::get_table() const noexcept -> const bitmap& {
  return wt_ptr->table;
}
