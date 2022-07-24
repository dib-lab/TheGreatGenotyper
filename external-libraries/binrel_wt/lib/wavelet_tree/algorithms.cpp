#include <brwt/wavelet_tree/algorithms.h>

#include "../generic_algorithms.h"          // int_binary_search
#include "bitmask_support.h"                // symbol_id-stuff
#include <brwt/wavelet_tree/wavelet_tree.h> // wavelet_tree class definition.
#include <algorithm>                        // min

namespace brwt {

using node_proxy = wavelet_tree::node_proxy;

// ==========================================
// Generic helpers
// ==========================================

template <typename T>
static const T& get_left(const std::pair<T, T>& p) {
  return p.first;
}

template <typename T>
static const T& get_right(const std::pair<T, T>& p) {
  return p.second;
}

// ==========================================
// node_proxy extensions
// ==========================================

static size_type inclusive_rank_0(const node_proxy& node,
                                  const index_type pos) noexcept {
  return node.rank_0(pos);
}

static size_type inclusive_rank_1(const node_proxy& node,
                                  const index_type pos) noexcept {
  return node.rank_1(pos);
}

static size_type exclusive_rank_0(const node_proxy& node,
                                  const index_type pos) noexcept {
  assert(pos >= 0 && pos <= node.size());
  if (pos == 0) {
    return 0;
  }
  return node.rank_0(pos - 1);
}

static size_type exclusive_rank_1(const node_proxy& node,
                                  const index_type pos) noexcept {
  assert(pos >= 0 && pos <= node.size());
  if (pos == 0) {
    return 0;
  }
  return node.rank_1(pos - 1);
}

// ==========================================
// index_range, node_proxy extensions
// ==========================================

static auto before_end(const index_range& range) noexcept {
  assert(size(range) > 0);
  return end(range) - 1;
}

static auto make_lhs_range(const index_range& range,
                           const node_proxy& node) noexcept {
  assert(!empty(range));
  assert(begin(range) >= 0 && end(range) <= node.size());

  return index_range{exclusive_rank_0(node, begin(range)),
                     inclusive_rank_0(node, before_end(range))};
}

static auto make_rhs_range(const index_range& range,
                           const node_proxy& node) noexcept {
  assert(!empty(range));
  assert(begin(range) >= 0 && end(range) <= node.size());

  return index_range{exclusive_rank_1(node, begin(range)),
                     inclusive_rank_1(node, before_end(range))};
}

static auto make_rhs_range_using_lhs(const index_range& range,
                                     const index_range& lhs_range) {
  assert(begin(range) >= begin(lhs_range));
  assert(end(range) >= end(lhs_range));
  return index_range{begin(range) - begin(lhs_range),
                     end(range) - end(lhs_range)};
}

static auto make_lhs_and_rhs_ranges(const index_range& range,
                                    const node_proxy& node) {
  assert(begin(range) >= 0 && end(range) <= node.size());

  const auto lhs = make_lhs_range(range, node);
  const auto rhs = make_rhs_range_using_lhs(range, lhs);
  return std::make_pair(lhs, rhs);
}

// ==========================================
// Wavelet Tree algorithms implementation
// ==========================================

size_type inclusive_rank(const wavelet_tree& wt, const symbol_id symbol,
                         const index_type pos) noexcept {
  return wt.rank(symbol, pos);
}

size_type exclusive_rank(const wavelet_tree& wt, const symbol_id symbol,
                         const index_type pos) noexcept {
  assert(pos >= 0 && pos <= wt.size());
  if (pos == 0) {
    return 0;
  }
  return inclusive_rank(wt, symbol, pos - 1);
}

size_type inclusive_rank(const wavelet_tree& wt,
                         const less_equal<symbol_id> cond,
                         index_type pos) noexcept {
  assert(pos >= 0 && pos < wt.size());

  const auto max_symbol = cond.max_value;
  size_type count = 0;
  auto node = wt.make_root();
  while (!node.is_leaf()) {
    if (node.is_lhs_symbol(max_symbol)) {
      pos = node.rank_0(pos) - 1;
      if (pos == -1) {
        return count;
      }
      node = node.make_lhs();
    } else {
      count += node.rank_0(pos);
      pos = node.rank_1(pos) - 1;
      if (pos == -1) {
        return count;
      }
      node = node.make_rhs();
    }
  }
  if (node.is_lhs_symbol(max_symbol)) {
    count += node.rank_0(pos);
  } else {
    count += (pos + 1);
  }
  return count;
}

size_type exclusive_rank(const wavelet_tree& wt,
                         const less_equal<symbol_id> cond,
                         index_type pos) noexcept {
  assert(pos >= 0 && pos <= wt.size());
  if (pos == 0) {
    return 0;
  }
  return inclusive_rank(wt, cond, pos - 1);
}

size_type inclusive_rank(const wavelet_tree& wt, between<symbol_id> cond,
                         const index_type pos) noexcept {
  return exclusive_rank(wt, cond, pos + 1);
}

size_type exclusive_rank(const wavelet_tree& wt, between<symbol_id> cond,
                         index_type end_pos) noexcept {
  if (cond.min_value == 0) {
    return exclusive_rank(wt, less_equal<symbol_id>{cond.max_value}, end_pos);
  }
  return exclusive_rank(wt, less_equal<symbol_id>{cond.max_value}, end_pos) -
         exclusive_rank(wt, less_equal<symbol_id>{prev(cond.min_value)},
                        end_pos);
}

size_type rank(const wavelet_tree& wt, const index_range range,
               const between<symbol_id> cond) noexcept {
  // TODO(Diego): It can be implemented faster.
  return exclusive_rank(wt, cond, end(range)) -
         exclusive_rank(wt, cond, begin(range));
}

namespace count_symbols_detail {

// TODO(diego): Optimization. Some overloads of count_symbols generate children
// nodes even when its respective ranges are empty.

static size_type count_symbols(node_proxy node, const index_range range) {
  assert(begin(range) >= 0 && end(range) <= node.size());

  if (empty(range)) {
    return 0;
  }

  const auto children_ranges = make_lhs_and_rhs_ranges(range, node);
  const auto& lhs_range = get_left(children_ranges);
  const auto& rhs_range = get_right(children_ranges);

  if (!node.is_leaf()) {
    if (empty(lhs_range)) {
      return count_symbols(node.make_rhs(), rhs_range);
    }
    if (empty(rhs_range)) {
      return count_symbols(node.make_lhs(), lhs_range);
    }
    const auto children = node.make_lhs_and_rhs();
    return count_symbols(get_left(children), lhs_range) +
           count_symbols(get_right(children), rhs_range);
  }

  return size_type{empty(lhs_range) ? 0 : 1} +
         size_type{empty(rhs_range) ? 0 : 1};
}

static size_type count_symbols(node_proxy node, index_range range,
                               const greater_equal<symbol_id> cond) noexcept {
  assert(begin(range) >= 0 && end(range) <= node.size());

  if (empty(range)) {
    return 0;
  }
  const auto min_symbol = cond.min_value;
  size_type count = 0;
  while (!node.is_leaf()) {
    const auto lhs_range = make_lhs_range(range, node);
    const auto rhs_range = make_rhs_range_using_lhs(range, lhs_range);

    if (node.is_lhs_symbol(min_symbol)) {
      const auto children = node.make_lhs_and_rhs(); // [lhs, rhs]
      count += count_symbols(get_right(children), rhs_range);
      node = get_left(children);
      range = lhs_range;
    } else {
      node = node.make_rhs();
      range = rhs_range;
    }

    if (empty(range)) {
      return count;
    }
  }
  assert(!empty(range));

  // At this point, the symbol of the rhs child is always >= min_symbol
  const auto rhs_range = make_rhs_range(range, node);
  if (!empty(rhs_range)) {
    count += 1;
  }
  if (node.is_lhs_symbol(min_symbol) && size(rhs_range) < size(range)) {
    count += 1; // the lhs range contains a valid symbol and is not empty.
  }
  return count;
}

static size_type count_symbols(node_proxy node, index_range range,
                               const less_equal<symbol_id> cond) noexcept {
  assert(begin(range) >= 0 && end(range) <= node.size());

  if (empty(range)) {
    return 0;
  }
  const auto max_symbol = cond.max_value;
  size_type count = 0;
  while (!node.is_leaf()) {
    const auto lhs_range = make_lhs_range(range, node);
    const auto rhs_range = make_rhs_range_using_lhs(range, lhs_range);

    if (node.is_rhs_symbol(max_symbol)) {
      const auto children = node.make_lhs_and_rhs(); // [lhs, rhs]
      count += count_symbols(get_left(children), lhs_range);
      node = get_right(children);
      range = rhs_range;
    } else {
      node = node.make_lhs();
      range = lhs_range;
    }

    if (empty(range)) {
      return count;
    }
  }
  assert(!empty(range));

  // At this point, the symbol of the lhs child is always <= max_symbol
  const auto lhs_range = make_lhs_range(range, node);
  if (!empty(lhs_range)) {
    count += 1;
  }
  if (node.is_rhs_symbol(max_symbol) && size(lhs_range) < size(range)) {
    count += 1; // the rhs range contains a valid symbol and is not empty.
  }
  return count;
}

static size_type count_symbols(node_proxy node, index_range range,
                               const symbol_id min_symbol,
                               const symbol_id max_symbol) noexcept {

  assert(begin(range) >= 0 && end(range) <= node.size());
  assert(min_symbol <= max_symbol);

  if (empty(range)) {
    return 0;
  }

  while (!node.is_leaf()) {
    if (node.is_lhs_symbol(max_symbol)) {
      range = make_lhs_range(range, node);
      node = node.make_lhs();
    } else if (node.is_rhs_symbol(min_symbol)) {
      range = make_rhs_range(range, node);
      node = node.make_rhs();
    } else {
      assert(node.is_lhs_symbol(min_symbol) && node.is_rhs_symbol(max_symbol));

      const auto lhs_range = make_lhs_range(range, node);
      const auto rhs_range = make_rhs_range_using_lhs(range, lhs_range);
      const auto children = node.make_lhs_and_rhs();

      return count_symbols(get_left(children), lhs_range,
                           greater_equal<symbol_id>{min_symbol}) +
             count_symbols(get_right(children), rhs_range,
                           less_equal<symbol_id>{max_symbol});
    }

    if (empty(range)) {
      return 0;
    }
  }
  assert(!empty(range));
  assert(node.is_leaf());

  if (node.is_lhs_symbol(max_symbol)) {
    return empty(make_lhs_range(range, node)) ? 0 : 1;
  }
  if (node.is_rhs_symbol(min_symbol)) {
    return empty(make_rhs_range(range, node)) ? 0 : 1;
  }
  assert(node.is_lhs_symbol(min_symbol) && node.is_rhs_symbol(max_symbol));

  const auto lhs_range = make_lhs_range(range, node);
  size_type result = 0;
  if (!empty(lhs_range)) {
    result += 1; // lhs range is not empty
  }
  if (size(lhs_range) < size(range)) {
    result += 1; // rhs range is not empty
  }
  return result;
}

} // namespace count_symbols_detail

size_type count_distinct_symbols(const wavelet_tree& wt,
                                 const index_range range,
                                 const less_equal<symbol_id> cond) noexcept {
  assert(begin(range) >= 0 && end(range) <= wt.size());
  assert(cond.max_value >= 0 && cond.max_value <= wt.max_symbol_id());

  return count_symbols_detail::count_symbols(wt.make_root(), range, cond);
}

size_type count_distinct_symbols(const wavelet_tree& wt,
                                 const index_range range,
                                 const greater_equal<symbol_id> cond) noexcept {
  assert(begin(range) >= 0 && end(range) <= wt.size());
  assert(cond.min_value >= 0 && cond.min_value <= wt.max_symbol_id());

  return count_symbols_detail::count_symbols(wt.make_root(), range, cond);
}

size_type count_distinct_symbols(const wavelet_tree& wt,
                                 const index_range range,
                                 const between<symbol_id> cond) noexcept {
  assert(begin(range) >= 0 && end(range) <= wt.size());
  assert(cond.min_value >= 0 && cond.max_value <= wt.max_symbol_id());
  namespace detail = count_symbols_detail;
  return detail::count_symbols(wt.make_root(), range, cond.min_value,
                               cond.max_value);
}

std::pair<symbol_id, index_type>
nth_element(const wavelet_tree& wt, index_range range, size_type nth) noexcept {

  assert(nth > 0 && nth <= size(range));
  const auto root_begin = begin(range);

  symbol_id symbol{};
  auto node = wt.make_root();
  while (!node.is_leaf()) {
    const auto lhs_range = make_lhs_range(range, node);
    const auto rhs_range = make_rhs_range_using_lhs(range, lhs_range);

    if (nth <= size(lhs_range)) {
      range = lhs_range;
      node = node.make_lhs();
      // nth -= 0
      // symbol |= 0;
    } else {
      range = rhs_range;
      node = node.make_rhs();
      nth -= size(lhs_range);
      symbol |= static_cast<symbol_id>(1u);
    }
    symbol <<= static_cast<symbol_id>(1);
  }
  {
    const auto lhs_range = make_lhs_range(range, node);
    if (nth > size(lhs_range)) {
      nth -= size(lhs_range);
      symbol |= static_cast<symbol_id>(1u);
    }
  }
  assert(nth > 0);

  const auto abs_nth = nth + exclusive_rank(wt, symbol, root_begin);
  return std::make_pair(symbol, wt.select(symbol, abs_nth));
}

index_type select(const wavelet_tree& wt, const between<symbol_id> cond,
                  const size_type nth) noexcept {
  // Find position such that:
  // inclusive_rank(wt, cond, pos) == nth;
  // exclusive_rank(wt, cond, pos) == nth - 1;
  auto pred = [&](const index_type pos) {
    return inclusive_rank(wt, cond, pos) < nth;
  };
  const auto pos = int_binary_search(index_type{0}, wt.size(), pred);
  if (pos == wt.size()) {
    return index_npos;
  }
  return pos;
}

namespace select_first_detail {

// index maps
static index_type make_lhs_start(const node_proxy& node,
                                 const index_type start) noexcept {
  return exclusive_rank_0(node, start);
}

static index_type make_rhs_start(const node_proxy& node,
                                 const index_type start) noexcept {
  return exclusive_rank_1(node, start);
}

static index_type select_first_0(const node_proxy& node,
                                 const index_type start) noexcept {
  // TODO(Jorge): Return node.select_next_0(pos) when available.
  assert(start >= 0 && start <= node.size());

  if (start == node.size()) {
    return index_npos;
  }
  if (!node.access(start)) {
    return start;
  }
  const auto nth = exclusive_rank_0(node, start) + 1;
  return node.select_0(nth);
}

static index_type select_first_1(const node_proxy& node,
                                 const index_type start) noexcept {
  // TODO(Jorge): Return node.select_next_1(pos) when available.
  assert(start >= 0 && start <= node.size());

  if (start == node.size()) {
    return index_npos;
  }
  if (node.access(start)) {
    return start;
  }
  const auto nth = exclusive_rank_1(node, start) + 1;
  return node.select_1(nth);
}

static index_type
leaf_select_first(const node_proxy& node, const index_type start,
                  const greater_equal<symbol_id> cond) noexcept {
  assert(node.is_leaf());

  if (node.is_lhs_symbol(cond.min_value)) {
    // All symbols are valid.
    return start;
  }
  return select_first_1(node, start);
}

static index_type leaf_select_first(const node_proxy& node,
                                    const index_type start,
                                    const less_equal<symbol_id> cond) noexcept {
  assert(node.is_leaf());

  if (node.is_rhs_symbol(cond.max_value)) {
    // All symbols are valid.
    return start;
  }
  return select_first_0(node, start);
}

static index_type leaf_select_first(const node_proxy& node,
                                    const index_type start,
                                    const between<symbol_id> cond) {
  assert(node.is_leaf());

  if (node.is_lhs_symbol(cond.max_value)) {
    return select_first_0(node, start);
  }
  if (node.is_rhs_symbol(cond.min_value)) {
    return select_first_1(node, start);
  }

  assert(node.is_lhs_symbol(cond.min_value) &&
         node.is_rhs_symbol(cond.max_value));
  return start;
}

// convenience functions

static std::pair<index_type, index_type>
make_lhs_and_rhs_start(const node_proxy& node,
                       const index_type start) noexcept {
  const auto lhs = exclusive_rank_0(node, start);
  const auto rhs = start - lhs;
  return {lhs, rhs};
}

static index_type remap_pos_from_lhs(const node_proxy& node,
                                     const index_type pos) noexcept {
  if (pos == index_npos) {
    return index_npos;
  }
  return node.select_0(pos + 1);
}

static index_type remap_pos_from_rhs(const node_proxy& node,
                                     const index_type pos) noexcept {
  if (pos == index_npos) {
    return index_npos;
  }
  return node.select_1(pos + 1);
}

// Computes the minimum between lhs and rhs, considering that lhs or rhs can
// be equal to index_npos. If lhs == index_npos, returns rhs. Otherwise, if
// rhs == index_npos, returns lhs. Otherwise, returns min(lhs, rhs).
static constexpr index_type min_index(const index_type lhs_pos,
                                      const index_type rhs_pos) noexcept {
  if (lhs_pos == index_npos) {
    return rhs_pos;
  }
  if (rhs_pos == index_npos) {
    return lhs_pos;
  }
  return std::min(lhs_pos, rhs_pos);
}

// ==========================
// select first implementation
// ==========================

static index_type select_first(const node_proxy& node, const index_type start,
                               const greater_equal<symbol_id> cond) noexcept {
  assert(start >= 0 && start <= node.size());

  if (start == node.size()) {
    return index_npos;
  }

  if (node.is_leaf()) {
    return leaf_select_first(node, start, cond);
  }

  if (node.is_rhs_symbol(cond.min_value)) {
    const auto rhs_first =
        select_first(node.make_rhs(), make_rhs_start(node, start), cond);
    return remap_pos_from_rhs(node, rhs_first);
  }

  const auto mapped_lhs_pos = [&] {
    const auto lhs_first =
        select_first(node.make_lhs(), make_lhs_start(node, start), cond);
    return remap_pos_from_lhs(node, lhs_first);
  }();
  const auto mapped_rhs_pos = select_first_1(node, start);

  return min_index(mapped_lhs_pos, mapped_rhs_pos);
}

static size_type select_first(const node_proxy& node, const index_type start,
                              const less_equal<symbol_id> cond) noexcept {
  assert(start >= 0 && start <= node.size());

  if (start == node.size()) {
    return index_npos;
  }

  if (node.is_leaf()) {
    return leaf_select_first(node, start, cond);
  }

  if (node.is_lhs_symbol(cond.max_value)) {
    const auto lhs_first =
        select_first(node.make_lhs(), make_lhs_start(node, start), cond);
    return remap_pos_from_lhs(node, lhs_first);
  }

  // TODO(jorge): consider optimizing by checking if the right child is empty.

  const auto mapped_lhs_pos = select_first_0(node, start);
  const auto mapped_rhs_pos = [&] {
    const auto rhs_first =
        select_first(node.make_rhs(), make_rhs_start(node, start), cond);
    return remap_pos_from_rhs(node, rhs_first);
  }();

  return min_index(mapped_lhs_pos, mapped_rhs_pos);
}

static index_type select_first(const node_proxy& node, const index_type start,
                               const between<symbol_id> cond) noexcept {
  assert(start >= 0 && start <= node.size());

  if (start == node.size()) {
    return index_npos;
  }

  if (node.is_leaf()) {
    return leaf_select_first(node, start, cond);
  }

  // Checks whether the interval [min_symbol, max_symbol] is fully covered by
  // either the left branch or the right branch.
  if (node.is_lhs_symbol(cond.max_value)) {
    const auto lhs_first =
        select_first(node.make_lhs(), make_lhs_start(node, start), cond);
    return remap_pos_from_lhs(node, lhs_first);
  }
  if (node.is_rhs_symbol(cond.min_value)) {
    const auto rhs_first =
        select_first(node.make_rhs(), make_rhs_start(node, start), cond);
    return remap_pos_from_rhs(node, rhs_first);
  }

  assert(node.is_lhs_symbol(cond.min_value) &&
         node.is_rhs_symbol(cond.max_value));

  // Merge results from both branches.
  const auto children = node.make_lhs_and_rhs();
  const auto children_start = make_lhs_and_rhs_start(node, start);

  const auto lhs_first =
      select_first(get_left(children), get_left(children_start),
                   greater_equal<symbol_id>{cond.min_value});
  const auto rhs_first =
      select_first(get_right(children), get_right(children_start),
                   less_equal<symbol_id>{cond.max_value});

  return min_index(remap_pos_from_lhs(node, lhs_first),
                   remap_pos_from_rhs(node, rhs_first));
}

} // end namespace select_first_detail

index_type select_first(const wavelet_tree& wt, index_type start,
                        between<symbol_id> cond) noexcept {
  assert(cond.min_value >= 0 && cond.min_value <= cond.max_value &&
         cond.max_value <= wt.max_symbol_id());
  return select_first_detail::select_first(wt.make_root(), start, cond);
}

} // end namespace brwt
