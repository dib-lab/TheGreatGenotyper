#include <brwt/binary_relation.h>

#include <brwt/bit_hacks.h>   // used_bits
#include <brwt/index_range.h> // index_range
#include <algorithm>          // max, for_each, sort, unique
#include <cassert>            // assert
#include <cstddef>            // size_t
#include <iterator>           // begin, end
#include <tuple>              // tie
#include <type_traits>        // is_enum, underlying_type_t
#include <utility>            // pair, make_pair, move
#include <vector>             // vector

namespace brwt {

using std::size_t;
using node_proxy = wavelet_tree::node_proxy;

using object_id = binary_relation::object_id;
using label_id = binary_relation::label_id;
using pair_type = binary_relation::pair_type;

// ==========================================
// Generic helpers
// ==========================================

// TODO(Diego): Consider to put exclusive_scan in a private header. It is also
// used in the wavelet tree implementation.
template <typename InputIt, typename OutputIt, typename T>
static void exclusive_scan(InputIt first, const InputIt last, OutputIt d_first,
                           T init) {
  for (; first != last; ++first) {
    *d_first++ = std::exchange(init, init + (*first));
  }
}

template <typename InputRange, typename UnaryFunction>
static UnaryFunction for_each(InputRange&& range, UnaryFunction f) {
  return std::for_each(std::begin(range), std::end(range), f);
}

template <typename T>
static constexpr auto to_underlying_type(const T value) noexcept {
  static_assert(std::is_enum<T>::value);
  return static_cast<std::underlying_type_t<T>>(value);
}

// ==========================================
// binary_relation implementation
// ==========================================

static constexpr auto as_symbol(const label_id label) noexcept {
  return static_cast<symbol_id>(label);
}

static constexpr auto as_label(const symbol_id symbol) noexcept {
  return static_cast<label_id>(symbol);
}

static constexpr object_id prev(const object_id x) noexcept {
  assert(to_underlying_type(x) != 0);
  return static_cast<object_id>(to_underlying_type(x) - 1);
}

static constexpr label_id prev(const label_id alpha) noexcept {
  assert(to_underlying_type(alpha) != 0);
  return static_cast<label_id>(to_underlying_type(alpha) - 1);
}

static constexpr auto between_symbols(const label_id min,
                                      const label_id max) noexcept {
  assert(min <= max);
  return between<symbol_id>{as_symbol(min), as_symbol(max)};
}

static constexpr optional<pair_type>
optional_pair(const optional<object_id>& obj, const label_id label) {
  if (obj) {
    return pair_type{*obj, label};
  }
  return nullopt;
}

[[deprecated("when there is no pair with object_id=x, this function does not "
             "work as expected")]] auto
binary_relation::map(const object_id x) const noexcept -> index_type {
  // Note that in some cases upper_bound(x) == 0 so this function would return
  // -1 (an invalid index).
  return upper_bound(x) - 1;
}

[[deprecated(
    "the function 'map' will be removed, so this function would be "
    "unclear. Prefer get_associated_object() member function instead")]] auto
binary_relation::unmap(const index_type wt_pos) const noexcept -> object_id {
  return get_associated_object(wt_pos);
}

/// Returns the position of the first element in the wavelet tree such that its
/// associated object has value >= x.
/// In the stored wavelet tree, the range [lower_bound(x), wt.size()) contains
/// the labels of all the 'pairs' with object >= x (where 'pairs' refer to the
/// pairs of the original sequence).
/// If no such element exists, return m_wtree.size().
///
auto binary_relation::lower_bound(const object_id x) const noexcept
    -> index_type {
  if (x == 0) {
    return 0;
  }
  const auto flag_pos = m_bitmap.select_1(x);
  return (flag_pos + 1) - x; // m_bitmap.rank_0(flag_pos)
}

/// Returns the position of the first element in the wavelet tree such that its
/// associated object has value > x.
/// In the stored wavelet tree, the range [0, upper_bound(x)) contains the
/// labels of all the 'pairs' with object <= x (where 'pairs' refer to the pairs
/// of the original sequence).
/// If no such element exists, return m_wtree.size().
///
auto binary_relation::upper_bound(const object_id x) const noexcept
    -> index_type {
  const auto flag_pos = m_bitmap.select_1(x + 1);
  return (flag_pos + 1) - (x + 1); // m_bitmap.rank_0(flag_pos);
}

/// Returns the range of elements in the wavelet tree with (object == x)
///
auto binary_relation::equal_range(const object_id x) const noexcept {
  return index_range{lower_bound(x), upper_bound(x)};
}

/// Returns the range of elements with: (object >= x && object <= y)
///
auto binary_relation::make_mapped_range(const object_id x,
                                        const object_id y) const noexcept {
  // TODO(Diego): Consider renaming this member function.
  assert(x <= y);
  return index_range{lower_bound(x), upper_bound(y)};
}

/// Returns the associated object of the element in the wavelet tree at
/// position \p wt_pos.
///
auto binary_relation::get_associated_object(const index_type wt_pos) const
    noexcept -> object_id {
  assert(wt_pos >= 0 && wt_pos < m_wtree.size());

  const auto bit_pos = m_bitmap.select_0(wt_pos + 1);
  return object_id{(bit_pos + 1) - (wt_pos + 1)}; // m_bitmap.rank_1(bit_pos)
}

namespace pairs_constructor_detail {
using std::vector;

static auto count_objects_frequency(const vector<pair_type>& pairs,
                                    const object_id max_object) {
  // TODO(Diego): Consider to replace the histogram with a compressed vector.
  const auto num_objects = static_cast<size_t>(max_object) + 1;
  std::vector<size_type> frequency(num_objects);
  for_each(pairs, [&](const pair_type& p) { ++frequency[static_cast<size_t>(p.object)]; });
  return frequency;
}

template <typename ForwardRange, typename T>
static void inplace_exclusive_scan(ForwardRange& range, const T init) {
  exclusive_scan(std::begin(range), std::end(range), std::begin(range), init);
}

static wavelet_tree make_wavelet_tree(const vector<pair_type>& pairs,
                                      const label_id max_label,
                                      vector<size_type>& objects_frequency) {
  int_vector seq(/*count=*/static_cast<size_type>(pairs.size()),
                 /*bpe=*/used_bits(static_cast<word_type>(max_label)));

  // The first step is constructing the sequence of labels ordered by their
  // associated object value. The inplace_exclusive_scan and the for_each are an
  // inline counting sort to do this.
  inplace_exclusive_scan(objects_frequency, size_type{0});
  for_each(pairs, [&](const pair_type& p) {
    const auto next_pos = objects_frequency[static_cast<size_t>(p.object)]++;
    seq[next_pos] = static_cast<uint64_t>(p.label);
  });

  assert(objects_frequency.back() == seq.size());

  // Then, we have to sort each object range by label value and remove
  // duplicates. The one thing we know is the end of each object range (which is
  // stored in objects_frequency as a byproduct of the counting sort).
  auto first = seq.begin();
  auto seq_end = seq.begin();
  for_each(objects_frequency, [&](size_type& freq) {
    // First, remove duplicates. Note that at this point, freq contains the
    // accumulated frequency of all pairs until the current object.
    const auto last = seq.begin() + freq;
    std::sort(first, last);
    const auto unique_end = std::unique(first, last);
    seq_end = std::copy(first, unique_end, seq_end);

    // Then, use freq to keep the number of distinct pairs associated to the
    // current object (necessary for constructing the bitmap).
    freq = std::distance(first, unique_end);

    // Then, update first. Note that the end of the current object range is the
    // begin of the next object range.
    first = last;
  });

  // Finally, erase unused elements (because removing of duplicates) and
  // construct the wavelet tree.
  seq.erase(seq_end, seq.end());
  return wavelet_tree(seq);
}

static bitmap make_bitmap(const vector<size_type>& objects_frequency,
                          const size_type num_pairs) {
  auto bit_seq = [&] {
    const auto num_objects = static_cast<size_type>(objects_frequency.size());
    return bit_vector(/*count=*/num_pairs + num_objects);
  }();

  size_type acc_count = 0;
  for_each(objects_frequency, [&](const size_type count) {
    acc_count += count;
    bit_seq.set(acc_count, true);
    acc_count += 1;
  });
  return bitmap(std::move(bit_seq));
}

} // end namespace pairs_constructor_detail

binary_relation::binary_relation(const std::vector<pair_type>& pairs) {
  if (pairs.empty()) {
    return;
  }

  object_id max_object{};
  label_id max_label{};
  for_each(pairs, [&](const pair_type& pair) {
    max_object = std::max(max_object, pair.object);
    max_label = std::max(max_label, pair.label);
  });

  namespace detail = pairs_constructor_detail;
  auto objects_frequency = detail::count_objects_frequency(pairs, max_object);

  m_wtree = detail::make_wavelet_tree(pairs, max_label, objects_frequency);

  // Now that objects_frequency has been updated to contain the frequencies
  // of objects after ignoring duplicates pairs (done by make_wavelet_tree), we
  // can construct the bitmap.
  const auto num_unique_pairs = m_wtree.size();
  m_bitmap = detail::make_bitmap(objects_frequency, num_unique_pairs);

  // TODO(Diego): Assert for m_wtree.sigma() when available.
  assert(m_wtree.size() == num_unique_pairs);
  assert(m_bitmap.num_zeros() == num_unique_pairs);
  assert(m_bitmap.num_ones() == max_object + 1);
}

auto binary_relation::rank(object_id max_object, label_id max_label) const
    noexcept -> size_type {
  return exclusive_rank(m_wtree, less_equal<symbol_id>{as_symbol(max_label)},
                        upper_bound(max_object));
}

auto binary_relation::rank(object_id min_object, object_id max_object,
                           label_id max_label) const noexcept -> size_type {
  if (min_object == 0) {
    return rank(max_object, max_label);
  }
  return rank(max_object, max_label) - rank(prev(min_object), max_label);
}

auto binary_relation::rank(object_id max_object, label_id min_label,
                           label_id max_label) const noexcept -> size_type {
  if (min_label == 0) {
    return rank(max_object, max_label);
  }
  const auto cond =
      between<symbol_id>{as_symbol(min_label), as_symbol(max_label)};
  return exclusive_rank(m_wtree, cond, upper_bound(max_object));
}

auto binary_relation::nth_element(const object_id x, const object_id y,
                                  const label_id alpha, size_type nth,
                                  label_major_order_t /*unused*/) const noexcept
    -> optional<pair_type> {
  assert(x <= y);
  assert(nth > 0);

  if (alpha > 0) {
    nth += rank(x, y, prev(alpha));
  }

  const auto range = make_mapped_range(x, y);
  if (range.size() < nth) {
    return nullopt;
  }

  symbol_id symbol{};
  index_type wt_pos{};
  std::tie(symbol, wt_pos) = brwt::nth_element(m_wtree, range, nth);
  return pair_type{get_associated_object(wt_pos), as_label(symbol)};
}

auto binary_relation::nth_element(const object_id x, const label_id alpha,
                                  const label_id beta, const size_type nth,
                                  object_major_order_t /*unused*/) const
    noexcept -> optional<pair_type> {

  const auto first = lower_bound(x);
  const auto cond = between<symbol_id>{as_symbol(alpha), as_symbol(beta)};
  const auto wt_pos = [&] {
    const auto abs_nth = nth + exclusive_rank(m_wtree, cond, first);
    return brwt::select(m_wtree, cond, abs_nth);
  }();
  if (wt_pos == index_npos) {
    return nullopt;
  }

  const auto fixed_object = get_associated_object(wt_pos);
  assert(fixed_object >= x);

  const auto abs_nth = nth - [&] {
    const auto range = index_range{first, lower_bound(fixed_object)};
    return brwt::rank(m_wtree, range, cond);
  }();

  // TODO(Diego): It seems that several lower_bound, upper_bound pairs are
  // recomputed when some member functions use others. Consider implementing a
  // lazy propagation mechanism to avoid this.
  return nth_element(fixed_object, fixed_object, alpha, abs_nth, lab_major);
}

auto binary_relation::lower_bound(const pair_type start,
                                  const label_id min_label,
                                  const label_id max_label,
                                  object_major_order_t /*unused*/) const
    noexcept -> optional<pair_type> {
  assert(start.label >= min_label && start.label <= max_label);
  assert(min_label <= max_label);

  // TODO(Diego): The equal_range of start.object is computed again in the
  // call to nth_element(lab_major), and the upper_bound of start.object is
  // computed again in the call to nth_element(obj_major). Find a solution to
  // avoid this :)

  if (min_label == max_label) {
    // min_label (or max_label) is a fixed label.
    return optional_pair(obj_select(start.object, min_label, 1), min_label);
  }

  if (brwt::rank(m_wtree, equal_range(start.object),
                 between_symbols(start.label, max_label)) > 0) {
    return nth_element(start.object, start.object, start.label, 1, lab_major);
  }
  if (start.object + 1 == object_alphabet_size()) {
    return nullopt;
  }

  const auto wt_pos = select_first(m_wtree, upper_bound(start.object),
                                   between_symbols(min_label, max_label));
  if (wt_pos == index_npos) {
    return nullopt;
  }
  return pair_type{get_associated_object(wt_pos),
                   as_label(m_wtree.access(wt_pos))};
}

// ===------------------------------------------===
//                  Object view
// ===------------------------------------------===

auto binary_relation::obj_exclusive_rank(const object_id x,
                                         const label_id fixed_label) const
    noexcept -> size_type {
  return brwt::exclusive_rank(m_wtree, as_symbol(fixed_label), lower_bound(x));
}

auto binary_relation::obj_rank(const object_id x,
                               const label_id fixed_label) const noexcept
    -> size_type {
  return brwt::exclusive_rank(m_wtree, as_symbol(fixed_label), upper_bound(x));
}

auto binary_relation::obj_exclusive_rank(const object_id x,
                                         const label_id min_label,
                                         const label_id max_label) const
    noexcept -> size_type {

  assert(min_label <= max_label);
  return brwt::exclusive_rank(m_wtree, between_symbols(min_label, max_label),
                              lower_bound(x));
}

auto binary_relation::obj_rank(const object_id x, const label_id min_label,
                               const label_id max_label) const noexcept
    -> size_type {
  assert(min_label <= max_label);
  return brwt::exclusive_rank(m_wtree, between_symbols(min_label, max_label),
                              upper_bound(x));
}

auto binary_relation::obj_select(const object_id object_start,
                                 const label_id fixed_label,
                                 const size_type nth) const noexcept
    -> optional<object_id> {
  assert(nth > 0);

  const auto abs_nth = nth + obj_exclusive_rank(object_start, fixed_label);
  const auto wt_pos = m_wtree.select(as_symbol(fixed_label), abs_nth);

  if (wt_pos == index_npos) {
    return nullopt;
  }
  return get_associated_object(wt_pos);
}

auto binary_relation::count_distinct_labels(const object_id x,
                                            const object_id y,
                                            const label_id alpha,
                                            const label_id beta) const noexcept
    -> size_type {
  const auto range = make_mapped_range(x, y);
  const auto cond = between<symbol_id>{as_symbol(alpha), as_symbol(beta)};
  return count_distinct_symbols(m_wtree, range, cond);
}

bool binary_relation::load(std::istream &in) noexcept {
  if (!in.good())
    return false;
  try {
    if (!m_wtree.load(in)) {
      return false;
    }
    return m_bitmap.load(in);
  } catch (...) {
    return false;
  }
}

void binary_relation::serialize(std::ostream &out) const {
  if (!out.good())
    throw std::ostream::failure("Bad stream");
  m_wtree.serialize(out);
  m_bitmap.serialize(out);
}

} // end namespace brwt
