#ifndef BRWT_WAVELET_TREE_WAVELET_TREE_H
#define BRWT_WAVELET_TREE_WAVELET_TREE_H

#include <brwt/bitmap.h>       // bitmap
#include <brwt/common_types.h> // symbol_id, size_type, index_type
#include <brwt/int_vector.h>   // int_vector
#include <cstddef>             // ptrdiff_t
#include <utility>             // pair
#include <iostream>            // iostream

namespace brwt {

/// \brief This class represents a wavelet tree.
///
/// A wavelet tree is used to manipulate sequences. It provides access, rank and
/// select on <tt>O(bits)</tt> time where \c bits is the number of bits used to
/// represent the symbols of the sequence.
///
/// This class also provide proxy nodes so that one can extend the wavelet tree
/// functionality easily regardless of the internal representation (which might
/// be optimized in some way).
///
class wavelet_tree {
public:
  class node_proxy;

public:
  /// \brief Constructs an empty wavelet tree.
  ///
  wavelet_tree() = default;

  /// \brief Constructs a wavelet tree from the given sequence.
  ///
  /// \param sequence The input sequence.
  ///
  /// \post <tt>get_bits_per_symbol() == sequence.get_bpe()</tt>
  ///
  /// \par Complexity
  /// Given:
  /// \li <tt>bpe = sequence.get_bpe()</tt>
  /// \li <tt>n = sequence.length()</tt>
  /// \li <tt>sigma = 2<sup>bpe</sup></tt>
  ///
  /// The time and space complexity to build the wavelet tree using this
  /// constructor is <tt>O(bpe * n + sigma)</tt>.
  ///
  explicit wavelet_tree(const int_vector& sequence);

  /// \brief Retrieves the symbol at the given position.
  ///
  /// \pre <tt>pos < size()</tt>
  ///
  symbol_id access(index_type pos) const noexcept;

  /// \brief Counts how many occurrences has a symbol up to the given position.
  ///
  /// \pre <tt>symbol <= max_symbol_id()</tt>
  /// \pre <tt>pos >= 0 && pos < size()</tt>
  ///
  size_type rank(symbol_id symbol, index_type pos) const noexcept;

  /// \brief Finds the position of the \e nth occurrence of the given symbol.
  ///
  /// \pre <tt>symbol <= max_symbol_id()</tt>
  /// \pre <tt>nth > 0</tt>
  ///
  /// \returns The position of the \e nth symbol if it exists. Otherwise returns
  /// <tt>-1</tt>.
  ///
  index_type select(symbol_id symbol, size_type nth) const noexcept;

  /// \brief Gets the size (or length) of the original sequence.
  ///
  size_type size() const noexcept {
    return seq_len;
  }

  /// \brief Gets the number of bits per symbol used in this wavelet tree.
  ///
  int get_bits_per_symbol() const noexcept;

  /// \brief Returns the maximum symbol id representable for this wavelet tree.
  ///
  symbol_id max_symbol_id() const noexcept;

  /// \brief Creates a proxy to the root node.
  ///
  /// The returned node proxy allows navigating through the wavelet tree.
  ///
  node_proxy make_root() const noexcept;

  /// \brief Loads wavelet_tree from input stream.
  ///
  bool load(std::istream &in) noexcept;

  /// \brief Serializes wavelet_tree into output stream.
  ///
  void serialize(std::ostream &out) const;

private:
  /// Representation of the wavelet tree without pointers.
  bitmap table{};

  /// The length of the original sequence.
  size_type seq_len{};

  /// The number of bits per symbol used in <tt>*this</tt>.
  size_type bits_per_symbol{};
};

/// \brief Proxy class to access the nodes of a wavelet tree.
///
/// Note that the node proxies allow exploring the wavelet tree regardless
/// of the internal structure.
///
class wavelet_tree::node_proxy {
public:
  /// \brief Constructs a proxy to the root node of the given wavelet tree.
  ///
  explicit node_proxy(const wavelet_tree& wt) noexcept;

  // internal bitmap access

  /// \brief Retrieves the specified bit from this node bitmap.
  ///
  bool access(size_type pos) const noexcept;

  /// \brief Invokes \c rank_0 on this node bitmap.
  ///
  size_type rank_0(index_type pos) const noexcept;

  /// \brief Invokes \c rank_1 on this node bitmap.
  ///
  size_type rank_1(index_type pos) const noexcept;

  /// \brief Invokes \c select_0 on this node bitmap.
  ///
  index_type select_0(size_type nth) const noexcept;

  /// \brief Invokes \c select_1 on this node bitmap.
  ///
  index_type select_1(size_type nth) const noexcept;

  /// \brief Retrieves the size of this node bitmap.
  ///
  size_type size() const noexcept {
    return range_size;
  }

  // Level information

  /// \brief Checks whether the node has no materialized children.
  ///
  bool is_leaf() const noexcept {
    return level_mask == static_cast<symbol_id>(1);
  }

  /// \brief Checks if the next bit of the symbol is handled by the left child.
  ///
  /// Note that this function can be called even on a leaf node (the node that
  /// manages the last division).
  ///
  bool is_lhs_symbol(symbol_id symbol) const noexcept {
    return (symbol & level_mask) == 0;
  }

  /// \brief Checks if the next bit of the symbol is handled by the left child.
  ///
  /// Note that this function can be called even on a leaf node (the node that
  /// manages the last division).
  ///
  bool is_rhs_symbol(symbol_id symbol) const noexcept {
    return !is_lhs_symbol(symbol);
  }

  // Navigation

  /// \brief Constructs a proxy to the left hand side child.
  ///
  /// \pre <tt>!is_leaf()</tt>
  ///
  node_proxy make_lhs() const noexcept;

  /// \brief Constructs a proxy to the right hand side child.
  ///
  /// \pre <tt>!is_leaf()</tt>
  ///
  node_proxy make_rhs() const noexcept;

  /// \brief Returns a pair containing proxies to the left child and to the
  /// right child.
  ///
  /// This function is equivalent to <tt>std::make_pair(make_lhs(),
  /// make_rhs())</tt>, except that it reduces the number of needed bitmap
  /// operations.
  ///
  /// \pre <tt>!is_leaf()</tt>
  ///
  /// \returns An \c std::pair where \c p.first refers to the left child and \c
  /// p.second refer to the second child.
  ///
  std::pair<node_proxy, node_proxy> make_lhs_and_rhs() const noexcept;

  /// \brief Checks if two node proxies refer to the same node.
  ///
  friend bool operator==(const node_proxy& lhs,
                         const node_proxy& rhs) noexcept {
    return lhs.wt_ptr == rhs.wt_ptr &&           //
           lhs.range_begin == rhs.range_begin && //
           lhs.level_mask == rhs.level_mask;
  }

  /// \brief Checks it twho proxies do not refer to the same node.
  ///
  friend bool operator!=(const node_proxy& lhs,
                         const node_proxy& rhs) noexcept {
    return !(lhs == rhs);
  }

private:
  // Memberwise constructor
  node_proxy(const wavelet_tree& wt_, index_type begin_, size_type size_,
             size_type ones_before_, symbol_id level_mask_) noexcept;

  // Absolute position information
  index_type begin() const noexcept;
  index_type end() const noexcept;
  size_type zeros_before() const noexcept;
  size_type ones_before() const noexcept;

  // Auxiliary methods
  size_type count_zeros() const noexcept;
  const bitmap& get_table() const noexcept;

private:
  const wavelet_tree* wt_ptr;
  index_type range_begin;
  size_type range_size;
  size_type num_ones_before; // equals to: get_table().rank_1(begin() - 1)
  symbol_id level_mask;
};

// ==========================================
// Extra inline definitions
// ==========================================

inline auto wavelet_tree::make_root() const noexcept -> node_proxy {
  return node_proxy(*this);
}

} // namespace brwt

#endif // BRWT_WAVELET_TREE_WAVELET_TREE_H
