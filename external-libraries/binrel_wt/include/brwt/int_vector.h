#ifndef BRWT_INT_VECTOR_H
#define BRWT_INT_VECTOR_H

#include <brwt/bit_vector.h>      // bit_vector
#include <brwt/detail/iterator.h> // random_access_iterator
#include <algorithm>              // equal
#include <cassert>                // assert
#include <cstddef>                // ptrdiff_t
#include <initializer_list>       // initializer_list
#include <iostream>               // iostream

namespace brwt {

class int_vector {
public:
  using value_type = bit_vector::block_type;
  using size_type = bit_vector::size_type;
  using difference_type = std::ptrdiff_t;

  /// \brief Class to provide an l-value reference to a particular element from
  /// the array.
  class reference {
  private:
    reference(int_vector& vec, size_type pos) noexcept
        : vector{vec}, elem_pos{pos} {}

  public:
    reference(const reference&) = default;
    ~reference() = default;

    reference& operator=(const value_type value) noexcept {
      vector.set_value(elem_pos, value);
      return *this;
    }
    reference& operator=(const reference& other) noexcept {
      vector.set_value(elem_pos, other);
      return *this;
    }
    operator value_type() const noexcept { // NOLINT
      return vector.get_value(elem_pos);
    }

  private:
    int_vector& vector;
    size_type elem_pos;
    friend int_vector;
  };

  using const_reference = value_type;

  using iterator = detail::random_access_iterator<int_vector, value_type,
                                                  reference, difference_type>;

  using const_iterator =
      detail::random_access_iterator<const int_vector, value_type,
                                     const_reference, difference_type>;

public:
  /// \brief Constructs an empty vector.
  ///
  /// \par TIme complexity
  /// Constant.
  ///
  int_vector() = default;

  /// \brief Constructs a fixed-size vector with the given size requirements.
  ///
  /// All elements of the vector are zeroed.
  ///
  /// \param count The number of elements to store.
  /// \param bpe The number of bits per element to use.
  ///
  /// \par Time complexity
  /// Linear in <tt>count * bpe</tt>.
  ///
  /// \throws std::domain_error if \p bits is greater than or equal to the
  /// number of bits of <tt>value_type</tt>.
  ///
  int_vector(size_type count, int bpe);

  /// \brief Constructs the sequence with the given initializer list.
  ///
  /// \post If \p ilist is empty, \c get_bpe() will be equal to 0. Otherwise, if
  /// the maximum element from \p ilist is 0, \c get_bpe() will be equal to 1.
  /// Otherwise, \c get_bpe() will be equal to the number of bits used by the
  /// maximum element (as defined by \c brwt::used_bits()).
  ///
  /// \par Time complexity
  /// Linear in <tt>ilist.size()</tt>.
  ///
  int_vector(std::initializer_list<value_type> ilist);

  /// \name Element access
  /// @{

  /// \brief Returns a proxy to the element at the given position.
  ///
  /// \param pos The position of the element to access.
  ///
  /// \par Time complexity
  /// Constant.
  ///
  reference operator[](size_type pos) noexcept {
    return reference(*this, pos);
  }

  /// \brief Returns the value of the element at the given position.
  ///
  /// \param pos The position of the element to access.
  ///
  /// \par Time complexity
  /// Constant.
  ///
  const_reference operator[](size_type pos) const noexcept {
    return get_value(pos);
  }

  /// \brief Returns a \c reference to the first element.
  ///
  /// \pre <tt>!empty()</tt>
  ///
  reference front() noexcept {
    assert(!empty());
    return reference(*this, 0);
  }

  /// \brief Returns the value of the first element.
  ///
  /// \pre <tt>!empty()</tt>
  ///
  const_reference front() const noexcept {
    assert(!empty());
    return get_value(0);
  }

  /// \brief Returns a reference to the last element.
  ///
  /// \pre <tt>!empty()</tt>
  ///
  reference back() noexcept {
    assert(!empty());
    return reference(*this, size() - 1);
  }

  /// \brief Returns the value of the last element.
  ///
  /// \pre <tt>!empty()</tt>
  ///
  const_reference back() const noexcept {
    assert(!empty());
    return get_value(size() - 1);
  }

  /// @}

  /// \name Capacity
  /// @{

  /// \brief Returns the number of elements of the array.
  ///
  /// \par Time complexity
  /// Constant.
  ///
  size_type size() const noexcept {
    return num_elems;
  }

  /// \brief Checks whether the sequence is empty.
  ///
  bool empty() const noexcept {
    return size() == 0;
  }

  /// \brief Returns the number of bits per element.
  ///
  /// \par Time complexity
  /// Constant.
  ///
  size_type get_bpe() const noexcept {
    return bits_per_element;
  }

  /// \brief Returns the number of allocated bytes.
  ///
  /// \par Time complexity
  /// Constant.
  ///
  size_type allocated_bytes() const noexcept {
    return bit_seq.allocated_bytes();
  }

  /// @}

  /// \name Iterators
  /// @{

  iterator begin() noexcept {
    return iterator(*this, 0);
  }
  const_iterator begin() const noexcept {
    return const_iterator(*this, 0);
  }
  const_iterator cbegin() const noexcept {
    return const_iterator(*this, 0);
  }

  iterator end() noexcept {
    return iterator(*this, size());
  }
  const_iterator end() const noexcept {
    return const_iterator(*this, size());
  }
  const_iterator cend() const noexcept {
    return const_iterator(*this, size());
  }

  /// @}

  /// \name Modifiers
  /// @{

  /// \brief Clears the contents.
  ///
  /// \par Time complexity
  /// Constant.
  ///
  void clear() noexcept {
    // bit_seq is leaved as is.
    num_elems = 0;
    bits_per_element = 0;
  }

  /// \brief Removes the element at \p pos.
  ///
  /// \pre The iterator \p pos must be valid and dereferenceable.
  ///
  /// \returns Iterator following the element removed.
  ///
  /// \par Time complexity
  /// Linear in <tt>std::distance(pos, end())</tt>.
  ///
  iterator erase(const_iterator pos) noexcept;

  /// \brief Removes the elements in the range <tt>[first, last)</tt>.
  ///
  /// \returns Iterator following the last element removed.
  ///
  /// \par Time complexity
  /// Linear in <tt>std::distance(last, end())</tt>. Note that no element in the
  /// range needs to be destroyed as they are built-in integers (and therefore
  /// are trivially destructible).
  ///
  iterator erase(const_iterator first, const_iterator last) noexcept;

  /// \brief Swaps the contents.
  ///
  void swap(int_vector& other) noexcept {
    using std::swap;
    swap(bit_seq, other.bit_seq);
    swap(num_elems, other.num_elems);
    swap(bits_per_element, other.bits_per_element);
  }

  /// \brief Loads from input stream.
  ///
  bool load(std::istream &in) noexcept;

  /// \brief Serializes into output stream.
  ///
  void serialize(std::ostream &out) const;

  /// @}
private:
  value_type get_value(const size_type pos) const noexcept {
    assert(pos >= 0 && pos < num_elems && "Out of range");
    return bit_seq.get_chunk(pos * bits_per_element, bits_per_element);
  }
  void set_value(size_type pos, value_type value) noexcept;

  size_type index_of(const_iterator pos) const noexcept;
  iterator non_const(const_iterator pos) noexcept;

private:
  bit_vector bit_seq;
  size_type num_elems{};
  size_type bits_per_element{};
};

// ==========================================
// Non-member functions
// ==========================================

/// \brief Swaps the contents of \p lhs and \p rhs.
///
/// \relates int_vector
///
inline void swap(int_vector& lhs, int_vector& rhs) noexcept {
  lhs.swap(rhs);
}

/// \brief Compares \p lhs and \p rhs for equality.
///
/// Two vectors are equal if they have the same size and the same contents. The
/// number of bits per element is not taken into account.
///
/// \relates int_vector
///
inline bool operator==(const int_vector& lhs, const int_vector& rhs) {
  return lhs.size() == rhs.size() &&
         std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

/// \brief Compares \p lhs and \p rhs for inequality.
///
/// \returns <tt>!(lhs == rhs)</tt>
///
/// \relates int_vector
///
inline bool operator!=(const int_vector& lhs, const int_vector& rhs) {
  return !(lhs == rhs);
}

/// \brief Swaps the values of the elements that \p lhs and \p rhs are referring
/// to.
///
/// \relates int_vector::reference
///
inline void swap(int_vector::reference lhs, int_vector::reference rhs) {
  using value_t = int_vector::value_type;
  value_t tmp = lhs;
  lhs = rhs;
  rhs = tmp;
}

/// \brief Swaps the values of the elements that \p lhs and \p rhs are referring
/// to.
///
/// \relates int_vector::reference
///
inline void swap(int_vector::reference lhs, int_vector::value_type& rhs) {
  using value_t = int_vector::value_type;
  value_t tmp = lhs;
  lhs = rhs;
  rhs = tmp;
}

/// \brief Swaps the values of the elements that \p lhs and \p rhs are referring
/// to.
///
/// \relates int_vector::reference
///
inline void swap(int_vector::value_type& lhs, int_vector::reference rhs) {
  using value_t = int_vector::value_type;
  value_t tmp = lhs;
  lhs = rhs;
  rhs = tmp;
}

} // namespace brwt

#endif // BRWT_INT_VECTOR_H
