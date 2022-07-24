/*----------------------------------------------------------------------------*\
|   Copyright (C) 2013 Nathaniel McClatchey                                    |
|   Released under the Boost Software License Version 1.0, which may be found  |
| at http://www.boost.org/LICENSE_1_0.txt                                      |
\*----------------------------------------------------------------------------*/

/*! @file interval_heap.hpp
//    interval_heap.hpp provides functions for creating interval heaps by
//  rearranging a range of elements. A valid interval heap has the following
//  properties:
//    Given elements A and B, A is less than B iff compare(A,B), where compare
//    is a binary functor describing a strict weak ordering.
//    Left and right bounds are those with even and odd indices, respectively.
//    A right bound is never less than its corresponding left bound.
//    A left bound is never less than the left bound of its parent.
//    A right bound is never less than the right bounds of its children.
//  @note Exception safety: All functions in interval_heap.hpp provide at
//  least the basic exception-safety guarantee. That is, no element will be lost
//  in the event of an exception. Elements may, however, no longer form an
//  interval heap.
//  @note Exception safety: Pre-C++11: If the comparison and swap don't throw,
//  neither do any of the interval-heap functions.
//  @note Exception safety: C++11: If the comparison, move-assign, and
//  move-construct, do not throw, neither do any of the interval-heap functions.
//  @note Thread safety: No static variables are modified by any operation.
//  @note Thread safety: Simultaneous read-write or write-write operations are
//  unsafe.
*/
#ifndef BOOST_HEAP_INTERVAL_HEAP_HPP_
#define BOOST_HEAP_INTERVAL_HEAP_HPP_

#ifndef __cplusplus
#error interval_heap.hpp requires a C++ compiler.
#endif

//  Grab std::swap and std::move (if available), while avoiding extra baggage.
#if (__cplusplus >= 201103L)
#include <utility>
#else
#include <algorithm>
#endif
//  Grab iterator data, for better use of templates.
#include <iterator>

//  Bulk-loading operation is well-suited to threading.
#ifndef BOOST_HEAP_INTERVAL_HEAP_USE_STD_THREAD
#define BOOST_HEAP_INTERVAL_HEAP_USE_STD_THREAD false
#endif
#if (BOOST_HEAP_INTERVAL_HEAP_USE_STD_THREAD == true)
#include <thread>
#include <future>
#ifndef BOOST_HEAP_INTERVAL_HEAP_AVAILABLE_THREADS
#define BOOST_HEAP_INTERVAL_HEAP_AVAILABLE_THREADS std::thread::hardware_concurrency()
#endif
#endif

namespace boost {
namespace heap {
//------------------------------User-Accessible---------------------------------
/// @brief Moves elements in [first,last) to form an interval heap.
/// @details This function moves elements in a range to from an interval-heap.
/// @param first,last A range of random-access iterators.
/// @param compare A comparison object.
/// @post [ @a first, @a last) is a valid interval heap.
/// @invariant No element is added to or removed from the range.
/// @note Complexity: O(n)
/// @note Exception safety: Basic if move/copy/swap provide the basic
///   guarantee.
/// @remark Threaded if <tt>BOOST_HEAP_INTERVAL_HEAP_USE_STD_THREAD</tt> is
///   defined as <tt>true</tt> (default <tt>false</tt>). Number of threads used
///   does not exceed <tt>BOOST_HEAP_INTERVAL_HEAP_AVAILABLE_THREADS</tt>.
template <typename Iterator, typename Compare>
void make_interval_heap (Iterator first, Iterator last, Compare compare);

/// @brief Expands the interval heap to include the element at last-1.
/// @details This function expands an interval heap from [ @a first, @a last -
/// 1) to [ @a first, @a last), maintaining the interval-heap property. This is
/// typically used to add a new element to the interval heap.
/// @param first,last A range of random-access iterators.
/// @param compare A comparison object.
/// @pre [ @a first, @a last - 1) is a valid interval heap.
/// @post [ @a first, @a last) is a valid interval heap.
/// @invariant No element is added to or removed from the range.
/// @note Complexity: O(log n)
/// @note Exception safety: Basic if move/copy/swap provide the basic
/// guarantee. Strong if they do not throw exceptions.
template <typename Iterator, typename Compare>
void push_interval_heap (Iterator first, Iterator last, Compare compare);

//!@{
/// @brief Moves a minimal element to the back of the range for popping.
/// @details This function moves a minimal element to the end of the range of
///   iterators. This is typically done so that the element can be efficiently
///   removed (by @a pop_back, for example), but can also be used to sort the
///   range in descending order.
/// @param first,last A range of random-access iterators.
/// @param compare A comparison object.
/// @pre [ @a first, @a last) is a valid interval heap.
/// @post A minimal element from the range has been moved to @a last - 1.
/// @post [ @a first, @a last - 1) is a valid interval heap.
/// @invariant No element is added to or removed from the range.
/// @note Complexity: O(log n)
/// @note Exception safety: Basic if move/copy/swap provide the basic
///   guarantee. Strong if they do not throw exceptions.
template <typename Iterator, typename Compare>
void pop_interval_heap_min (Iterator first, Iterator last, Compare compare);

/// @brief Moves a maximal element to the back of the range for popping.
/// @details This function moves a maximal element to the end of the range of
///   iterators. This is typically done so that the element can be efficiently
///   removed (by @a pop_back, for example), but can also be used to sort the
///   range in ascending order.
/// @param first,last A range of random-access iterators.
/// @param compare A comparison object.
/// @pre [ @a first, @a last) is a valid interval heap.
/// @post A maximal element from the range has been moved to @a last - 1.
/// @post [ @a first, @a last - 1) is a valid interval heap.
/// @invariant No element is added to or removed from the range.
/// @note Complexity: O(log n)
/// @note Exception safety: Basic if move/copy/swap provide the basic
///   guarantee. Strong if they do not throw exceptions.
template <typename Iterator, typename Compare>
void pop_interval_heap_max (Iterator first, Iterator last, Compare compare);

/// @brief Moves an element to the back of the range for popping.
/// @details This function moves a specified element to the end of the range of
///   iterators. This is typically done so that the element can be efficiently
///   removed (by @a pop_back, for example).
/// @param first,last A range of random-access iterators.
/// @param index The offset, from @a first, of the element to be removed.
/// @param compare A comparison object.
/// @pre [ @a first, @a last) is a valid interval heap.
/// @pre @a index is in the range [ 0, @a last - @a first)
/// @post [ @a first, @a last - 1) is a valid interval heap.
/// @post The element originally at @a first + @a index has been moved to
///   @a last - 1.
/// @invariant No element is added to or removed from the range.
/// @note Complexity: O(log n)
/// @note Exception safety: Basic if move/copy/swap provide the basic
///   guarantee. Strong if they do not throw exceptions.
template <typename Iterator, typename Compare>
void pop_interval_heap (Iterator first, Iterator last,
          typename std::iterator_traits<Iterator>::difference_type index,
          Compare compare);
//!@}

/// @brief Restores the interval-heap property violated by a single element at
///   an arbitrary position.
/// @param first,last A range of random-access iterators.
/// @param index The offset, from @a first, of the offending element.
/// @param compare A comparison object.
/// @pre @a index is in the range [ 0, @a last - @a first).
/// @pre [ @a first, @a last) is a valid interval heap, if the element at
/// @a first + @a index is ignored.
/// @post [ @a first, @a last) is a valid interval heap.
/// @invariant No element is added to or removed from the range.
/// @note Complexity: O(log n)
/// @note Exception safety: Basic if move/copy/swap provide the basic
///   guarantee. Strong if they do not throw exceptions.
template <typename Iterator, typename Compare>
void update_interval_heap (Iterator first, Iterator last,
           typename std::iterator_traits<Iterator>::difference_type index,
           Compare compare);

/// @brief Sorts an interval heap in ascending order.
/// @details This function takes an interval heap and sorts its elements in
///   ascending order.
/// @param first,last A range of random-access iterators.
/// @param compare A comparison object.
/// @pre [ @a first, @a last) is a valid interval heap.
/// @post [ @a first, @a last) is sorted in ascending order.
/// @invariant No element is added to or removed from the range.
/// @note Complexity: O(n log n)
/// @note Exception safety: Basic
template <typename Iterator, typename Compare>
void sort_interval_heap (Iterator first, Iterator last, Compare compare);

/// @brief Finds the largest subrange that qualifies as an interval heap.
/// @note Complexity: O(n)
/// @note Exception safety: No-throw
template <typename Iterator, typename Compare>
Iterator is_interval_heap_until (Iterator first, Iterator last,Compare compare);

/// @brief Checks whether the range is a valid interval heap.
/// @note Complexity: O(n)
/// @note Exception safety:No-throw
template <typename Iterator, typename Compare>
bool is_interval_heap (Iterator first, Iterator last, Compare compare);

} //  Namespace heap
} //  Namespace boost

//  Pull in the implementation.
#include "interval_heap.inl"

#endif // BOOST_HEAP_INTERVAL_HEAP_HPP_
