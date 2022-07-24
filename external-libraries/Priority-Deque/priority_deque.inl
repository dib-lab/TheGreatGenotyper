#include "priority_deque.hpp"

namespace boost {
namespace container {
//-------------------------------Constructors----------------------------------|
//----------------------------Default Constructor------------------------------|
template <typename T, typename S, typename C>
priority_deque<T, S, C>::priority_deque (const C& comp, const S& seq)
  : sequence_(seq), compare_(comp)
{
  heap::make_interval_heap(sequence_.begin(), sequence_.end(), compare_);
}
#if (__cplusplus >= 201103L)
template <typename T, typename S, typename C>
priority_deque<T, S, C>::priority_deque (const C& comp, S&& seq)
  : sequence_(std::move(seq)), compare_(comp)
{
  heap::make_interval_heap(sequence_.begin(), sequence_.end(), compare_);
}
#endif

//---------------------------Create from Iterators-----------------------------|
template <typename T, typename S, typename C>
template <typename InputIterator>
priority_deque<T, S, C>::priority_deque (InputIterator first,InputIterator last,
                                         const C& comp, const S& seq)
: sequence_(seq), compare_(comp)
{
  sequence_.insert(sequence_.end(), first, last);
  try {
    heap::make_interval_heap(sequence_.begin(), sequence_.end(), compare_);
  } catch (...) {
    sequence_.erase(sequence_.end() - std::distance(first, last),
                    sequence_.end());
    throw;  //  Re-throw the current exception.
  }
}
#if (__cplusplus >= 201103L)
template <typename T, typename S, typename C>
template <typename InputIterator>
priority_deque<T, S, C>::priority_deque (InputIterator first,InputIterator last,
                                         const C& comp, S&& seq)
: sequence_(std::move(seq)), compare_(comp)
{
  sequence_.insert(sequence_.end(), first, last);
  try {
    heap::make_interval_heap(sequence_.begin(), sequence_.end(), compare_);
  } catch (...) {
    sequence_.erase(sequence_.end() - std::distance(first, last),
                    sequence_.end());
    throw;  //  Re-throw the current exception.
  }
}
#endif

//-----------------------------Restricted Access-------------------------------|
//------------------------------Insert / Emplace-------------------------------|
template <typename T, typename Sequence, typename Compare>
void priority_deque<T, Sequence, Compare>::push (const value_type& value) {
  sequence_.push_back(std::move(value));
  try {
    heap::push_interval_heap(sequence_.begin(), sequence_.end(), compare_);
  } catch (...) {
    sequence_.pop_back();
    throw;  //  Re-throw the current exception.
  }
}
#if (__cplusplus >= 201103L)
template <typename T, typename Sequence, typename Compare>
void priority_deque<T, Sequence, Compare>::push (value_type&& value) {
  sequence_.push_back(std::move(value));
  try {
    heap::push_interval_heap(sequence_.begin(), sequence_.end(), compare_);
  } catch (...) {
    sequence_.pop_back();
    throw;  //  Re-throw the current exception.
  }
}

template <typename T, typename Sequence, typename Compare>
template<typename... Args>
void priority_deque<T, Sequence, Compare>::emplace (Args&&... args) {
  sequence_.emplace_back(std::forward<Args>(args)...);
  try {
    heap::push_interval_heap(sequence_.begin(), sequence_.end(), compare_);
  } catch (...) {
    sequence_.pop_back();
    throw;  //  Re-throw the current exception.
  }
}
#endif

//---------------------------Observe Maximum/Minimum---------------------------|
template <typename T, typename Sequence, typename Compare>
inline typename priority_deque<T, Sequence, Compare>::const_reference
  priority_deque<T, Sequence, Compare>::maximum  (void) const
{
  BOOST_CONTAINER_PRIORITY_DEQUE_ASSERT(!empty(),
    "Empty priority deque has no maximal element. Reference undefined.");
  const_iterator it = sequence_.begin() + 1;
  return (it == sequence_.end()) ? sequence_.front() : *it;
}

template <typename T, typename Sequence, typename Compare>
typename priority_deque<T, Sequence, Compare>::const_reference
  priority_deque<T, Sequence, Compare>::minimum  (void) const
{
  BOOST_CONTAINER_PRIORITY_DEQUE_ASSERT(!empty(),
    "Empty priority deque has no minimal element. Reference undefined.");
  return sequence_.front();
}
//---------------------------Remove Maximum/Minimum----------------------------|
template <typename T, typename Sequence, typename Compare>
void priority_deque<T, Sequence, Compare>::pop_minimum (void) {
  BOOST_CONTAINER_PRIORITY_DEQUE_ASSERT(!empty(),
    "Empty priority deque has no maximal element. Removal impossible.");
  heap::pop_interval_heap_min(sequence_.begin(), sequence_.end(), compare_);
  try {
    sequence_.pop_back();
  } catch (...) {
//  If pop_back has a strong guarantee, this will restore the heap.
    heap::push_interval_heap(sequence_.begin(), sequence_.end(), compare_);
    throw;  //  Re-throw the current exception.
  }
}

template <typename T, typename Sequence, typename Compare>
void priority_deque<T, Sequence, Compare>::pop_maximum (void) {
  BOOST_CONTAINER_PRIORITY_DEQUE_ASSERT(!empty(),
    "Empty priority deque has no minimal element. Removal undefined.");
  heap::pop_interval_heap_max(sequence_.begin(), sequence_.end(), compare_);
  try {
    sequence_.pop_back();
  } catch (...) {
//  If pop_back has a strong guarantee, this will restore the heap.
    heap::push_interval_heap(sequence_.begin(), sequence_.end(), compare_);
    throw;  //  Re-throw the current exception.
  }
}

//--------------------------Whole-Deque Operations-----------------------------|
template <typename T, typename S, typename C>
void priority_deque<T, S, C>::clear (void) {
  sequence_.clear();
}
//-----------------------------------Merge-------------------------------------|
template <typename T, typename S, typename C>
template <typename InputIterator>
void priority_deque<T, S, C>::merge (InputIterator first, InputIterator last) {
  sequence_.insert(sequence_.end(), first, last);
  try {
    heap::make_interval_heap(sequence_.begin(), sequence_.end(), compare_);
  } catch (...) {
    sequence_.erase(sequence_.end() - std::distance(first, last),
                    sequence_.end());
    throw;  //  Re-throw the current exception.
  }
}

//----------------------------Swap Specialization------------------------------|
template <typename T, typename S, typename C>
inline void priority_deque<T, S, C>::swap (priority_deque<T, S, C>& other) {
  using std::swap;

  swap(compare_, other.compare_);
  sequence_.swap(other.sequence_);
}

template <typename T, typename S, typename C>
inline void swap (priority_deque<T, S, C>& deque1,
                  priority_deque<T, S, C>& deque2)
{
  deque1.swap(deque2);
}

//---------------------------Random-Access Mutators----------------------------|
template <typename T, typename S, typename C>
void priority_deque<T, S, C>::update (const_iterator random_it,
                                      const value_type& value)
{
  const difference_type ind = random_it - begin();
  BOOST_CONTAINER_PRIORITY_DEQUE_ASSERT((0 <= ind) &&
    (ind < end() - begin()), "Iterator out of bounds; can't set element.");
//  Providing the strong guarantee would require saving a copy.
  *(sequence_.begin() + ind) = value;
  heap::update_interval_heap(sequence_.begin(),sequence_.end(),ind,compare_);
}
#if (__cplusplus >= 201103L)
//  Move-enabled function provides strong exception-safety guarantee.
template <typename T, typename S, typename C>
void priority_deque<T, S, C>::update (const_iterator random_it,
                                      value_type&& value)
{
  using namespace std;

  const difference_type ind = random_it - begin();
  BOOST_CONTAINER_PRIORITY_DEQUE_ASSERT((0 <= ind) &&
    (ind < end() - begin()), "Iterator out of bounds; can't set element.");
//  Providing the strong guarantee would require saving a copy.
  swap(*(sequence_.begin() + ind), value);
  try {
    heap::update_interval_heap(sequence_.begin(),sequence_.end(),ind,compare_);
  } catch (...) {
//  Restore the original value.
    swap(*(sequence_.begin() + ind), value);
    throw;  //  Re-throw current exception.
  }
}
#endif

template <typename T, typename Sequence, typename Compare>
void priority_deque<T, Sequence, Compare>::erase (const_iterator random_it) {
  const difference_type index = random_it - begin();
  BOOST_CONTAINER_PRIORITY_DEQUE_ASSERT((0 <= index) &&
    (index < end() - begin()), "Iterator out of bounds; can't erase element.");
//  Move the element pointed to by random_it to the end of the sequence.
  heap::pop_interval_heap(sequence_.begin(), sequence_.end(), index,compare_);
  try {
//  Remove the (moved) element.
    sequence_.pop_back();
  } catch (...) {
//  If pop_back failed, restore the heap property.
    heap::push_interval_heap(sequence_.begin(), sequence_.end(), compare_);
    throw;  //  Re-throw the current exception.
  }
}

//----------------------- General Container Information -----------------------|
template <typename T, typename Sequence, typename Compare>
bool priority_deque<T, Sequence, Compare>::empty (void) const
{
  return sequence_.empty();
}

template <typename T, typename Sequence, typename Compare>
typename priority_deque<T, Sequence, Compare>::size_type
  priority_deque<T, Sequence, Compare>::size (void) const
{
  return sequence_.size();
}

template <typename T, typename Sequence, typename Compare>
typename priority_deque<T, Sequence, Compare>::size_type
  priority_deque<T, Sequence, Compare>::max_size (void) const
{
  return sequence_.max_size();
}

} //  Namespace boost::container
} //  Namespace boost
