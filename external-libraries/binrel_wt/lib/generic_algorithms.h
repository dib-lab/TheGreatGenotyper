#ifndef BINREL_WT_LIB_GENERIC_ALGORITHMS_H // NOLINT
#define BINREL_WT_LIB_GENERIC_ALGORITHMS_H

namespace brwt {

/// \brief Binary searches for the first integer that does not satisfy the
/// predicate.
///
/// \returns The first integer that does not satisfy the predicate, or \p b if
/// no such integer is found.
///
/// It is guaranteed that the algorithm never calls <tt>pred(b)</tt>.
///
template <typename T, typename Pred>
T int_binary_search(T a, T b, Pred pred) {
  while (a != b) {
    const T mid = a + (b - a) / 2;
    if (pred(mid)) {
      a = mid + 1;
    } else {
      b = mid;
    }
  }
  return a;
}

} // end namespace brwt

#endif // Header guard
