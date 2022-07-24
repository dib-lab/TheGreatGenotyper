Priority-Deque

==============



This project is my implementation of a priority deque (double-ended priority queue) as a C++ class template.

Documentation can be generated automatically by using Doxygen.

Goals:
* Interface similar to STL.
  * Status: Complete.
  * priority_deque implements the priority_queue interface.
  * Additional functions, for double-ended access and random-access among others, are in STL style.
* Interface in Boost style.
  * Status: Incomplete.
* Performance competitive with STL's std::priority_queue.
  * Status: Very close.
  * Benchmarks indicate ~20% slower performance in push, ~5% slower in pop. (GCC, -O3)
  * Threaded execution can be enabled for bulk-load and merging operations. Performance, as a ratio of time for priority_deque operation to time for std::priority_queue operation, follows:
    * 1.228 : 1 (1 CPU core)
    * 0.874 : 1 (2 CPU cores)
    * 0.545 : 1 (4 CPU cores)
* Provide the strongest exception-safety guarantee possible without compromising performance.
  * Status: Complete.
  * Most operations provide the strong (rollback) exception-safety guarantee if comparison throws.
  * Bulk-loading operations cannot provide this guarantee without severely compromising performance.
  * The random-access update operation can only provide this guarantee by copying an element or using move semantics.
* Provide thread safety, if possible.
  * Status: Mostly impractical.
  * Functions have no side-effects besides the obvious (instances of the deque may be altered), but a single deque should not be accessed by multiple threads simultaneously. Multiple reads are okay; simultaneous read-write and write-write induce race conditions.
  * Automatic locking via mutex crippled performance, and has been removed.