#include "doctest.h"
#include <brwt/bit_hacks.h>

#include <cstdint> // uint32_t, uint64_t

/// TODO(Diego): Test the count_leading_zeros overloads.
/// TODO(Diego): Test the count_trailing_zeros overloads.
/// TODO(Diego): Test used_bits

using uint = unsigned int;
using ulong = unsigned long;
using ullong = unsigned long long;
using std::uint32_t;
using std::uint64_t;

TEST_SUITE("bit_hacks");

TEST_CASE("pop_count") {
  using brwt::pop_count;

  SUBCASE("unsigned overload") {
    CHECK(pop_count(0x1111u) == 4);
    CHECK(pop_count(0x1212u) == 4);
    CHECK(pop_count(0x1F2Fu) == 10);
    CHECK(pop_count(0xFFFFu) == 16);
  }
  SUBCASE("unsigned long overload") {
    CHECK(pop_count(0x1111'1111ul) == 8);
    CHECK(pop_count(0x1214'2814ul) == 8);
    CHECK(pop_count(0x12F2'439Eul) == 15);
    CHECK(pop_count(0xFFFF'FFFFul) == 32);
  }
  SUBCASE("unsigned long long overload") {
    CHECK(pop_count(0x1111'1111'1111'1111ull) == 16);
    CHECK(pop_count(0x1248'1428'1144'2488ull) == 16);
    CHECK(pop_count(0x14C0'2491'A9B3'2390ull) == 23);
    CHECK(pop_count(0xFFFF'FFFF'FFFF'FFFFull) == 64);
  }

  // check constexpr
  static_assert(pop_count(0x1111u) == 4, "");
  static_assert(pop_count(0x12F2'439Eul) == 15, "");
  static_assert(pop_count(0x14C0'2491'A9B3'2390ull) == 23, "");
  // check noexcept
  static_assert(noexcept(pop_count(uint{})), "");
  static_assert(noexcept(pop_count(ulong{})), "");
  static_assert(noexcept(pop_count(ullong{})), "");
}

TEST_CASE("lsb_mask") {
  using brwt::lsb_mask;

  CHECK(lsb_mask<unsigned>(1) == 0b0000'0001);
  CHECK(lsb_mask<unsigned>(3) == 0b0000'0111);
  CHECK(lsb_mask<unsigned>(5) == 0b0001'1111);
  CHECK(lsb_mask<unsigned>(7) == 0b0111'1111);

  CHECK(lsb_mask<uint64_t>(62) == 0x3FFF'FFFF'FFFF'FFFF);
  CHECK(lsb_mask<uint64_t>(63) == 0x7FFF'FFFF'FFFF'FFFF);

  // check constexpr
  static_assert(lsb_mask<uint>(8) == 0xFF, "");
  static_assert(lsb_mask<ulong>(8) == 0xFF, "");
  static_assert(lsb_mask<ullong>(8) == 0xFF, "");

  // check noexcept
  static_assert(noexcept(lsb_mask<uint>(4)), "");
  static_assert(noexcept(lsb_mask<ulong>(4)), "");
  static_assert(noexcept(lsb_mask<ullong>(4)), "");
}

TEST_CASE("[rank_0][integer][full_range]") {
  using brwt::rank_0;

  static_assert(rank_0<uint32_t>(0xFFFF'FFFF) == 0, "");
  static_assert(rank_0<uint32_t>(0x0000'0000) == 32, "");
  static_assert(rank_0<uint32_t>(0x0101'1F0F) == 21, "");
  static_assert(rank_0<uint32_t>(0x0E21'9239) == 20, "");
  static_assert(rank_0<uint32_t>(0x0001'0000) == 31, "");
  static_assert(rank_0<uint32_t>(0x0F00'0000) == 28, "");

  static_assert(rank_0<uint64_t>(0x0000'0000'0000'0000) == 64, "");
  static_assert(rank_0<uint64_t>(0x00FF'0F0F'0F0F'FF72) == 28, "");
  static_assert(rank_0<uint64_t>(0x1211'1128'4281'1488) == 48, "");

  static_assert(noexcept(rank_0(uint{})), "");
  static_assert(noexcept(rank_0(ulong{})), "");
  static_assert(noexcept(rank_0(ullong{})), "");
}

TEST_CASE("[rank_1][integer][full_range]") {
  using brwt::rank_1;

  static_assert(rank_1<uint32_t>(0x0000'0000) == 0, "");
  static_assert(rank_1<uint32_t>(0xFFFF'FFFF) == 32, "");
  static_assert(rank_1<uint32_t>(0x0101'7ED7) == 14, "");
  static_assert(rank_1<uint32_t>(0xFFFF'1122) == 20, "");
  static_assert(rank_1<uint32_t>(0x0F01'090F) == 11, "");
  static_assert(rank_1<uint32_t>(0x0020'0001) == 2, "");

  static_assert(rank_1<uint64_t>(0xFFFF'FFFF'FFFF'FFFF) == 64, "");
  static_assert(rank_1<uint64_t>(0xFF00'EEFA'2313'5123) == 32, "");
  static_assert(rank_1<uint64_t>(0x4233'891A'1241'1213) == 21, "");

  static_assert(noexcept(rank_1(uint{})), "");
  static_assert(noexcept(rank_1(ulong{})), "");
  static_assert(noexcept(rank_1(ullong{})), "");
}

TEST_CASE("[rank_0][integer][end_position]") {
  using brwt::rank_0;

  static_assert(rank_0<uint32_t>(0b0000, 0) == 1, "");
  static_assert(rank_0<uint32_t>(0b0000, 1) == 2, "");
  static_assert(rank_0<uint32_t>(0b0000, 2) == 3, "");
  static_assert(rank_0<uint32_t>(0b0011, 0) == 0, "");
  static_assert(rank_0<uint32_t>(0b0011, 1) == 0, "");
  static_assert(rank_0<uint32_t>(0b0011, 2) == 1, "");
  static_assert(rank_0<uint32_t>(0b0011, 3) == 2, "");

  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 0) == 1, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 1) == 2, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 2) == 2, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 10) == 7, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 20) == 14, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 30) == 21, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 40) == 29, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 50) == 36, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 61) == 43, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 62) == 44, "");
  static_assert(rank_0<uint64_t>(0x1234'1241'4123'1514, 63) == 45, "");

  // 1234'1241'4123'1514
  // end 0001'0010'0011'0100 <-    (48..63: 11, total: 45)
  //     0001'0010'0100'0001 <-    (32..47: 12, total: 34)
  //     0100'0001'0010'0011 <-    (16..31: 11, total: 22)
  //     0001'0101'0001'0100 begin ( 0..15: 11, total: 11)

  static_assert(noexcept(rank_0(0u, 0)), "");
  static_assert(noexcept(rank_0(0ul, 0)), "");
  static_assert(noexcept(rank_0(0ull, 0)), "");
}

TEST_CASE("[rank_1][integer][end_position]") {
  using brwt::rank_1;

  static_assert(rank_1<uint32_t>(0b1111, 0) == 1, "");
  static_assert(rank_1<uint32_t>(0b1111, 1) == 2, "");
  static_assert(rank_1<uint32_t>(0b1111, 2) == 3, "");
  static_assert(rank_1<uint32_t>(0b1111, 3) == 4, "");
  static_assert(rank_1<uint32_t>(0b1010, 0) == 0, "");
  static_assert(rank_1<uint32_t>(0b1010, 1) == 1, "");
  static_assert(rank_1<uint32_t>(0b1010, 2) == 1, "");
  static_assert(rank_1<uint32_t>(0b1010, 3) == 2, "");

  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 0) == 1, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 1) == 2, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 2) == 2, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 10) == 4, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 20) == 9, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 30) == 12, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 40) == 15, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 50) == 17, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 61) == 20, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 62) == 20, "");
  static_assert(rank_1<uint64_t>(0x8940'1258'4123'5983, 63) == 21, "");

  // 8940'1258'4123'5983
  // end 1000'1001'0100'0000 <-    (48..63: 4, total: 21)
  //     0001'0010'0101'0100 <-    (32..47: 5, total: 17)
  //     0100'0001'0010'0011 <-    (16..31: 5, total: 12)
  //     0101'1001'1000'0011 begin ( 0..15: 7, total: 7)

  static_assert(noexcept(rank_1(0u, 0)), "");
  static_assert(noexcept(rank_1(0ul, 0)), "");
  static_assert(noexcept(rank_1(0ull, 0)), "");
}

TEST_CASE("is_power_of_two") {
  using brwt::is_power_of_two;

  static_assert(is_power_of_two(1), "");
  static_assert(is_power_of_two(2), "");
  static_assert(is_power_of_two(4), "");
  static_assert(is_power_of_two(8), "");
  static_assert(is_power_of_two(16), "");
  static_assert(is_power_of_two(32), "");
  static_assert(is_power_of_two(64), "");
  static_assert(is_power_of_two(128), "");
  static_assert(is_power_of_two(256), "");
  static_assert(is_power_of_two(512), "");

  static_assert(!is_power_of_two(3), "");
  static_assert(!is_power_of_two(5), "");
  static_assert(!is_power_of_two(6), "");
  static_assert(!is_power_of_two(7), "");
  static_assert(!is_power_of_two(9), "");
  static_assert(!is_power_of_two(10), "");
  static_assert(!is_power_of_two(15), "");
  static_assert(!is_power_of_two(31), "");
  static_assert(!is_power_of_two(500), "");
}

TEST_SUITE_END();
