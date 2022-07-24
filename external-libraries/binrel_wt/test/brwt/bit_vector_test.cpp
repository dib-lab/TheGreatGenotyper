#include "doctest.h"
#include <brwt/bit_vector.h>

#include <algorithm> // min
#include <ostream>   // ostream
#include <string>    // string

// TODO(Diego): Test size() instead of length(). Remove length(). Test
// num_blocks(). Test data().

using brwt::bit_vector;
static_assert(bit_vector::bits_per_block >= 64, "");

namespace brwt {
static std::ostream& operator<<(std::ostream& os, const bit_vector& v) {
  for (auto i = v.length() - 1; i >= 0; --i) {
    os << (v.get(i) ? '1' : '0');
  }
  return os;
}

static bool operator==(const bit_vector& lhs, const bit_vector& rhs) noexcept {
  if (lhs.length() != rhs.length()) {
    return false;
  }
  for (bit_vector::size_type i = 0; i < lhs.length(); ++i) {
    if (lhs.get(i) != rhs.get(i)) {
      return false;
    }
  }
  return true;
}
} // end namespace brwt

TEST_SUITE("bit_vector");

TEST_CASE("bit_vector::bit_vector()") {
  CHECK(bit_vector() == bit_vector(0));

  static_assert(noexcept(bit_vector()), "");
}

TEST_CASE("bit_vector::bit_vector(size_type)") {
  CHECK(bit_vector(0) == bit_vector(0, 0));
  CHECK(bit_vector(32) == bit_vector(32, 0));
  CHECK(bit_vector(139) == bit_vector(139, 0));
}

TEST_CASE("bit_vector::bit_vector(size_type, block_type)") {
  SUBCASE("empty vector") {
    const bit_vector v(/*count=*/0, /*value=*/0b1101);
    CHECK(v.length() == 0);
  }
  SUBCASE("Initial value is zero") {
    const bit_vector v(/*count=*/130, /*value=*/0);
    CHECK(v.length() == 130);
    CHECK(v.get_chunk(/*pos=*/0, /*count=*/50) == 0);
    CHECK(v.get_chunk(/*pos=*/50, /*count=*/50) == 0);
    CHECK(v.get_chunk(/*pos=*/100, /*count=*/30) == 0);
  }
  SUBCASE("Initial value is truncated") {
    const bit_vector v(/*count=*/3, /*value=*/0b1101);
    CHECK(v.length() == 3);
    CHECK(v.get_chunk(/*pos=*/0, /*count=*/3) == 0b0101);
    CHECK(v.get_chunk(/*pos=*/0, /*count=*/2) == 0b0001);
  }
  SUBCASE("Extra bits are filled with 0") {
    SUBCASE("In small vectors") {
      const bit_vector v(/*count=*/10, /*value=*/0b1001101);
      CHECK(v.length() == 10);
      CHECK(v.get_chunk(/*pos=*/0, /*count=*/7) == 0b1001101);
      CHECK(v.get_chunk(/*pos=*/0, /*count=*/10) == 0b1001101);
    }
    SUBCASE("And big vectors") {
      const bit_vector v(/*count=*/873, /*value=*/0xFFFF'EEEE'AABB'CDCD);
      CHECK(v.length() == 873);
      CHECK(v.get_chunk(/*pos=*/0, /*count=*/64) == 0xFFFF'EEEE'AABB'CDCD);
      for (int pos = 64; pos < v.length(); pos += 64) {
        const auto count = std::min(64, static_cast<int>(v.length() - pos));
        CHECK(v.get_chunk(pos, count) == 0);
      }
    }
  }
}

TEST_CASE("bit_vector::bit_vector(const std::string&)") {
  CHECK(bit_vector("") == bit_vector());
  CHECK(bit_vector("100111") == bit_vector(6, 0b100111));
  CHECK(bit_vector("11111100111") == bit_vector(11, 0b11111100111));

  const bit_vector v(std::string(210, '1'));
  CHECK(v.length() == 210);
  CHECK(v.get_chunk(/*pos=*/0, /*count=*/64) == 0xFFFF'FFFF'FFFF'FFFF);
  CHECK(v.get_chunk(/*pos=*/64, /*count=*/64) == 0xFFFF'FFFF'FFFF'FFFF);
  CHECK(v.get_chunk(/*pos=*/128, /*count=*/64) == 0xFFFF'FFFF'FFFF'FFFF);
  CHECK(v.get_chunk(/*pos=*/192, /*count=*/18) == 0x3FFFF);
}

TEST_CASE("bit_vector::length()") {
  CHECK(bit_vector().length() == 0);
  CHECK(bit_vector(0).length() == 0);
  CHECK(bit_vector(32).length() == 32);
  CHECK(bit_vector(197).length() == 197);

  static_assert(noexcept(bit_vector().length()), "");
}

TEST_CASE("bit_vector::allocated_bytes()") {
  CHECK(bit_vector(0).allocated_bytes() == 0);
  CHECK(bit_vector(1).allocated_bytes() >= 8);
  CHECK(bit_vector(64).allocated_bytes() >= 8);
  CHECK(bit_vector(65).allocated_bytes() >= 16);
  CHECK(bit_vector(700).allocated_bytes() >= 88);

  static_assert(noexcept(bit_vector().allocated_bytes()), "");
}

TEST_CASE("bit_vector::get(size_type)") {
  const bit_vector v(100, 0xFF00'4FF4'FF11'33AA);
  CHECK(!v.get(99));
  CHECK(!v.get(0));
  CHECK(v.get(1));
  CHECK(!v.get(2));
  CHECK(v.get(3));
  CHECK(!v.get(55));
  CHECK(v.get(56));

  static_assert(noexcept(bit_vector().get(0)), "");
}

TEST_CASE("bit_vector::set(size_type, bool)") {
  bit_vector v(200, 0b1000'1111);
  REQUIRE(v.get(7));
  REQUIRE(!v.get(6));
  REQUIRE(v.get(1));

  v.set(6, true);
  v.set(7, false);
  v.set(0, false);
  v.set(199, true);

  CHECK(v.get(6));
  CHECK(!v.get(7));
  CHECK(!v.get(8));
  CHECK(!v.get(0));
  CHECK(!v.get(198));
  CHECK(v.get(199));

  static_assert(noexcept(bit_vector().set(0, false)), "");
}

TEST_CASE("bit_vector::get_chunk") {
  const bit_vector v(100, 0xFFFF'1EEE'1428'1192);
  CHECK(v.get_chunk(0, 16) == 0x1192);
  CHECK(v.get_chunk(16, 16) == 0x1428);
  CHECK(v.get_chunk(32, 16) == 0x1EEE);
  CHECK(v.get_chunk(48, 16) == 0xFFFF);
  CHECK(v.get_chunk(64, 16) == 0);

  CHECK(v.get_chunk(8, 16) == 0x2811);
  CHECK(v.get_chunk(24, 36) == 0xFFF'1EEE'14);

  static_assert(noexcept(bit_vector().get_chunk(0, 10)), "");
}

TEST_CASE("bit_vector::set_chunk") {
  bit_vector v(200, 0xFFFF'0033'5192'1001);
  REQUIRE(v.get_chunk(64, 64) == 0);
  REQUIRE(v.get_chunk(150, 50) == 0);
  REQUIRE(v.get_chunk(0, 32) == 0x5192'1001);
  REQUIRE(v.get_chunk(32, 32) == 0xFFFF'0033);

  v.set_chunk(16, 16, 0x1ACC);
  v.set_chunk(72, 32, 0x1489'1289);
  v.set_chunk(48, 5, 0b11010);
  v.set_chunk(53, 3, 0);
  v.set_chunk(60, 20, 0x1428'9);
  v.set_chunk(120, 64, 0x2214'1242'4412'6342);

  CHECK(v.get_chunk(0, 32) == 0x1ACC'1001);
  CHECK(v.get_chunk(0, 64) == 0x9F1A'0033'1ACC'1001);

  CHECK(v.get_chunk(48, 32) == 0x1428'9F1A);
  CHECK(v.get_chunk(64, 64) == 0x4200'0014'8912'1428);
  CHECK(v.get_chunk(120, 64) == 0x2214'1242'4412'6342);
  CHECK(v.get_chunk(128, 64) == 0x0022'1412'4244'1263);
  CHECK(v.get_chunk(160, 40) == 0x00'0022'1412);

  static_assert(noexcept(bit_vector().set_chunk(0, 8, 0xFF)), "");
}

TEST_CASE("bit_vector::get_block") {
  constexpr auto bpb = bit_vector::bits_per_block;

  bit_vector v(3 * bpb, 0xFF11'FF22);
  const auto& cv = v;

  CHECK(cv.get_block(0) == 0xFF11'FF22);
  CHECK(cv.get_block(1) == 0);
  CHECK(cv.get_block(2) == 0);

  v.set(bpb, true);
  v.set(2 * bpb + 2, true);
  CHECK(cv.get_block(0) == 0xFF11'FF22);
  CHECK(cv.get_block(1) == 1);
  CHECK(cv.get_block(2) == 4);

  static_assert(noexcept(bit_vector().get_block(0)), "");
}

TEST_CASE("bit_vector::set_block") {
  constexpr auto bpb = bit_vector::bits_per_block;
  bit_vector v(4 * bpb);
  v.set_block(0, 0xFF19'1984);
  v.set_block(2, 0x3492'4238);
  v.set_block(3, 0x4750'1434);

  CHECK(v.get_block(0) == 0xFF19'1984);
  CHECK(v.get_block(1) == 0);
  CHECK(v.get_block(2) == 0x3492'4238);
  CHECK(v.get_block(3) == 0x4750'1434);

  v.set_block(3, 0);
  v.set_block(1, 0x1489'1232);
  v.set_block(0, 0x4881);

  CHECK(v.get_block(0) == 0x4881);
  CHECK(v.get_block(1) == 0x1489'1232);
  CHECK(v.get_block(2) == 0x3492'4238);
  CHECK(v.get_block(3) == 0);

  static_assert(noexcept(bit_vector().set_block(0, 0xFFFF)), "");
}

TEST_CASE("bit_vector::serialize_and_load") {
  constexpr auto bpb = bit_vector::bits_per_block;
  bit_vector v(4 * bpb);
  v.set_block(0, 0xFF19'1984);
  v.set_block(2, 0x3492'4238);
  v.set_block(3, 0x4750'1434);
  const std::string file_name = "../test/vector_dump_test";
  std::ofstream outstream(file_name);
  v.serialize(outstream);
  outstream.close();

  std::ifstream instream(file_name);
  CHECK(v_loaded.load(instream) == true);
  instream.close();

  CHECK(v_loaded.get_block(0) == 0xFF19'1984);
  CHECK(v_loaded.get_block(1) == 0);
  CHECK(v_loaded.get_block(2) == 0x3492'4238);
  CHECK(v_loaded.get_block(3) == 0x4750'1434);
}

TEST_SUITE_END();
