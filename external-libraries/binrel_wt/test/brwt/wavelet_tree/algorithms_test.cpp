#include "doctest.h"
#include <brwt/wavelet_tree.h>

using brwt::index_npos;
using brwt::index_type;
using brwt::int_vector;
using brwt::symbol_id;
using brwt::wavelet_tree;

static constexpr symbol_id operator"" _sym(const unsigned long long value) {
  return static_cast<symbol_id>(value);
}

TEST_SUITE("wavelet_tree/algorithms");

TEST_CASE("[select_first] small alphabet") {
  int_vector vec = {0, 2, 2, 1, 2, 3, 1, 3, 2, 1, 3, 0,
                    0, 1, 2, 0, 1, 0, 0, 0, 3, 3, 2, 1};

  const auto wt = wavelet_tree(vec);

  const auto select_first = [&](index_type start, symbol_id min,
                                symbol_id max) {
    return brwt::select_first(wt, start, brwt::between<symbol_id>{min, max});
  };

  // seq = 0221 2313 2130 0120 1000 3321

  CHECK(select_first(/*start=*/0, 0_sym, 1_sym) == 0);
  CHECK(select_first(/*start=*/0, 1_sym, 1_sym) == 3);

  CHECK(select_first(/*start=*/0, 0_sym, 2_sym) == 0);
  CHECK(select_first(/*start=*/0, 1_sym, 2_sym) == 1);
  CHECK(select_first(/*start=*/0, 2_sym, 2_sym) == 1);

  CHECK(select_first(/*start=*/0, 0_sym, 3_sym) == 0);
  CHECK(select_first(/*start=*/0, 1_sym, 3_sym) == 1);
  CHECK(select_first(/*start=*/0, 2_sym, 3_sym) == 1);
  CHECK(select_first(/*start=*/0, 3_sym, 3_sym) == 5);

  CHECK(select_first(/*start=*/11, 0_sym, 3_sym) == 11);
  CHECK(select_first(/*start=*/11, 1_sym, 3_sym) == 13);
  CHECK(select_first(/*start=*/11, 2_sym, 3_sym) == 14);
  CHECK(select_first(/*start=*/11, 3_sym, 3_sym) == 20);

  SUBCASE("If no element satisfies the conditions returns index_npos") {
    CHECK(select_first(/*start=*/20, 0_sym, 0_sym) == index_npos);
    CHECK(select_first(/*start=*/21, 0_sym, 0_sym) == index_npos);
    CHECK(select_first(/*start=*/22, 0_sym, 0_sym) == index_npos);
    CHECK(select_first(/*start=*/23, 0_sym, 0_sym) == index_npos);
  }
}

TEST_CASE("[select_first] medium alphabet") {
  int_vector vec = {0, 2, 2, 6, 5, 3, 6, 3, 2, 1, 3, 2,
                    0, 4, 5, 3, 7, 0, 0, 0, 3, 3, 2, 1};

  const auto wt = wavelet_tree(vec);

  const auto select_first = [&](index_type start, symbol_id min,
                                symbol_id max) {
    return brwt::select_first(wt, start, brwt::between<symbol_id>{min, max});
  };

  // seq = 0226 5363 2132 0453 7000 3321

  CHECK(select_first(/*start=*/0, 0_sym, 1_sym) == 0);
  CHECK(select_first(/*start=*/0, 1_sym, 1_sym) == 9);

  CHECK(select_first(/*start=*/0, 4_sym, 4_sym) == 13);

  CHECK(select_first(/*start=*/0, 5_sym, 7_sym) == 3);
  CHECK(select_first(/*start=*/0, 6_sym, 7_sym) == 3);
  CHECK(select_first(/*start=*/0, 7_sym, 7_sym) == 16);

  CHECK(select_first(/*start=*/9, 1_sym, 1_sym) == 9);
  CHECK(select_first(/*start=*/10, 1_sym, 1_sym) == 23);
  CHECK(select_first(/*start=*/11, 2_sym, 2_sym) == 11);
  CHECK(select_first(/*start=*/12, 2_sym, 2_sym) == 22);

  SUBCASE("If no element satisfies the conditions returns index_npos") {
    CHECK(select_first(/*start=*/17, 7_sym, 7_sym) == index_npos);
    CHECK(select_first(/*start=*/17, 6_sym, 7_sym) == index_npos);
    CHECK(select_first(/*start=*/17, 5_sym, 7_sym) == index_npos);
    CHECK(select_first(/*start=*/17, 4_sym, 7_sym) == index_npos);
    CHECK(select_first(/*start=*/17, 4_sym, 6_sym) == index_npos);
    CHECK(select_first(/*start=*/17, 4_sym, 5_sym) == index_npos);
    CHECK(select_first(/*start=*/17, 4_sym, 4_sym) == index_npos);

    CHECK(select_first(/*start=*/7, 6_sym, 6_sym) == index_npos);
    CHECK(select_first(/*start=*/8, 6_sym, 6_sym) == index_npos);
    CHECK(select_first(/*start=*/9, 6_sym, 6_sym) == index_npos);
    CHECK(select_first(/*start=*/10, 6_sym, 6_sym) == index_npos);
    CHECK(select_first(/*start=*/11, 6_sym, 6_sym) == index_npos);
    CHECK(select_first(/*start=*/12, 6_sym, 6_sym) == index_npos);
    CHECK(select_first(/*start=*/13, 6_sym, 6_sym) == index_npos);
  }
}

TEST_CASE("[select_first] few nodes") {
  int_vector vec = {0, 2, 2, 1};

  const auto wt = wavelet_tree(vec);

  const auto select_first = [&](index_type start, symbol_id min,
                                symbol_id max) {
    return brwt::select_first(wt, start, brwt::between<symbol_id>{min, max});
  };

  // seq = 0221

  CHECK(select_first(/*start=*/0, 0_sym, 1_sym) == 0);
  CHECK(select_first(/*start=*/1, 0_sym, 1_sym) == 3);
  CHECK(select_first(/*start=*/1, 0_sym, 2_sym) == 1);
  CHECK(select_first(/*start=*/1, 0_sym, 0_sym) == index_npos);

  CHECK(select_first(/*start=*/0, 1_sym, 1_sym) == 3);
  CHECK(select_first(/*start=*/0, 2_sym, 2_sym) == 1);
  CHECK(select_first(/*start=*/0, 3_sym, 3_sym) == index_npos);

  CHECK(select_first(/*start=*/3, 2_sym, 2_sym) == index_npos);
  CHECK(select_first(/*start=*/3, 3_sym, 3_sym) == index_npos);
}

TEST_CASE("[select_first] one node") {
  int_vector vec = {0, 0, 0, 1, 0, 0, 0, 0, 0};

  const auto wt = wavelet_tree(vec);

  const auto select_first = [&](index_type start, symbol_id min,
                                symbol_id max) {
    return brwt::select_first(wt, start, brwt::between<symbol_id>{min, max});
  };

  // seq = 0001 0000 0

  CHECK(select_first(/*start=*/0, 0_sym, 0_sym) == 0);
  CHECK(select_first(/*start=*/1, 0_sym, 0_sym) == 1);
  CHECK(select_first(/*start=*/0, 0_sym, 1_sym) == 0);
  CHECK(select_first(/*start=*/0, 1_sym, 1_sym) == 3);
  CHECK(select_first(/*start=*/5, 1_sym, 1_sym) == index_npos);
}

TEST_CASE("[select_first] one node, one symbol") {
  int_vector vec = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};

  const auto wt = wavelet_tree(vec);

  const auto select_first = [&](index_type start, symbol_id min,
                                symbol_id max) {
    return brwt::select_first(wt, start, brwt::between<symbol_id>{min, max});
  };

  CHECK(select_first(/*start=*/0, 0_sym, 0_sym) == 0);
  CHECK(select_first(/*start=*/0, 0_sym, 1_sym) == 0);
  CHECK(select_first(/*start=*/0, 1_sym, 1_sym) == index_npos);
  CHECK(select_first(/*start=*/4, 1_sym, 1_sym) == index_npos);
  CHECK(select_first(/*start=*/8, 1_sym, 1_sym) == index_npos);
}

TEST_SUITE_END();
