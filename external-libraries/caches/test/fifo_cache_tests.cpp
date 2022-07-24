#include "cache.hpp"
#include "fifo_cache_policy.hpp"
#include <gtest/gtest.h>

template <typename Key, typename Value>
using fifo_cache_t =
    typename caches::fixed_sized_cache<Key, Value,
                                       caches::FIFOCachePolicy<Key>>;

TEST(FIFOCache, Simple_Test)
{
    fifo_cache_t<int, int> fc(2);

    fc.Put(1, 10);
    fc.Put(2, 20);

    EXPECT_EQ(fc.Size(), 2u);
    EXPECT_EQ(*fc.TryGet(1), 10);
    EXPECT_EQ(*fc.TryGet(2), 20);

    fc.Put(1, 30);
    EXPECT_EQ(fc.Size(), 2u);
    EXPECT_EQ(*fc.TryGet(1), 30);

    fc.Put(3, 30);
    EXPECT_FALSE(fc.TryGet(1));
    EXPECT_EQ(*fc.TryGet(2), 20);
    EXPECT_EQ(*fc.TryGet(3), 30);
}

TEST(FIFOCache, Missing_Value)
{
    fifo_cache_t<int, int> fc(2);

    fc.Put(1, 10);

    EXPECT_EQ(fc.Size(), 1u);
    EXPECT_EQ(*fc.TryGet(1), 10);
    EXPECT_FALSE(fc.TryGet(2));
    EXPECT_THROW(fc.Get(2), std::range_error);
}

TEST(FIFOCache, Sequence_Test)
{
    constexpr unsigned int TEST_SIZE = 10;
    fifo_cache_t<std::string, int> fc(TEST_SIZE);

    for (size_t i = 0; i < TEST_SIZE; ++i)
    {
        fc.Put(std::to_string('0' + i), i);
    }

    EXPECT_EQ(fc.Size(), TEST_SIZE);

    for (ssize_t i = 0; i < TEST_SIZE; ++i)
    {
        EXPECT_EQ(*fc.TryGet(std::to_string('0' + i)), i);
    }

    // replace a half
    for (size_t i = 0; i < TEST_SIZE / 2; ++i)
    {
        fc.Put(std::to_string('a' + i), i);
    }

    EXPECT_EQ(fc.Size(), TEST_SIZE);

    for (size_t i = 0; i < TEST_SIZE / 2; ++i)
    {
        EXPECT_FALSE(fc.TryGet(std::to_string('0' + i)));
        EXPECT_THROW(fc.Get(std::to_string('0' + i)), std::range_error);
    }

    for (ssize_t i = 0; i < TEST_SIZE / 2; ++i)
    {
        EXPECT_EQ(*fc.TryGet(std::to_string('a' + i)), i);
    }

    for (ssize_t i = TEST_SIZE / 2; i < TEST_SIZE; ++i)
    {
        EXPECT_EQ(*fc.TryGet(std::to_string('0' + i)), i);
    }
}
