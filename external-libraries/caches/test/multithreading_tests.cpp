#include "cache.hpp"
#include "cache_policy.hpp"
#include "fifo_cache_policy.hpp"
#include "lfu_cache_policy.hpp"
#include "lru_cache_policy.hpp"

#include <thread>

#include <gtest/gtest.h>

template <typename Key, typename Value>
using lfu_cache_t =
    typename caches::fixed_sized_cache<Key, Value, caches::LFUCachePolicy<Key>>;

template <typename Key, typename Value>
using lru_cache_t =
    typename caches::fixed_sized_cache<Key, Value, caches::LRUCachePolicy<Key>>;

template <typename Key, typename Value>
using fifo_cache_t =
    typename caches::fixed_sized_cache<Key, Value,
                                       caches::FIFOCachePolicy<Key>>;

template <typename Key, typename Value>
using nopolicy_cache_t = typename caches::fixed_sized_cache<Key, Value>;

template <class Policy> class CacheTest : public ::testing::Test
{
};
typedef ::testing::Types<lfu_cache_t<int, size_t>, lru_cache_t<int, size_t>,
                         fifo_cache_t<int, size_t>,
                         nopolicy_cache_t<int, size_t>>
    CacheTypes;

TYPED_TEST_SUITE(CacheTest, CacheTypes);

TYPED_TEST(CacheTest, Multithreaded)
{
    TypeParam cache(90);

    std::vector<std::thread> threads;
    for (size_t i = 0; i < 4; ++i)
    {
        threads.emplace_back([&cache]() {
            for (size_t j = 0; j < 20000; ++j)
            {
                if (auto fetch = cache.TryGet(j % 100))
                {
                    EXPECT_EQ(j % 100, *fetch);
                }
                else
                {
                    cache.Put(j % 100, j % 100);
                }
            }
        });
    }

    for (auto &thread : threads)
    {
        thread.join();
    }
}
