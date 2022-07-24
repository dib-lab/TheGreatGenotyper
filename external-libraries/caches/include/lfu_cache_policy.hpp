#ifndef LFU_CACHE_POLICY_HPP
#define LFU_CACHE_POLICY_HPP

#include "cache_policy.hpp"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <mutex>
#include <unordered_map>

namespace caches
{
template <typename Key> class LFUCachePolicy : public ICachePolicy<Key>
{
  public:
    using lfu_iterator = typename std::multimap<std::size_t, Key>::iterator;
    using touch_guard = typename std::unique_lock<std::mutex>;

    LFUCachePolicy() = default;
    ~LFUCachePolicy() override = default;

    LFUCachePolicy(const LFUCachePolicy &other)
        : frequency_storage(other.frequency_storage),
          lfu_storage(other.lfu_storage)
    {
    }

    void Insert(const Key &key) override
    {
        constexpr std::size_t INIT_VAL = 1;
        // all new value initialized with the frequency 1
        lfu_storage[key] = frequency_storage.emplace_hint(
            frequency_storage.cbegin(), INIT_VAL, key);
    }

    void Touch(const Key &key) override
    {
        touch_guard lock{touch_op};

        // get the previous frequency value of a key
        auto elem_for_update = lfu_storage[key];
        auto updated_elem =
            std::make_pair(elem_for_update->first + 1, elem_for_update->second);
        // update the previous value
        frequency_storage.erase(elem_for_update);
        lfu_storage[key] = frequency_storage.emplace_hint(
            frequency_storage.cend(), std::move(updated_elem));
    }

    void Erase(const Key &key) override
    {
        auto it = lfu_storage.find(key);
        assert(it != lfu_storage.end());

        frequency_storage.erase(it->second);
        lfu_storage.erase(it);
    }

    void Clear() override
    {
        frequency_storage.clear();
        lfu_storage.clear();
    }

    const Key &ReplCandidate() const override
    {
        // at the beginning of the frequency_storage we have the
        // least frequency used value
        return frequency_storage.cbegin()->second;
    }

  private:
    std::multimap<std::size_t, Key> frequency_storage;
    std::unordered_map<Key, lfu_iterator> lfu_storage;
    std::mutex touch_op;
};
} // namespace caches

#endif // LFU_CACHE_POLICY_HPP
