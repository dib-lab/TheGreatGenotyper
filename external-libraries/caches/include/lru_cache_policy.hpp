#ifndef LRU_CACHE_POLICY_HPP
#define LRU_CACHE_POLICY_HPP

#include "cache_policy.hpp"
#include <list>
#include <mutex>
#include <unordered_map>

namespace caches
{
template <typename Key> class LRUCachePolicy : public ICachePolicy<Key>
{
  public:
    using lru_iterator = typename std::list<Key>::iterator;
    using touch_guard = typename std::unique_lock<std::mutex>;

    LRUCachePolicy() = default;
    ~LRUCachePolicy() = default;

    LRUCachePolicy(const LRUCachePolicy &other)
        : lru_queue(other.lru_queue), key_finder(other.key_finder)
    {
    }

    void Insert(const Key &key) override
    {
        lru_queue.emplace_front(key);
        key_finder[key] = lru_queue.begin();
    }

    void Touch(const Key &key) override
    {
        auto it = key_finder[key];

        touch_guard lock{touch_op};

        // move the touched element at the beginning of the lru_queue
        lru_queue.splice(lru_queue.begin(), lru_queue, it);
    }

    void Erase(const Key &) override
    {
        // remove the least recently used element
        key_finder.erase(lru_queue.back());
        lru_queue.pop_back();
    }

    void Clear() override
    {
        key_finder.clear();
        lru_queue.clear();
    }

    // return a key of a displacement candidate
    const Key &ReplCandidate() const override
    {
        return lru_queue.back();
    }

  private:
    std::list<Key> lru_queue;
    std::unordered_map<Key, lru_iterator> key_finder;
    std::mutex touch_op;
};
} // namespace caches

#endif // LRU_CACHE_POLICY_HPP
