#ifndef ORDERED_MAP_HPP_
#define ORDERED_MAP_HPP_

#include <map>
#include <vector>

template <typename Key, typename Value> class orderedMap {
public:
  orderedMap() = default;
  ~orderedMap() = default;

  void insert(const Key &key, const Value &value)
  {
    if (map_.find(key) == map_.end()) {
      map_[key] = Value{};
      keys_.push_back(key);
    }
    map_.insert(key, value);
  }

  Value &operator[](const Key &key)
  {
    if (map_.find(key) == map_.end()) {
      map_[key] = Value{};
      keys_.push_back(key);
    }
    return map_[key];
  }

  const Value &operator[](const Key &key) const { return map_.at(key); }

  const std::vector<Key> &keys() const { return keys_; }

  bool contains(const Key &key) const { return map_.find(key) != map_.end(); }

  auto find(const Key &key) const { return map_.find(key); }

  auto find(const Key &key) { return map_.find(key); }

  auto begin() const { return map_.begin(); }

  auto end() const { return map_.end(); }

  auto begin() { return map_.begin(); }

  auto end() { return map_.end(); }

  auto size() const { return map_.size(); }

  void clear()
  {
    map_.clear();
    keys_.clear();
  }

private:
  std::map<Key, Value> map_{};
  std::vector<Key> keys_{};
};

#endif