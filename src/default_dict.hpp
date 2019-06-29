#ifndef SRC_DEFAULT_DICT_HPP_
#define SRC_DEFAULT_DICT_HPP_

#include <unordered_map>

// A special class for
//
// Inheriting from STL containers is frowned upon.
// However, private inheritance mitigates the problems, in particular it's not
// possible to have a pointer to the base class.
// https://stackoverflow.com/a/19355666/467327
// https://stackoverflow.com/a/2035044/467327
// I would have just had a polymorphic function to do things, but we need to
// store a default value.

template <
    class Key, class Value,
    class Comparator = typename std::unordered_map<Key, Value>::key_compare,
    class Allocator = typename std::unordered_map<Key, Value>::allocator_type>
class DefaultDict
    : private std::unordered_map<Key, Value, Comparator, Allocator> {
 private:
  const Value default_value_;

 public:
  explicit DefaultDict(Value default_value) : default_value_(default_value){};
  using std::unordered_map<Key, Value, Comparator, Allocator>::size;

  Value& at(const Key &key) {
    auto search = this->find(key);
    if (search == this->end()) {
      // TODO this seems like a problem. Since we are returning a reference we
      // could change default_value.
      return default_value_;
    }
    return std::unordered_map<Key, Value, Comparator, Allocator>::at(key);
  }

  void increment(const Key &key, const Value &value) {
    auto search = this->find(key);
    if (search == this->end()) {
      assert(this->insert(std::make_pair(key, value)).second);
    } else {
      search->second += value;
    }
  };

};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("DefaultDict") { auto d = DefaultDict<int, int>(); }

#endif  // DOCTEST_LIBRARY_INCLUDED


#endif  // SRC_DEFAULT_DICT_HPP_
