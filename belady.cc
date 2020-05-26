#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class belady {
  const size_t cache_size;
  map <string, size_t> block_indexes;
  size_t complete = 1;
  size_t current = 0;

 public:
  belady(size_t cache_size)
      : cache_size(cache_size)
  {}
  void access(const string &a) {
    if (block_indexes[a] < complete) {
      current++;
      block_indexes[a] = current;
    } else if (block_indexes[a] == current) {
      // Nothing */
    } else {
      block_indexes[a] = current;
      size_t temp = current;
      size_t count = 0;
      while (1) {
        size_t number_equal = count_if(block_indexes.cbegin(),
                                       block_indexes.cend(),
                                       [temp](pair<string, size_t> p) {
                                         return temp==p.second;
                                       });
        //cout << "    number_equal=" << number_equal << endl;
        count += number_equal;
        assert(count <= cache_size);
        if (count == cache_size || temp == complete) {
          complete = temp;
          return;
        } else if (count < cache_size) {
          count--;
          temp--;
        } else {
          assert(0);
        }
      }
    }
  }
  void print() const {
    cout << "complete=" << complete
         << " current=" << current
         << " {";
    size_t count = 0;
    for (const auto &pair : block_indexes) {
      if (count > 0) cout << " ";
      cout << pair.first << ":" << pair.second;
      count++;
    }
    cout << "}" << endl;
  }
  size_t get_misses() const {
    return current;
  }
};

int main(int argc [[maybe_unused]], const char *arg[] [[maybe_unused]]) {
  vector<string> trace = {"a", "b", "c", "d", "e", "b", "c", "f", "a", "b", "g", "d"};
  struct belady b(3);
  b.print();
  for (const string &a : trace) {
    b.access(a);
    cout << a << " ";
    b.print();
  }
  cout << "misses=" << b.get_misses() << endl;
  return 0;
}
