#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

class belady {
  const size_t cache_size;
  std::unordered_map <std::string, size_t> block_indexes;
  size_t complete = 1;
  size_t current = 0;

 public:
  belady(size_t cache_size)
      : cache_size(cache_size)
  {}
  void access(const std::string &a) {
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
        size_t number_equal = static_cast<size_t>(std::count_if(
            block_indexes.cbegin(), block_indexes.cend(),
            [temp](std::pair<std::string, size_t> p) {
              return temp==p.second;
            }));
        //std::cout << "    number_equal=" << number_equal << endl;
        count += number_equal;
        // This assert *will* fail if cache_size==1.
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
    std::cout << "complete=" << complete
              << " current=" << current
              << " {";
    // Print the block_indexes grouped by same index.
    std::map<size_t, std::set<std::string>> invert;
    for (const auto& [addr, index] : block_indexes) {
      invert[index].insert(addr);
    }
    size_t icount = 0;
    for (const auto& [index, addresses] : invert) {
      if (icount > 0) std::cout << " ";
      std::cout << "[" << index << "]={";
      size_t acount = 0;
      for (const std::string& address : addresses) {
        if (acount > 0) std::cout << ", ";
        std::cout << address;
        ++acount;
      }
      std::cout << "}";
      ++icount;
    }
    if (0) {
      size_t count = 0;
      for (const auto &pair : block_indexes) {
        if (count > 0) std::cout << " ";
        std::cout << pair.first << ":" << pair.second;
        count++;
      }
      std::cout << "}" << std::endl;
    }
    std::cout << std::endl;
  }
  size_t get_misses() const {
    return current;
  }
};

#if 0
struct tinfo {
  string address;
  size_t original_index;
  size_t next_access;
};
size_t opt(const vector<string> &trace, size_t cache_size) {
  vector<tinfo> infos;
  size_t i = 0;
  for (const string& a : trace) {
    infos.push_back({a, i++, SIZE_MAX});
  }
  stable_sort(infos.begin(), infos.end(),
              [](const tinfo &a, const tinfo &b) {
                return a.address < b.address;
              });
  for (size_t i = 0; i + 1 < infos.size(); ++i) {
    if (infos[i].address == infos[i+1].address) {
      infos[i].next_access = infos[i+1].original_index;
    }
  }
  sort(infos.begin(), infos.end(),
       [](const tinfo &a, const tinfo &b) {
         return a.original_index < b.original_index;
       });
  auto cmp_access = [](const tinfo &a, const tinfo &b) {
    return a.next_access > b.next_access;
  };
  priority_queue<tinfo, vector<tinfo>, decltype(cmp_access)> queue(cmp_access);
  for (const tinfo& ti : infos) {

  }
}
#endif

int main(int argc [[maybe_unused]], [[maybe_unused]] const char *arg[]) {
  const std::vector<std::string> trace =
      {"a", "b", "c", "d", "e", "b", "c", "f", "a", "b", "g", "d"};
  for (size_t c = 2; c < 8; ++c) {
    belady b(c);
    b.print();
    for (const std::string &a : trace) {
      b.access(a);
      std::cout << a << " ";
      b.print();
    }
    std::cout << "c=" << c << " misses=" << b.get_misses() << std::endl;
  }
  return 0;
}
