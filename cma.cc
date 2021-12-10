// Critical Marker Algorithm from BilardiEkPa17
//
// This is the algorithm that runs in time O(nm) where n is the number of
// distinct addresses and m is the length of the trace.  This algorithm is
// mostly useful for explaining, since it's not a fast algorithm.

#include <cassert>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

template <class T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &v) {
  os << "{";
  bool first = true;
  for (const T &t : v) {
    if (!first) os << ", ";
    first = false;
    os << t;
  }
  return os << "}";
}

class Cma {
 public:
  size_t Access(const std::string &t) {
    const std::optional<size_t> depth = FindDepth(t);
    if (!depth.has_value()) {
      lru_stack_.insert(lru_stack_.begin(), t);
      IncrementAllCriticalMarkers();
      critical_markers_.push_back(0);
      RecomputeBeta();
      return lru_stack_.size();
    } else if (depth == 0) {
      // state remains the same
      return 1;
    } else {
      const size_t z = FindZ(*depth);
      const size_t depth_opt = FindDepthOpt(z);
      RotateLruStack(*depth);
      for (size_t i = 1; i < critical_markers_.size(); ++i) {
        if (critical_markers_[i] < depth) ++critical_markers_[i];
      }
      critical_markers_[depth_opt] = 1;
      std::cout << "Just before computing beta:" << std::endl << *this << std::endl;
      RecomputeBeta();
      return depth_opt;
    }
  }
  friend std::ostream& operator<<(std::ostream&os, const Cma &cma) {
    os << "L=" << cma.lru_stack_ << std::endl;
    os << "M=" << cma.critical_markers_ << std::endl;
    return os << "B=" << cma.beta_;
  }
 private:
  void RotateLruStack(size_t depth) {
    while (depth > 0) {
      std::swap(lru_stack_[depth-1], lru_stack_[depth]);
      depth--;
    }
  }
  std::optional<size_t> FindDepth(const std::string &t) {
    for (size_t i = 0; i < lru_stack_.size(); ++i) {
      if (lru_stack_[i] == t) {
        return i;
      }
    }
    return {};
  }
  void IncrementAllCriticalMarkers() {
    for (size_t &v : critical_markers_) ++v;
  }
  size_t FindZ(size_t depth) {
    size_t result = 0;
    for (size_t i = 1; i < depth; ++i) {
      if (beta_[i] == 0) result = i;
    }
    return result;
  }
  size_t FindDepthOpt(size_t z) {
    for (size_t i = 2; i < critical_markers_.size(); ++i) {
      if (critical_markers_[i] >= z) return i;
    }
    assert(0);
  }
  void RecomputeBeta() {
    beta_.resize(lru_stack_.size());
    for (size_t i = 0; i < beta_.size(); ++i) {
      beta_[i] = ComputeBetaI(i);
    }
  }
  size_t ComputeBetaI(size_t i) {
    size_t count = 0;
    for (size_t c =  1; c < critical_markers_.size(); ++c) {
      // Critical_markers_ is 1 less than what the paper expects, but so is i.
      if (critical_markers_[c] < i) ++count;
    }
    //std::cout << "C=" << critical_markers_ << std::endl;
    //std::cout << "count=" << count << " i=" << i << std::endl;
    // i is one less than one the paper expects, so we subtract by i rather than
    // i-1.
    assert(count >= i);
    return count - i;
  }

  // All these arrays are indexed from 0, but the paper indexes from 1.
  //
  // These vectors, when viewed as stacks have the "top" of the stack at element
  // 0 and the bottom at, e.g., lru_stack_.back();
  std::vector<std::string> lru_stack_;
  // Critical markers are 1 less than one the paper uses.  (Where the paper sets
  // the critical marker to 1, we set it to 0).
  std::vector<size_t> critical_markers_;
  std::vector<size_t> beta_;
};

int main() {
  const std::vector<std::string> trace =
      {"a", "b", "c", "d", "e", "b", "c", "f", "a", "b", "g", "d"};
  Cma cma;
  for (const std::string &a : trace) {
    std::cout << "Accessing " << a << std::endl;
    size_t c = cma.Access(a);
    std::cout << "x=" << a << " c=" << c << std::endl << cma << std::endl;
    std::cout << std::endl;
  }
}
