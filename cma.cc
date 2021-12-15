// Critical Marker Algorithm from BilardiEkPa17
//
// This is the algorithm that runs in time O(nm) where n is the number of
// distinct addresses and m is the length of the trace.  This algorithm is
// mostly useful for explaining, since it's not a fast algorithm.

#include <cassert>
#include <iostream>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "ost.h"

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

std::ostream& operator<<(std::ostream &os, const std::tuple<size_t, size_t> &v) {
  auto [ c, lru_depth ] = v;
  return os << "{" << c << ", " << lru_depth << "}";
}

std::ostream& operator<<(std::ostream& os, const std::optional<size_t> &v) {
  if (v == std::nullopt)
    return os << "nullopt";
  else
    return os << *v;
}

template <class T, class V>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<T, V> &map) {
  os << "{";
  bool first = true;
  for (const auto &[a, b] : map) {
    if (!first) os << ", ";
    os << "{" << a << ", " << b << "}";
  }
  return os << "}";
}

size_t ComputeBetaI(const std::vector<size_t> &critical_markers, size_t i) {
  size_t count = 0;
  for (size_t c =  1; c < critical_markers.size(); ++c) {
    // Critical_markers_ is 1 less than what the paper expects, but so is i.
    if (critical_markers[c] <= i) ++count;
  }
  //std::cout << "C=" << critical_markers_ << std::endl;
  //std::cout << "count=" << count << " i=" << i << std::endl;
  // i is one less than one the paper expects, so we subtract by i rather than
  // i-1.
  assert(count >= i);
  return count - i;
}

void ComputeBeta(const std::vector<size_t> &critical_markers,
                 /*out*/ std::vector<size_t> &beta) {
  beta.resize(critical_markers.size());
  for (size_t i = 0; i < beta.size(); ++i) {
    beta[i] = ComputeBetaI(critical_markers, i);
  }
}

class Cma {
 public:
  std::tuple</*OPT depth*/size_t, /*LRU depth*/size_t> Access(const std::string &t) {
    const std::optional<size_t> depth = FindDepth(t);
    most_recent_z_ = std::nullopt;
    most_recent_dopt_ = std::nullopt;
    if (!depth.has_value()) {
      lru_stack_.insert(lru_stack_.begin(), t);
      IncrementAllCriticalMarkers();
      critical_markers_.push_back(1);
      RecomputeBeta();
      return {lru_stack_.size(), lru_stack_.size()};
    } else if (depth == 0) {
      // state remains the same
      return {0, 0};
    } else {
      std::cout << "D=" << *depth << std::endl;
      const size_t z = FindZ(*depth);
      std::cout << "z=" << z << std::endl;
      const size_t depth_opt = FindDepthOpt(z);
      std::cout << "depth_opt=" << depth_opt << std::endl;
      RotateLruStack(*depth);
      for (size_t i = 1; i < critical_markers_.size(); ++i) {
        if (critical_markers_[i] < depth) ++critical_markers_[i];
      }
      critical_markers_[depth_opt] = 1;
      //std::cout << "Just before computing beta:" << std::endl << *this << std::endl;
      RecomputeBeta();
      //std::cout << "Just after  computing beta:" << std::endl << *this << std::endl;
      most_recent_z_ = z;
      most_recent_dopt_ = depth_opt;
      return {depth_opt, *depth};
    }
  }
  const std::vector<std::string> &GetL() const { return lru_stack_; }
  const std::vector<size_t> &GetM() const { return critical_markers_; }
  const std::vector<size_t> &GetBeta() const { return beta_; }
  friend std::ostream& operator<<(std::ostream&os, const Cma &cma) {
    os << "L=" << cma.lru_stack_ << std::endl;
    os << "M=" << cma.critical_markers_ << std::endl;
    return os << "B=" << cma.beta_;
  }
  std::optional<size_t> most_recent_z_;
  std::optional<size_t> most_recent_dopt_;
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
  size_t FindZ(size_t depth) const {
    size_t result = 0;
    for (size_t i = 0; i < depth; ++i) {
      if (beta_[i] == 0) result = i;
    }
    return result;
  }
  size_t FindDepthOpt(size_t z) const {
    for (size_t i = 1; i < critical_markers_.size(); ++i) {
      if (critical_markers_[i] > z) return i;
    }
    assert(0);
  }
  void RecomputeBeta() {
    ComputeBeta(critical_markers_, beta_);
  }

  // All these arrays are indexed from 0, but the paper indexes from 1.
  //
  // These vectors, when viewed as stacks have the "top" of the stack at element
  // 0 and the bottom at, e.g., lru_stack_.back();
  std::vector<std::string> lru_stack_;
  // Critical markers are are the same as what the paper says, but the indexes
  // are less.  For the paper's M(1) we store in critical_markers_[0].
  std::vector<size_t> critical_markers_;
  std::vector<size_t> beta_;
};


class TimeStamp {
 public:
  explicit TimeStamp(size_t t) :time_(t) {}
  TimeStamp() :TimeStamp(0) {}
  // Prefix increment
  TimeStamp& operator++() {
    ++time_;
    return *this;
  }
  friend bool operator<(const TimeStamp &a, const TimeStamp &b) {
    // Sort so that bigger timestamps are first.
    return b.time_ < a.time_;
  }
  friend std::ostream& operator<<(std::ostream&os, const TimeStamp &ts) {
    return os << ts.time_;
  }
 private:
  size_t time_;
};


class Csa {
 public:
  std::tuple</*OPT depth*/size_t, /*LRU depth*/size_t> Access(const std::string &t) {
    std::cout << std::endl << "CSA access " << t << std::endl;
    std::cout << *this << std::endl;
    const std::optional<size_t> depth = FindDepth(t);
    most_recent_z_ = std::nullopt;
    most_recent_dopt_ = std::nullopt;
    if (!depth.has_value()) {
      lru_stack_.InsertOrAssign(timestep_, t);
      stack_positions_[t] = timestep_;
      ++timestep_;
      IncrementAllCriticalMarkers();
      critical_markers_.push_back(1);
      RecomputeBeta();
      size_t size = lru_stack_.Size();
      return {size, size};
    } else if (depth == 0) {
      // state remains the same
      return {0, 0};
    } else {
      std::cout << "D=" << *depth << std::endl;
      const size_t z = FindZ(*depth);
      std::cout << "z=" << z << std::endl;
      const size_t depth_opt = FindDepthOpt(z);
      std::cout << "depth_opt=" << depth_opt << std::endl;
      RotateLruStack(*depth);
      ++timestep_;
      for (size_t i = 1; i < critical_markers_.size(); ++i) {
        if (critical_markers_[i] < depth) ++critical_markers_[i];
      }
      critical_markers_[depth_opt] = 1;
      //std::cout << "Just before computing beta:" << std::endl << *this << std::endl;
      RecomputeBeta();
      //std::cout << "Just after  computing beta:" << std::endl << *this << std::endl;
      most_recent_z_ = z;
      most_recent_dopt_ = depth_opt;
      return {depth_opt, *depth};
    }
  }
  const OrderStatisticTree<TimeStamp, std::string> &GetL() const { return lru_stack_; }
  const std::vector<size_t> &GetM() const { return critical_markers_; }
  const std::vector<size_t> &GetBeta() const { return beta_; }
  friend std::ostream& operator<<(std::ostream&os, const Csa &csa) {
    os << "L=" << csa.lru_stack_ << std::endl;
    os << "S=" << csa.stack_positions_ << std::endl;
    os << "M=" << csa.critical_markers_ << std::endl;
    return os << "B=" << csa.beta_;
  }
  std::optional<size_t> most_recent_z_;
  std::optional<size_t> most_recent_dopt_;
 private:
  void RotateLruStack(size_t depth) {
    std::optional<std::pair<TimeStamp, std::string>> pair = lru_stack_.Select(depth);
    assert(pair.has_value());
    stack_positions_[pair->second] = timestep_;
    lru_stack_.Erase(pair->first);
    lru_stack_.InsertOrAssign(timestep_, pair->second);
  }
  std::optional<size_t> FindDepth(const std::string &t) {
    auto it = stack_positions_.find(t);
    if (it == stack_positions_.end()) return {};
    auto [rank, pair] = lru_stack_.Rank(it->second);
    assert(pair.second == t);
    return rank;
  }
  void IncrementAllCriticalMarkers() {
    for (size_t &v : critical_markers_) ++v;
  }
  size_t FindZ(size_t depth) const {
    size_t result = 0;
    for (size_t i = 0; i < depth; ++i) {
      if (beta_[i] == 0) result = i;
    }
    return result;
  }
  size_t FindDepthOpt(size_t z) const {
    for (size_t i = 1; i < critical_markers_.size(); ++i) {
      if (critical_markers_[i] > z) return i;
    }
    assert(0);
  }
  void RecomputeBeta() {
    ComputeBeta(critical_markers_, beta_);
  }

  // For each string in lru_stack_, what is it's rank in lru_stack_.
  std::unordered_map<std::string, TimeStamp> stack_positions_;
  // All these arrays are indexed from 0, but the paper indexes from 1.
  //
  // This is sorted by the negative of the access time, so that the top of the
  // stack is first.
  OrderStatisticTree<TimeStamp, std::string> lru_stack_;
  // Critical markers are are the same as what the paper says, but the indexes
  // are less.  For the paper's M(1) we store in critical_markers_[0].
  std::vector<size_t> critical_markers_;
  std::vector<size_t> beta_;
  TimeStamp timestep_;
};

void CheckHelper(bool c, std::string_view expr) {
  if (!c) {
    std::cerr << "Failed check: " << expr << std::endl;
    exit(1);
  }
}
#define CHECK(x) CheckHelper(x, #x)

template <class T>
static void CheckEqHelper(const T &a, const T &b, std::string_view astring,
                          std::string_view bstring, std::string_view file,
                          int line) {
  if (a != b) {
    std::cerr << "At " << file << ":" << line << " Failed EQ check: "
              << astring << "=" << a << " not equal to " << bstring << "=" << b
              << std::endl;
    exit(1);
  }
}

#define CHECK_EQ(x, y) CheckEqHelper(x, y, #x, #y, __FILE__, __LINE__)

static void TestBetaFor(const std::vector<size_t> &m,
                        const std::vector<size_t> &expected_beta) {
  std::vector<size_t> beta;
  ComputeBeta(m, beta);
  if (beta != expected_beta) {
    std::cerr << "For M=" << m << " beta=" << beta << " expected beta="
              << expected_beta << std::endl;
    exit(1);
  }
}

static void TestBeta() {
  // These expected beta values are calculated by hand.
  TestBetaFor({1}, {0});

  TestBetaFor({0}, {0});
  TestBetaFor({1, 0}, {1, 0});
  TestBetaFor({2, 1, 0}, {1, 1, 0});
  TestBetaFor({3, 2, 1, 0}, {1, 1, 1, 0});
  TestBetaFor({4, 3, 2, 1, 0}, {1, 1, 1, 1, 0});
  TestBetaFor({4, 3, 2, 1, 1}, {0, 1, 1, 1, 0});
}

class CsaAndCma {
 public:
  std::tuple</*OPT depth*/size_t, /*LRU depth*/size_t> Access(const std::string &t) {
    auto cma_result = cma_.Access(t);
    auto csa_result = csa_.Access(t);
    std::cout << "cma result = " << cma_result << std::endl;
    std::cout << "csa result = " << csa_result << std::endl;
    assert(cma_result == csa_result);
    return csa_result;
  }
  const std::vector<std::string> &GetL() const {
    const std::vector<std::string> &cma_result = cma_.GetL();
    const OrderStatisticTree<TimeStamp, std::string> &csa_result = csa_.GetL();
    assert(cma_result.size() == csa_result.Size());
    for (size_t i = 0; i < cma_result.size(); ++i) {
      std::optional<std::pair<TimeStamp, std::string>> ith = csa_result.Select(i);
      assert(ith.has_value());
      assert(cma_result[i] == ith->second);
    }
    return cma_result;
  }
  const std::vector<size_t> &GetM() const { return cma_.GetM(); }
  const std::vector<size_t> &GetBeta() const { return cma_.GetBeta(); }
  std::optional<size_t> GetMostRecentZ() const { return cma_.most_recent_z_; }
  std::optional<size_t> GetMostRecentDopt() const { return cma_.most_recent_dopt_; }
  friend std::ostream& operator<<(std::ostream&os, const CsaAndCma &both) {
    return os << both.cma_;
  }

 private:
  Cma cma_;
  Csa csa_;
};

using TI = std::tuple<size_t, size_t>;
using VI = std::vector<size_t>;

static void TestABCDEB() {
  CsaAndCma cma;
  CHECK_EQ(cma.Access("a"), TI({1ul, 1ul}));
  CHECK_EQ(cma.GetL(), {"a"});
  CHECK_EQ(cma.GetM(), {1});
  CHECK_EQ(cma.GetBeta(), {0});

  CHECK_EQ(cma.Access("b"), TI({2ul, 2ul}));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"b", "a"}));
  CHECK_EQ(cma.GetM(), VI({2, 1}));
  CHECK_EQ(cma.GetBeta(), VI({0, 0}));

  CHECK_EQ(cma.Access("c"), TI({3ul, 3ul}));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"c", "b", "a"}));
  CHECK_EQ(cma.GetM(), VI({3, 2, 1}));
  CHECK_EQ(cma.GetBeta(), VI({0, 0, 0}));

  CHECK_EQ(cma.Access("d"), TI({4ul, 4ul}));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"d", "c", "b", "a"}));
  CHECK_EQ(cma.GetM(), VI({4, 3, 2, 1}));
  CHECK_EQ(cma.GetBeta(), VI({0, 0, 0, 0}));

  CHECK_EQ(cma.Access("e"), TI({5ul, 5ul}));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"e", "d", "c", "b", "a"}));
  CHECK_EQ(cma.GetM(), VI({5, 4, 3, 2, 1}));
  CHECK_EQ(cma.GetBeta(), VI({0, 0, 0, 0, 0}));

  std::cout << cma << std::endl;
  CHECK_EQ(cma.Access("b"), TI({1ul, 3ul}));
  // In the original paper, this would be 3, but since we are using 1-based
  // indexing it's two.
  CHECK_EQ(cma.GetMostRecentZ(), std::optional<size_t>(2ul));
  CHECK_EQ(cma.GetMostRecentDopt(), std::optional<size_t>(1ul));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"b", "e", "d", "c", "a"}));
  CHECK_EQ(cma.GetM(), VI({5, 1, 3, 3, 2}));
}

int main() {
  TestBeta();
  TestABCDEB();
  const std::vector<std::pair<std::string, size_t>> trace =
      {{"a", 1}, // OPT=a
       {"b", 2}, // OPT=b a
       {"c", 3}, // OPT=c b a  (c must be in cache 0, but is still before a)
       // OPT=b c a   (c is next accessed after b)
       {"d", 4}, // OPT=d b c a   (d is in cache0, but b beats c, and then c beats a)
       {"e", 5}, // OPT=e b c a d
       {"b", 1}, // OPT=b e c a d  (b moves up, and then e doesn't have to fight c)
       // Something goes wrong in my CMA implementation here.
       {"c", 2}, // OPT=c b e a d  (c moves up, b beats e, e doesn't have to fight a)
       {"f", 6}, // OPT=f b e a d c (f moves up, c moves out (but could have had e move out)
       {"a", 3}, // OPT=a b e f d c
       {"b", 1}, // OPT=b a e f d c
       {"g", 7}, // OPT=g a e f d c b
       {"d", 4}  // OPT=d a e f g c b
      };
  Cma cma;
  for (auto &[address, depth] : trace) {
    std::cout << "Accessing " << address << std::endl;
    auto [c, lru_depth] = cma.Access(address);
    std::cout << "a=" << address << " c=" << c << std::endl << cma << std::endl;
    std::cout << std::endl;
    assert(depth == SIZE_MAX || depth == c);
  }
}
