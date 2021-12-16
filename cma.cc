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


class SizetComparesBackward {
 public:
  explicit SizetComparesBackward(size_t t) :v_(t) {}
  SizetComparesBackward() :SizetComparesBackward(0) {}
  // Prefix increment
  SizetComparesBackward& operator++() {
    ++v_;
    return *this;
  }
  // Postfix increment
  SizetComparesBackward operator++(int) {
    SizetComparesBackward temp = *this;
    ++*this;
    return temp;
  }
  friend bool operator<(const SizetComparesBackward &a, const SizetComparesBackward &b) {
    // Sort so that bigger timestamps are first.
    return b.v_ < a.v_;
  }
  friend std::ostream& operator<<(std::ostream&os, const SizetComparesBackward &ts) {
    return os << ts.v_;
  }
 private:
  size_t v_;
};

using TimeStamp = SizetComparesBackward;

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
      size_t size = lru_stack_.Size();
      gamma_.InsertOrAssign(timestep_, size);
      gamma_inverted_.insert({size, timestep_});
      stack_positions_[t] = timestep_;
      ++timestep_;
      IncrementAllCriticalMarkers();
      critical_markers_.push_back(1);
      RecomputeBeta();
      beta_diff_.InsertOrAssign(beta_counter_++, 0);
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
      {
        auto it = gamma_inverted_.find(depth_opt);
        assert(it != gamma_inverted_.end());
        TimeStamp gamma_at = it->second;
        it->second = timestep_;  // set gamma_inverted[*depth] = timestep_
        gamma_.Erase(gamma_at);
        gamma_.InsertOrAssign(timestep_, depth_opt);
      }
      ++timestep_;
      for (size_t i = 1; i < critical_markers_.size(); ++i) {
        if (critical_markers_[i] < depth) ++critical_markers_[i];
      }
      critical_markers_[depth_opt] = 1;
      std::cout << "Just before computing beta:" << std::endl << *this << std::endl;
      RecomputeBeta();
      std::cout << "Just after  computing beta:" << std::endl << *this << std::endl;
      std::cout << "depth_opt=" << depth_opt << std::endl;
      std::cout << "beta_diff_=" << beta_diff_ << std::endl;
      {
        const auto [k, v] = beta_diff_.Select(*depth);
        beta_diff_.InsertOrAssign(k, v+1);
        std::cout << "+1 beta_diff_=" << beta_diff_ << std::endl;
      }
      {
        size_t M_star = critical_markers_[depth_opt - 1];
        std::cout << "M_star=" << M_star << std::endl;
        const auto [k, v] = beta_diff_.Select(M_star - 1);
        std::cout << "select found k=" << k << std::endl;
        beta_diff_.InsertOrAssign(k, v-1);
        std::cout << "-1 beta_diff_=" << beta_diff_ << std::endl;
      }
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
    os << "G=" << csa.gamma_ << std::endl;
    os << "Ginv=" << csa.gamma_inverted_ << std::endl;
    os << "S=" << csa.stack_positions_ << std::endl;
    os << "M=" << csa.critical_markers_ << std::endl;
    return os << "B=" << csa.beta_;
  }
  void Validate() const {
    for (size_t i = 0; i < beta_.size(); ++i) {
      if (beta_[i] != static_cast<size_t>(beta_diff_.SelectPrefix(i))) {
        std::cout << "just before assertion failed:" << std::endl << *this << std::endl;
        std::cout << "beta_[" << i << "]=" << beta_[i] << std::endl;
        std::cout << "beta_diff_ prefix (" << i << ")=" << beta_diff_.SelectPrefix(i) << std::endl;
        assert(beta_[i] == static_cast<size_t>(beta_diff_.SelectPrefix(i)));
      }
    }
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
    auto [rank, pairptr] = lru_stack_.Rank(it->second);
    assert(pairptr);
    assert(pairptr->second == t);
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

  // In BilardiEkPa17 the critical stack algorithm is described two ways: as a
  // set of stacks and as a set of trees.
  //
  //   A is the address trace.  (The paper refers to a tree T^{A}, but all
  //     that's needed is a hash table.)
  //
  //   Lambda is the LRU stack (the tree version is called T^\Lambda in the
  //     paper).
  //
  //   Gamma is the sequence of capacities in the order in which they became
  //     most recently critical.  For example, suppose we have a trace that the
  //     last 4 accesses with critical stack distances of 3, 7, 3, 5.  That
  //     means that a cache of size 3 has a hit (but 4 would miss), then 7 has a
  //     hit (but 8 would miss), then 3 has a hit (but 4 would miss) then 5 has
  //     a hit (but 6 would miss.)  Then the first 3 elements of Gamma would be
  //     {5,3,7}.
  //
  //   Beta is the *imbalance function* (the tree version is called T^\Beta).
  //     The paper essentially describes a prefix sum tree in which Beta(i) is
  //     calculated by taking the sum from j=0 to i of all the values in
  //     T^\Beta.
  //
  // In our implementation:
  //
  //   stack_positions_ is approximately A from the paper.  stack_positions_
  //     tells us, for each address, the rank of the address in the lru stack.
  //
  //   lru_stack_ is the stack of addresses, corresponding.  We need to be able
  //     to erase the ith element from the top of the stack and insert an
  //     element at the top of the stack. We use an order-statistic tree to
  //     implement the stack.
  //
  //   gamma_ is the stack of critical distances.  It has the same operations
  //     needed on lru_stack_ and also uses an order-statistic tree.
  //
  //   beta_diff_ is a prefix sum tree.  The tree is indexed by a counter
  //     (beta_counter_) which increases but compares backwards so that larger
  //     counters sort earlier in the tree.  The leaves of the trees contained
  //     signed integers.  Those integers are kind of the derivative of beta_,
  //     that is, prefix_sum(beta_diff_, i).sum == beta_[i+1] = beta_[i].  (Keep
  //     in mind that we don't represent beta_ explicitly).  But we need a
  //     little bit more from beta_.  We need to be able to compute the minimum
  //     from the paper's equation (2):
  //
  //     Z(j) = max { i > j : Beta(i) == 0 }             (Equation 2)
  //
  //     We can compute Z(j) using a prefix sum in the beta_diff_ tree.  To do
  //     this the beta_diff_ prefix sums contain the following components:
  //
  //       sum:  The sum of the beta_diff_ leaves in the subtree.
  //
  //       relative_depth: If you start at 0 at the left edge of the subtree,
  //         how deep does the partial sum get?
  //
  //       distance_from_min: How many values since the rightmost instance of
  //         the deepest point.
  //
  //     When combining two prefix sums a and b we have
  //
  //       combined.sum = a.sum + b.sum
  //
  //       combined.relative_depth =
  //           min(a.relative_depth, a.sum + b.relative_depth)
  //
  //       combined.distance_from_min =
  //           (a.relative_depth > a.sum + b.relative_depth)
  //           ? a.distance_from_min + b.subtree_size
  //           : b.distance_from_min.
  //
  //     Given those prefix sums, we can calculate Z(j) as follows:
  //
  //       1) Construct the prefix_sum for the elements 0 to j.
  //
  //       2a) If the relative_depth equals 0, then we want the
  //           distance_from_min..
  //
  //       2b) Otherwise, there are no 0 (which I don't think can happen).

  // For each string in lru_stack_, what is it's rank in lru_stack_.
  std::unordered_map<std::string, TimeStamp> stack_positions_;
  // All these arrays are indexed from 0, but the paper indexes from 1.
  //
  // This is sorted by the negative of the access time, so that the top of the
  // stack is first.
  OrderStatisticTree<TimeStamp, std::string> lru_stack_;
  OrderStatisticTree<TimeStamp, size_t> gamma_;
  std::unordered_map<size_t, TimeStamp> gamma_inverted_;
  TimeStamp timestep_; // Used to index items in lru_stack_ and gamma_.
  // Critical markers are are the same as what the paper says, but the indexes
  // are less.  For the paper's M(1) we store in critical_markers_[0].
  std::vector<size_t> critical_markers_;
  std::vector<size_t> beta_; // to be deprecated
  PrefixTree<SizetComparesBackward, ptrdiff_t> beta_diff_;
  SizetComparesBackward beta_counter_;
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
    csa_.Validate();
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
