// Critical Marker Algorithm from BilardiEkPa17
//
// This is the algorithm that runs in time O(nm) where n is the number of
// distinct addresses and m is the length of the trace.  This algorithm is
// mostly useful for explaining, since it's not a fast algorithm.

#include <cassert>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "ost.h"

static bool verbose = false;
static const bool expensive = true;

template <class T, class U>
std::ostream& operator<<(std::ostream &os, const std::pair<T, U> &v);

template <class T>
std::ostream& operator<<(std::ostream &os, const std::optional<T> &v) {
  if (v) return os << *v;
  else return os << "nullopt";
}

template <class Container>
std::ostream& PrintContainer(std::ostream &os, const Container &container) {
  os << "{";
  bool first = true;
  for (const auto &t : container) {
    if (!first) os << ", ";
    first = false;
    os << t;
  }
  return os << "}";
}

template <class T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &v) {
  return PrintContainer(os, v);
}

template <class T, class V>
std::ostream& operator<<(std::ostream &os, const std::map<T, V> &v) {
  return PrintContainer(os, v);
}

template <class T, class U>
std::ostream& operator<<(std::ostream &os, const std::pair<T, U> &v) {
  const auto& [a, b] = v;
  return os << "{" << a << ", " << b << "}";
}

template <class T, class V>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<T, V> &map) {
  return PrintContainer(os, map);
}

std::ostream& operator<<(std::ostream& os, const std::optional<size_t> &v) {
  if (v == std::nullopt)
    return os << "nullopt";
  else
    return os << *v;
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
  std::pair</*OPT depth*/std::optional<size_t>,
            /*LRU depth*/std::optional<size_t>> Access(const std::string &t) {
    const std::optional<size_t> depth = FindDepth(t);
    most_recent_z_ = std::nullopt;
    most_recent_dopt_ = std::nullopt;
    if (!depth.has_value()) {
      lru_stack_.insert(lru_stack_.begin(), t);
      IncrementAllCriticalMarkers();
      critical_markers_.push_back(1);
      RecomputeBeta();
      return {std::nullopt, std::nullopt};
    } else if (depth == 0) {
      // state remains the same
      return {0, 0};
    } else {
      if (verbose) std::cout << "D=" << *depth << std::endl;
      const size_t z = FindZ(*depth);
      if (verbose) std::cout << "z=" << z << std::endl;
      const size_t depth_opt = FindDepthOpt(z);
      if (verbose) std::cout << "depth_opt=" << depth_opt << std::endl;
      RotateLruStack(*depth);
      for (size_t i = 1; i < critical_markers_.size(); ++i) {
        if (critical_markers_[i] <= depth) ++critical_markers_[i];
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

struct BetaPrefix {
  static constexpr bool IsPrintable = true;
  BetaPrefix() :BetaPrefix(0) {}
  explicit BetaPrefix(ptrdiff_t increment)
      :sum(increment), min(increment),
       distance_from_min(increment <= 0 ? 0 : 1), tree_size(1) {}
  BetaPrefix operator+(const BetaPrefix &b) const {
    BetaPrefix result;
    result.sum = sum + b.sum;
    result.min = std::min(min, sum + b.min);
    result.distance_from_min =
        (min < sum + b.min)
        ? (distance_from_min + b.tree_size)
        : b.distance_from_min;
    result.tree_size = tree_size + b.tree_size;
    return result;
  }

  friend std::ostream& operator<<(std::ostream& out, const BetaPrefix& n) {
    return out << "{s=" << n.sum << " min=" << n.min << " dmin="
               << n.distance_from_min << " ts=" << n.tree_size << "}";
  }
  bool operator==(const BetaPrefix &other) {
    return sum == other.sum && min == other.min
        && distance_from_min == other.distance_from_min
        && tree_size == other.tree_size;
  }

  ptrdiff_t sum;
  ptrdiff_t min;
  size_t    distance_from_min;
  size_t    tree_size;
};

void TestBetaPrefix() {
  // Beta is {0,1,1,0,0}
  const BetaPrefix a(0);
  const BetaPrefix b(1);
  const BetaPrefix c(0);
  const BetaPrefix d(-1);
  const BetaPrefix e(0);

  const BetaPrefix sum1 = a + b;
  assert(sum1.sum == 1);
  assert(sum1.min == 0);
  assert(sum1.distance_from_min == 1);
  assert(sum1.tree_size == 2);

  const BetaPrefix sum2 = sum1 + c;
  assert(sum2.sum == 1);
  assert(sum2.min == 0);
  assert(sum2.distance_from_min == 2);
  assert(sum2.tree_size == 3);

  const BetaPrefix sum3 = sum2 + d;
  assert(sum3.sum == 0);
  assert(sum3.min == 0);
  assert(sum3.distance_from_min == 0);
  assert(sum3.tree_size == 4);

  const BetaPrefix sum4 = sum3 + e;
  assert(sum4.sum == 0);
  assert(sum4.min == 0);
  assert(sum4.distance_from_min == 0);
  assert(sum4.tree_size == 5);

  // Now check associativity.
  assert((a + b) + c == sum2);
  assert(a + (b + c) == sum2);
  assert(((a + b) + c) + d == sum3);
  assert((a + (b + c)) + d == sum3);
  assert((a + b) + (c + d) == sum3);
  assert(a + (b + (c + d)) == sum3);
  assert(a + ((b + c) + d) == sum3);


}


class Csa {
 public:
  std::pair</*OPT depth*/std::optional<size_t>,
            /*LRU depth*/std::optional<size_t>> Access(const std::string &t) {
    if (expensive) Validate();
    if (verbose) {
      std::cout << std::endl << "CSA access " << t << std::endl;
      std::cout << *this << std::endl;
    }
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
      return {std::nullopt, std::nullopt};
    } else if (depth == 0) {
      // state remains the same
      return {0, 0};
    } else {
      if (verbose) std::cout << "D=" << *depth << std::endl;
      const size_t z = FindZ(*depth);
      const size_t depth_opt = FindDepthOpt(z);
      if (verbose) std::cout << "depth_opt=" << depth_opt << std::endl;
      size_t M_star = critical_markers_[depth_opt];
      if (verbose) std::cout << "M*=" << M_star << std::endl;
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
        if (critical_markers_[i] <= depth) ++critical_markers_[i];
      }
      critical_markers_[depth_opt] = 1;
      if (verbose) std::cout << "Just before computing beta:" << std::endl << *this << std::endl;
      RecomputeBeta();
      if (verbose) {
        std::cout << "Just after  computing beta:" << std::endl << *this << std::endl;
        std::cout << "depth_opt=" << depth_opt << std::endl;
        std::cout << "beta_diff_=" << beta_diff_ << std::endl;
      }
      // BilardiEkPa17 Equation 19:
      //   Decrement beta_[i] for all i > M*
      //   where M* = critical_markers_[depth_opt - 1].
      {
        const auto [k, v] = beta_diff_.Select(M_star);
        if (verbose) std::cout << "select found k=" << k << std::endl;
        beta_diff_.InsertOrAssign(k, v-1);
        if (verbose) std::cout << "-1 beta_diff_=" << beta_diff_ << std::endl;
      }
      // BilardiEkPa17 EQuation 20:
      //   Shift Beta[i] to Beta[i+1] for i < D
      //   Increment Beta[i] for all i > D
      //   where D is the LRU stack depth, represented by *depth here.
      {
        beta_diff_.InsertOrAssign(beta_counter_++, 0);
        if (verbose) std::cout << "Shifted beta_diff_=" << beta_diff_ << std::endl;
        if (verbose) std::cout << "Incrementing beta_diff_[" << *depth << "] Delta=" << *depth + 1 << std::endl;
        const auto [k1, v1] = beta_diff_.Select(*depth + 1);
        beta_diff_.Erase(k1);
        if (verbose) std::cout << "erased slot " << " beta_diff_=" << beta_diff_ << std::endl;
        if (*depth + 1 < beta_diff_.Size()) {
          const auto [k, v] = beta_diff_.Select(*depth + 1);
          beta_diff_.InsertOrAssign(k, v1 + v + 1);
          if (verbose) std::cout << "+1 slot " << *depth << " beta_diff_=" << beta_diff_ << std::endl;
        } else {
          if (verbose) std::cout << "Just let the tail fall off:" << beta_diff_ << std::endl;
        }
        {
          auto final_sum = beta_diff_.SelectPrefix(beta_diff_.Size() - 1);
          assert(final_sum.sum == 0);
        }
        assert(beta_diff_.Size() == lru_stack_.Size());
      }
      most_recent_z_ = z;
      most_recent_dopt_ = depth_opt;
      if (expensive) Validate();
      return {depth_opt, *depth};
    }
  }
  const OrderStatisticTree<TimeStamp, std::string> &GetL() const { return lru_stack_; }
  const std::vector<size_t> &GetM() const { return critical_markers_; }
  const std::vector<size_t> &GetBeta() const { return beta_; }
  std::vector<size_t> GetGamma() const {
    std::vector<size_t> result;
    for (size_t i = 0; i < gamma_.Size(); ++i) {
      auto [timestamp, crit] = gamma_.Select(i);
      result.push_back(crit);
    }
    return result;
  }
  friend std::ostream& operator<<(std::ostream&os, const Csa &csa) {
    os << "L=" << csa.lru_stack_ << std::endl;
    os << "G=" << csa.gamma_ << std::endl;
    os << "Ginv=" << csa.gamma_inverted_ << std::endl;
    os << "Phi should be={";
    for (size_t i = 0; i < csa.gamma_.Size(); ++i) {
      if (i != 0) os << ", ";
      auto [time, capacity] = csa.gamma_.Select(i);
      assert(capacity > 0);
      os << csa.critical_markers_[capacity-1];
    }
    std::cout << "}" << std::endl;
    os << "S=" << csa.stack_positions_ << std::endl;
    os << "M=" << csa.critical_markers_ << std::endl;
    os << "B=" << csa.beta_ << std::endl;
    os << "Bd=" << csa.beta_diff_;
    return os;
  }
  void Validate() const {
    for (size_t i = 0; i < beta_.size(); ++i) {
      if (beta_[i] != static_cast<size_t>(beta_diff_.SelectPrefix(i).sum)) {
        std::cout << "just before assertion failed:" << std::endl << *this << std::endl;
        std::cout << "beta_[" << i << "]=" << beta_[i] << std::endl;
        std::cout << "beta_diff_ prefix (" << i << ")=" << beta_diff_.SelectPrefix(i) << std::endl;
        assert(beta_[i] == static_cast<size_t>(beta_diff_.SelectPrefix(i).sum));
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
    // Note: depth + 1 = the paper's Delta_t.
    if (verbose) {
      std::cout << "FindZ(" << depth << ") (delta=" << depth + 1 << ")" << std::endl;
      std::cout << "beta=" << beta_ << std::endl;
      std::cout << "beta_diff=" << beta_diff_ << std::endl;
      std::cout << "beta_diff internals:" << std::endl;
      beta_diff_.PrintTreeStructure(std::cout) << std::endl;
    }
    size_t result = 0;
    for (size_t i = 0; i < depth; ++i) {
      if (beta_[i] == 0) result = i;
    }
    if (verbose) std::cout << "FindZ(" << depth << ")" << "=" << result << std::endl;
    BetaPrefix prefix = beta_diff_.SelectPrefix(depth - 1);
    if (verbose) std::cout << "beta_prefix=" << prefix << std::endl;
    assert(result == prefix.tree_size - prefix.distance_from_min - 1);
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
  PrefixTree<SizetComparesBackward, ptrdiff_t, BetaPrefix> beta_diff_;
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
  std::pair</*OPT depth*/std::optional<size_t>,
            /*LRU depth*/std::optional<size_t>> Access(const std::string &t) {
    auto cma_result = cma_.Access(t);
    auto csa_result = csa_.Access(t);
    //std::cout << "cma result = " << cma_result << std::endl;
    //std::cout << "csa result = " << csa_result << std::endl;
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
  const std::vector<size_t> &GetBeta() const {
    auto& csa_beta = csa_.GetBeta();
    auto& cma_beta = cma_.GetBeta();
    CHECK_EQ(csa_beta, cma_beta);
    return cma_beta;
  }
  std::vector<size_t> GetGamma() const {
    return csa_.GetGamma();
  }
  std::optional<size_t> GetMostRecentZ() const { return cma_.most_recent_z_; }
  std::optional<size_t> GetMostRecentDopt() const { return cma_.most_recent_dopt_; }
  friend std::ostream& operator<<(std::ostream&os, const CsaAndCma &both) {
    return os << both.cma_;
  }

 private:
  Cma cma_;
  Csa csa_;
};

using TI = std::pair<std::optional<size_t>, std::optional<size_t>>;
using VI = std::vector<size_t>;
using VS = std::vector<std::string>;

static void TestABCABC() {
  CsaAndCma cnc;
  CHECK_EQ(cnc.Access("a"), TI({{}, {}}));
  CHECK_EQ(cnc.GetL(), {"a"});
  CHECK_EQ(cnc.GetM(), {1});
  CHECK_EQ(cnc.GetBeta(), {0});
  CHECK_EQ(cnc.GetGamma(), {1});
  // Phi should be {1}

  CHECK_EQ(cnc.Access("b"), TI({{}, {}}));
  CHECK_EQ(cnc.GetL(), VS({"b", "a"}));
  CHECK_EQ(cnc.GetM(), VI({2, 1}));
  CHECK_EQ(cnc.GetBeta(), VI({0, 0}));
  CHECK_EQ(cnc.GetGamma(), VI({2, 1}));
  // Phi should be {1, 2}

  CHECK_EQ(cnc.Access("c"), TI({{}, {}}));
  CHECK_EQ(cnc.GetL(), VS({"c", "b", "a"}));
  CHECK_EQ(cnc.GetM(), VI({3, 2, 1}));
  CHECK_EQ(cnc.GetBeta(), VI({0, 0, 0}));
  CHECK_EQ(cnc.GetGamma(), VI({3, 2, 1}));
  // Phi should be {1, 2, 3}

  CHECK_EQ(cnc.Access("a"), TI({{1}, {2}}));
  CHECK_EQ(cnc.GetL(), VS({"a", "c", "b"}));
  CHECK_EQ(cnc.GetM(), VI({3, 1, 2}));
  CHECK_EQ(cnc.GetBeta(), VI({0, 0, 0}));
  CHECK_EQ(cnc.GetGamma(), VI({1, 3, 2}));
  // Phi should be {3, 2, 1}

  CHECK_EQ(cnc.Access("b"), TI({{2}, {2}}));
  CHECK_EQ(cnc.GetL(), VS({"b", "a", "c"}));
  CHECK_EQ(cnc.GetM(), VI({3, 2, 1}));
  CHECK_EQ(cnc.GetBeta(), VI({0, 0, 0}));
  CHECK_EQ(cnc.GetGamma(), VI({2, 1, 3}));
  // Phi should be {2, 3, 1}

  CHECK_EQ(cnc.Access("c"), TI({{1}, {2}}));
  CHECK_EQ(cnc.GetL(), VS({"c", "b", "a"}));
  CHECK_EQ(cnc.GetM(), VI({3, 1, 2}));
  CHECK_EQ(cnc.GetBeta(), VI({0, 0, 0}));
  CHECK_EQ(cnc.GetGamma(), VI({1, 2, 3}));
  // Phi should be {3, 1, 2}

  CHECK_EQ(cnc.Access("a"), TI({{2}, {2}}));
  CHECK_EQ(cnc.GetL(), VS({"a", "c", "b"}));
  CHECK_EQ(cnc.GetM(), VI({3, 2, 1}));
  CHECK_EQ(cnc.GetBeta(), VI({0, 0, 0}));
  CHECK_EQ(cnc.GetGamma(), VI({2, 1, 3}));
  // Phi should be {2, 3, 1}

}


static void TestABCDEB() {
  CsaAndCma cma;
  CHECK_EQ(cma.Access("a"), TI({{}, {}}));
  CHECK_EQ(cma.GetL(), {"a"});
  CHECK_EQ(cma.GetM(), {1});
  CHECK_EQ(cma.GetBeta(), {0});

  CHECK_EQ(cma.Access("b"), TI({{}, {}}));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"b", "a"}));
  CHECK_EQ(cma.GetM(), VI({2, 1}));
  CHECK_EQ(cma.GetBeta(), VI({0, 0}));

  CHECK_EQ(cma.Access("c"), TI({{}, {}}));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"c", "b", "a"}));
  CHECK_EQ(cma.GetM(), VI({3, 2, 1}));
  CHECK_EQ(cma.GetBeta(), VI({0, 0, 0}));

  CHECK_EQ(cma.Access("d"), TI({{}, {}}));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"d", "c", "b", "a"}));
  CHECK_EQ(cma.GetM(), VI({4, 3, 2, 1}));
  CHECK_EQ(cma.GetBeta(), VI({0, 0, 0, 0}));

  CHECK_EQ(cma.Access("e"), TI({{}, {}}));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"e", "d", "c", "b", "a"}));
  CHECK_EQ(cma.GetM(), VI({5, 4, 3, 2, 1}));
  CHECK_EQ(cma.GetBeta(), VI({0, 0, 0, 0, 0}));

  //std::cout << cma << std::endl;
  CHECK_EQ(cma.Access("b"), TI({1ul, 3ul}));
  // In the original paper, this would be 3, but since we are using 1-based
  // indexing it's two.
  CHECK_EQ(cma.GetMostRecentZ(), std::optional<size_t>(2ul));
  CHECK_EQ(cma.GetMostRecentDopt(), std::optional<size_t>(1ul));
  CHECK_EQ(cma.GetL(), std::vector<std::string>({"b", "e", "d", "c", "a"}));
  CHECK_EQ(cma.GetM(), VI({5, 1, 4, 3, 2}));
}

class Opt {
 public:
  void Access(std::string a) {
    auto prev = prev_access_.find(a);
    if (prev != prev_access_.end()) {
      trace_[prev->second].second = timestep_;
    }
    trace_.push_back({a, std::nullopt});
    prev_access_[a] = timestep_;
    ++timestep_;

  }
  std::vector<std::optional<size_t>> GetStackDepths() const {
    std::vector<std::optional<size_t>> result;
    std::vector<std::pair<std::string, std::optional<size_t>>> opt_stack;
    for (size_t time = 0; time < trace_.size(); ++time) {
      auto [address, next] = trace_[time];
      auto it = opt_stack.begin();
      size_t count = 0;
      for (; it != opt_stack.end(); ++it, ++count) {
        if (it->first == address) break;
      }
      if (it == opt_stack.end()) {
        ++count;
        result.push_back(std::nullopt);
      } else {
        result.push_back(count);
        opt_stack.erase(it);
      }
      opt_stack.insert(opt_stack.begin(), {address, next});
      // One pass of insertion sort.  But the item in location 0 must stay
      if (0) std::cout << "size=" << opt_stack.size() << " count=" << count << std::endl;
      for (size_t i = 1; i + 1 < opt_stack.size() && i < count; ++i) {
        if (0) std::cout << " consider swap " << i << " with " << i + 1 << std::endl;
        if ((opt_stack[i].second
             && opt_stack[i+1].second
             && *opt_stack[i].second > *opt_stack[i+1].second)
            || (!opt_stack[i].second
                && opt_stack[i+1].second)) {
          if (0) std::cout << " swap  " << i << " with " << i + 1 << std::endl;
          std::swap(opt_stack[i], opt_stack[i+1]);
        }
      }
      if (verbose) std::cout << "OPT Access " << address << " " << result.back()
                             << " stack=" << opt_stack << std::endl;
    }
    return result;
  }
 private:
  // The second element of the pair is the next access time.
  std::vector<std::pair<std::string, std::optional<size_t>>> trace_;
  size_t timestep_ = 0;
  std::unordered_map<std::string, size_t> prev_access_;
};


void randomtest() {
  for (size_t trial = 0; trial < 10; ++trial) {
    std::vector<std::string> addresses = {"a", "b", "c"};
    std::vector<std::string> trace;
    Opt opt;
    Csa csa;
    Cma cma;
    std::vector<std::optional<size_t>> csa_results, cma_results;
    for (size_t i = 0; i < (verbose ? 20 : 80); ++i) {
      if (i % 20 == 0) {
        addresses.push_back(std::to_string(addresses.size()));
      }
      std::string a = addresses[static_cast<size_t>(random()) % addresses.size()];;
      if (verbose) std::cout << "a(" << a << ")" << std::endl;
      opt.Access(a);
      {
        auto [opt_depth, lru_depth] = csa.Access(a);
        csa_results.push_back(opt_depth);
      }
      {
        auto [opt_depth, lru_depth] = cma.Access(a);
        cma_results.push_back(opt_depth);
        //std::cout << cma << std::endl;
      }
    }
    std::vector<std::optional<size_t>> opt_results = opt.GetStackDepths();
    if (verbose) {
      std::cout << "cma results = " << cma_results << std::endl;
      std::cout << "csa results = " << csa_results << std::endl;
      std::cout << "opt results = " << opt_results << std::endl;
    }
    assert(cma_results == opt_results);
    assert(csa_results == opt_results);
  }
}

int main() {
  TestBeta();
  TestBetaPrefix();
  TestABCABC();
  return 0;
  TestABCDEB();
  const std::vector<std::pair<std::string, std::optional<size_t>>> trace =
      {{"a", {}}, // OPT=a
       {"b", {}}, // OPT=b a
       {"c", {}}, // OPT=c b a  (c must be in cache 0, but is still before a)
       // OPT=b c a   (c is next accessed after b)
       {"d", {}}, // OPT=d b c a   (d is in cache0, but b beats c, and then c beats a)
       {"e", {}}, // OPT=e b c a d
       {"b", 1}, // OPT=b e c a d  (b moves up, and then e doesn't have to fight c)
       // Something goes wrong in my CMA implementation here.
       {"c", 2}, // OPT=c b e a d  (c moves up, b beats e, e doesn't have to fight a)
       {"f", {}}, // OPT=f b c a d e (f moves up, e moves out (but could have had c move out)
       {"a", 3}, // OPT=a b f c d e
       {"b", 1}, // OPT=b a f c d e
       {"g", {}},// OPT=g b a f d c e
       {"d", 4}  // OPT=d g b a f c e
      };
  Opt opt;
  for (auto &[address, depth] : trace) {
    opt.Access(address);
  }
  const std::vector<std::optional<size_t>> opt_result = opt.GetStackDepths();
  assert(opt_result.size() == trace.size());
  verbose=true;
  {
    Csa csa;
    size_t time = 0;
    for (auto &[address, depth] : trace) {
      //std::cout << "Accessing " << address << std::endl;
      auto [c, lru_depth] = csa.Access(address);
      assert(c == opt_result[time]);
      //std::cout << "a=" << address << " c=" << c << std::endl << csa << std::endl;
      // std::cout << std::endl;
      assert(depth == SIZE_MAX || depth == c);
      csa.Validate();
      ++time;
    }
  }
  verbose=false;
  randomtest();
  {
    Cma cma;
    const std::vector<std::string> trace =
        {"b", "a", "d", "a", "c", "b", "e", "f", "e", "c"};
    for (size_t i = 0; i < trace.size(); ++i) {
      cma.Access(trace[i]);
      std::cout << trace[i] << i << std::endl;
      std::cout << cma << std::endl << std::endl;
    }
  }
}
