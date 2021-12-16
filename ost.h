#ifndef OST_H_
#define OST_H_

class NullStruct {
};

// A combiner allows us to add together a bunch of V's to get an Value.  We use
// this to get prefix sums.

// The nullcombiner doesn't add anything.
template<class V>
struct NullCombiner {
 public:
  static constexpr bool IsPrintable = false;
  NullCombiner() {}
  explicit NullCombiner([[maybe_unused]] const V& v) {}
  using Value = void;
  void GetValue() const {}
  NullCombiner operator+([[maybe_unused]] const NullCombiner b) const {
    return *this;
  }
  void Check([[maybe_unused]] const V& v,
             [[maybe_unused]] const NullCombiner *a,
             [[maybe_unused]] const NullCombiner *b) const {
  }
  friend std::ostream& operator<<(std::ostream& out,
                                  [[maybe_unused]] const NullCombiner& n) {
    return out;
  }
};

// An order statistic tree with no prefix summing.
template<class K, class V, class Combiner = NullCombiner<V>>
    class OrderStatisticTree {
  struct Node;
 public:
  size_t Size() const {
    return SubtreeSize(root_);
  }
  using ValueType = std::pair<const K, V>;
  using RankResult = std::pair<size_t, const ValueType*>;
  // Returns the number of key-value pairs that are less than k.  If there is a
  // key-value pair matching k, the second element is a pointer to a const pair.
  // The pair is const since, to modify it, you must call InsertOrAssign to get
  // the prefix sums updated properly.
  const RankResult Rank(const K& k) const {
    Check();
    return Rank(root_, k);
  }
  void InsertOrAssign(const K& k, const V& v) {
    Check();
    root_ = Insert(root_, k, v);
    Check();
  }
  void Erase(const K& k) {
    Check();
    root_ = Erase(root_, k);
    Check();
  }
  // Requires idx < Size().  Returns a reference to the const key-value pair
  // with rank idx.
  const std::pair<const K, V>& Select(size_t idx) const {
    Check();
    assert(idx < Size());
    const Node* n = Select(root_, idx);
    return n->pair;
  }
  // Returns the combined value the first idx values.  If idx is too big,
  // returns the sum of all the values.
  const Combiner SelectPrefix(size_t idx) const {
    Check();
    if (idx >= Size()) return root_ ? root_->subtree_value : Combiner();
    return SelectPrefix(root_, idx);
  }
  void Check() const {
    Check(root_);
  }
  std::ostream& PrintTreeStructure(std::ostream& out) const {
    out << std::endl << "{";
    return PrintNodeStructure(out, root_, true, 0) << "}" << std::endl;
  }
 private:
  static std::ostream& PrintNodeStructure(std::ostream& out,
                                          const Node *node,
                                          bool need_newline,
                                          size_t depth) {
    if (node) {
      if (node->left) {
        PrintNodeStructure(out, node->left, need_newline, depth+1);
        need_newline = true;
      }
      if (need_newline) out << std::endl << std::string(depth, ' ');
      out << node->k << ":" << node->v << ":" << node->subtree_value;
      if (node->right) {
        PrintNodeStructure(out, node->right, true, depth+1);
      }
    }
    return out;
  }
  friend std::ostream& operator<<(std::ostream& out, const OrderStatisticTree& tree) {
    out << "{";
    for (size_t i = 0; i < tree.Size(); ++i) {
      if (i > 0) out << ", ";
      auto [k, v] = tree.Select(i);
      out << "{" << k << ", " << v;
      if constexpr (Combiner::IsPrintable) {
        out << ", " << tree.SelectPrefix(i);
      }
      out << "}";
      //if (i > 0) out << ", ";
      //out << tree.SelectPrefix(i);
    }
    return out << "}";
  }

#if 0
  friend std::ostream& operator<<(std::ostream& out,
                                  const OrderStatisticTree& tree) {
    return out << "{TREE " << tree.root_ << "}";
  }
  friend std::ostream& operator<<(std::ostream& out, const Node* node) {
    if (node == nullptr) return out;
    return out << node->left << (node->left ? " " : "")
               << node->k << ":" << node->v
               << (node->right ? " " : "") << node->right;
    if (0) {
      return out << "[size=" << node->subtree_size << " L=" << node->left
                 << " K=" << node->k << " V=" << node->v << " R=" << node->right
                 << "]";
    }
  }
#endif
  Node* MaybeRebalance(Node *n) {
    assert(n);
    size_t left_size = 1 + SubtreeSize(n->left);
    size_t right_size = 1 + SubtreeSize(n->right);
    if (left_size * 2 <= right_size && right_size * 2 <= left_size) {
      return n;
    }
    tmp_.clear();
    tmp_.reserve(SubtreeSize(n));
    Tree2Vector(n);
    return Vector2Tree(tmp_.begin(), tmp_.end());
  }
  void Tree2Vector(Node *n) {
    if (n) {
      Tree2Vector(n->left);
      tmp_.push_back(n);
      Tree2Vector(n->right);
    }
  }
  using Viterator = typename std::vector<Node*>::iterator;
  static Node* Vector2Tree(Viterator begin, Viterator end) {
    if (begin >= end) return nullptr;
    Viterator mid = begin + (end - begin) / 2;
    Node *n = *mid;
    n->left = Vector2Tree(begin, mid);
    n->right = Vector2Tree(mid + 1, end);
    RecomputeSubtreeSummary(n);
    Check(n);
    return n;
  }
  static RankResult Rank(const Node *n, const K& k) {
    assert(n != nullptr);
    if (k < GetKey(n)) return Rank(n->left, k);
    else if (GetKey(n) < k) {
      auto [size, pair] = Rank(n->right, k);
      return {1 + SubtreeSize(n->left) + size, pair};
    } else {
      return {SubtreeSize(n->left), &n->pair};
    }
  }
  Node* Insert(Node* n, const K& k, const V& v) {
    if (n == nullptr) {
      n = new Node(k, v);
    } else {
      n = MaybeRebalance(n);
      if (GetKey(n) < k) {
        UpdateRight(n, Insert(n->right, k, v));
        Check(n);
      } else  if (k < GetKey(n)) {
        Check(n);
        Node *newleft = Insert(n->left, k, v);
        UpdateLeft(n, newleft);
        Check(n);
      } else {
        GetVal(n) = v;
        RecomputeCombination(n);
        Check(n);
      }
    }
    Check(n);
    return n;
  }
  static Node* UnlinkRightMost(Node *n,
                               Node **removed_node) {
    assert(n);
    if (n->right) {
      UpdateRight(n, UnlinkRightMost(n->right, removed_node));
      return n;
    } else {
      *removed_node = n;
      return n->left;
    }
  }
  static Node* Erase(Node* n, const K& k) {
    if (n == nullptr) return nullptr;
    else if (k < GetKey(n)) {
      UpdateLeft(n, Erase(n->left, k));
      return n;
    } else if (GetKey(n) < k) {
      UpdateRight(n, Erase(n->right, k));
      return n;
    } else {
      return DeleteNode(n);
    }
  }
  static Node* DeleteNode(Node *n) {
    assert(n);
    if (n->left) {
      Node *new_n;
      Node *new_left = UnlinkRightMost(n->left, &new_n);
      UpdateLeftAndRight(new_n, new_left, n->right);
      delete n;
      Check(new_n);
      return new_n;
    } else {
      Node* result = n->right;
      delete n;
      return result;
    }
  }
  static const Node * Select(const Node *n, size_t idx)
  {
    assert(n);
    assert(idx < SubtreeSize(n));
    if (SubtreeSize(n->left) == idx) {
      return n;
    } else if (SubtreeSize(n->left) > idx) {
      return Select(n->left, idx);
    } else {
      return Select(n->right, idx - SubtreeSize(n->left) - 1ul);
    }
  }
  static Combiner SelectPrefix(const Node *n, size_t idx) {
    if (n == nullptr) return Combiner();
    if (SubtreeSize(n->left) == idx) {
      if (n->left) {
        return n->left->subtree_value + Combiner(GetVal(n));
      } else {
        return Combiner(GetVal(n));
      }
    } else if (SubtreeSize(n->left) > idx) {
      return SelectPrefix(n->left, idx);
    } else {
      return (n->left ? n->left->subtree_value + Combiner(GetVal(n)) : Combiner(GetVal(n)))
          + SelectPrefix(n->right, idx - SubtreeSize(n->left) - 1ul);
    }
  }
  static void UpdateLeftAndRight(Node *n, Node *left, Node *right) {
    assert(n);
    n->left = left;
    n->right = right;
    RecomputeSubtreeSummary(n);
  }
  static void UpdateLeft(Node *n, Node *left) {
    assert(n);
    n->left = left;
    RecomputeSubtreeSummary(n);
    Check(n);
  }
  static void UpdateRight(Node *n, Node *right) {
    assert(n);
    n->right = right;
    RecomputeSubtreeSummary(n);
  }
  static void RecomputeCombination(Node *n) {
    if (n->left) {
      n->subtree_value = n->left->subtree_value + Combiner(GetVal(n));
    } else {
      n->subtree_value = Combiner(GetVal(n));
    }
    if (n->right) {
      n->subtree_value = n->subtree_value + n->right->subtree_value;
    }
  }
  static void RecomputeSubtreeSummary(Node *n) {
    n->subtree_size = 1ul + SubtreeSize(n->left) + SubtreeSize(n->right);
    RecomputeCombination(n);
  }
  static size_t SubtreeSize(const Node *n) {
    return n ? n->subtree_size : 0;
  }
  static const bool check_expensive = false;
  static void Check(const Node *n) {
    if (!n) return;
    assert(n->subtree_size == 1 + SubtreeSize(n->left) + SubtreeSize(n->right));
#if 0
    n->subtree_value.Check(GetVal(n),
                           n->left ? &n->left->subtree_value : nullptr,
                           n->right ? &n->right->subtree_value : nullptr);
#endif
    if (n->left) {
      assert(GetKey(n->left) < GetKey(n));
      if (check_expensive) Check(n->left);
    }
    if (n->right) {
      assert(GetKey(n) < GetKey(n->right));
      if (check_expensive) Check(n->right);
    }
  }
  static const K& GetKey(const Node *n) {
    return n->pair.first;
  }
  static const V& GetVal(const Node *n) {
    return n->pair.second;
  }
  static V& GetVal(Node *n) {
    return n->pair.second;
  }
  struct Node {
    size_t subtree_size = 1;
    Node* left = nullptr;
    Node* right = nullptr;
    std::pair<const K, V> pair;
    Combiner subtree_value;
    Node(const K& k, const V&v) :pair({k, v}), subtree_value(v) {}
  };
  Node *root_ = nullptr;
  // Reuse tmp_ to avoid memory allocations during rebalance.
  std::vector<Node*> tmp_;
};

template<class V>
struct AddCombiner {
  static constexpr bool IsPrintable = true;
  AddCombiner() {}
  explicit AddCombiner(const V& v) :sum(v) {}
  AddCombiner operator+(const AddCombiner b) const {
    return AddCombiner(sum + b.sum);
  }
  friend std::ostream& operator<<(std::ostream& out, const AddCombiner& n) {
    return out << n.sum;
  }

  V sum = 0;
};

template<class K, class V, class Combiner=AddCombiner<V>>
using PrefixTree = OrderStatisticTree<K, V, Combiner>;
#endif  // OST_H_
