#include <iostream>
template <size_t N>
class StackStorage {
 public:
  StackStorage() : top(memory) {};
  StackStorage(const StackStorage<N>&) = delete;
  StackStorage<N>& operator=(const StackStorage<N>&) = delete;

  template <typename T>
  T* allocate(size_t count) {
    std::align(alignof(T), sizeof(T), reinterpret_cast<void*&>(top), space);
    char* old_top = top;
    top += count * sizeof(T);
    return reinterpret_cast<T*>(old_top);
  }

 private:
  char memory[N];
  char* top;
  size_t space = N;
};
template <typename T, size_t N>
class StackAllocator {
 public:
  using value_type = T;

  StackAllocator()
      : stack_storage(new StackStorage<N>)
      , use_heap_memory(true) {}

  explicit StackAllocator(StackStorage<N>& other) : stack_storage(&other) {
  }
  template <typename U>
  StackAllocator(const StackAllocator<U, N>& other)
      : stack_storage(other.get_storage()) {
  }

  template <typename U>
  StackAllocator<T, N>& operator=(const StackAllocator<U, N>& other) {
    stack_storage = other.stack_storage;
  }

  T* allocate(size_t count) {
    return stack_storage->template allocate<T>(count);
  }

  void deallocate(T*, size_t) {}

  struct propagate_on_container_copy_assignment : std::true_type {};

  template <typename U>
  struct rebind {
    using other = StackAllocator<U, N>;
  };

  StackStorage<N>* get_storage() const { return stack_storage; }
  ~StackAllocator() {
    if (use_heap_memory) {
      delete stack_storage;
    }
  }
 private:
  StackStorage<N>* stack_storage;
  size_t space = N;
  bool use_heap_memory = false;
};

template <typename T, typename U, size_t N>
bool operator==(const StackAllocator<T, N>& first,
                const StackAllocator<U, N>& second) {
  return first.get_storage() == second.get_storage();
}

template <typename T, typename U, size_t N>
bool operator!=(const StackAllocator<T, N>& first,
                const StackAllocator<U, N>& second) {
  return !(first == second);
}

template <typename T, typename Alloc = std::allocator<T>>
class List {
 private:

  struct BaseNode {
    BaseNode* prev;
    BaseNode* next;
  };
  struct Node : BaseNode {
    T value;
  };

  template <bool IsConst>
  class base_iterator;
  BaseNode fake_node;
  base_iterator<false> start;
  base_iterator<false> finish;
  size_t sz;
  using NodeAlloc = std::allocator_traits<Alloc>::template rebind_alloc<Node>;
  using NodeTraits = std::allocator_traits<NodeAlloc>;
  [[no_unique_address]] NodeAlloc alloc;
  void Destroy() {
    for (size_t i = 0; i < sz; ++i) {
      BaseNode* delete_node = fake_node.next;
      fake_node.next = delete_node->next;
      NodeTraits::destroy(alloc, &static_cast<Node*>(delete_node)->value);
      NodeTraits::deallocate(alloc, static_cast<Node*>(delete_node), 1);
    }
    sz = 0;
  }

  void push_back_default() {
    Node* new_node = NodeTraits::allocate(alloc, 1);
    try {
      NodeTraits::construct(alloc, &new_node->value);
      BaseNode* end_node = fake_node.prev;
      end_node->next = static_cast<BaseNode*>(new_node);
      new_node->next = &fake_node;
      new_node->prev = end_node;
      fake_node.prev = static_cast<BaseNode*>(new_node);
      start.set_iter_node(fake_node.next);
      finish.set_iter_node(&fake_node);
      ++sz;
    } catch (...) {
      NodeTraits::deallocate(alloc, new_node, 1);
      throw;
    }
  }

 public:
  using iterator = base_iterator<false>;
  using const_iterator = base_iterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  iterator begin() { return start; };
  iterator end() { return finish; };
  const_iterator begin() const { return start; };
  const_iterator end() const { return finish; };
  const_iterator cbegin() const { return start; };
  const_iterator cend() const { return finish; }
  reverse_iterator rbegin() { return std::reverse_iterator(finish); };
  reverse_iterator rend() { return std::reverse_iterator(start); };
  const_reverse_iterator rbegin() const {
    return std::reverse_iterator(finish);
  };
  const_reverse_iterator rend() const {
    return std::reverse_iterator(start);
  };
  const_reverse_iterator crbegin() const {
    return std::reverse_iterator(finish);
  };
  const_reverse_iterator crend() const { return std::reverse_iterator(start); };
  List()
      : fake_node{&fake_node, &fake_node}
      , start(&fake_node)
      , finish(&fake_node)
      , sz(0) {}

  explicit List(const Alloc& other_alloc)
      : fake_node{&fake_node, &fake_node}
      , start(&fake_node)
      , finish(&fake_node)
      , sz(0)
      , alloc(other_alloc) {}

  List(size_t count, const T& value, const Alloc& other_alloc)
      : fake_node{&fake_node, &fake_node}
      , start(&fake_node)
      , finish(&fake_node)
      , sz(0)
      , alloc(other_alloc) {
    try {
      for (size_t i = 0; i < count; ++i) {
        push_back(value);
      }
    } catch (...) {
      Destroy();
      throw;
    }
  }

  explicit List(size_t count, const Alloc& other_alloc)
      : fake_node{&fake_node, &fake_node}
      , start(&fake_node)
      , finish(&fake_node)
      , sz(0)
      , alloc(other_alloc) {
    try {
      for (size_t i = 0; i < count; ++i) {
        push_back_default();
      }
    } catch (...) {
      Destroy();
      throw;
    }
  }

  explicit List(size_t count)
      : fake_node{&fake_node, &fake_node}
      , start(&fake_node)
      , finish(&fake_node)
      , sz(0) {
    try {
      for (size_t i = 0; i < count; ++i) {
        push_back_default();
      }
    } catch (...) {
      Destroy();
      throw;
    }
  }

  List(const List<T, Alloc>& other)
      : fake_node{&fake_node, &fake_node}
      , start(&fake_node)
      , finish(&fake_node)
      , sz(0)
      , alloc(NodeTraits::select_on_container_copy_construction(other.alloc)) {
    try {
      for (List<T, Alloc>::const_iterator it = other.begin(); it != other.end();
           ++it) {
        push_back(*it);
      }
    } catch (...) {
      Destroy();
      throw;
    }
  }

  List<T, Alloc>& operator=(const List<T, Alloc>& other) {
    try {
      if (this != &other) {
        if (std::is_base_of_v<std::true_type,
                              typename std::allocator_traits<Alloc>::propagate_on_container_copy_assignment>) {
          alloc = other.alloc;
        }
        if (sz <= other.sz) {
          List<T, Alloc>::iterator it = begin();
          List<T, Alloc>::const_iterator other_it = other.begin();
          for (size_t i = 0; i < sz; ++i, ++it, ++other_it) {
            *it = *other_it;
          }
          for (; other_it != other.end(); ++other_it) {
            push_back(*other_it);
          }
        } else {
          size_t delta = sz - other.sz;
          for (size_t i = 0; i < delta; ++i) {
            pop_back();
          }
          List<T, Alloc>::iterator it = begin();
          List<T, Alloc>::const_iterator other_it = other.begin();
          for (size_t i = 0; i < sz; ++i, ++it, ++other_it) {
            *it = *other_it;
          }
        }
      }
    } catch (...) {
      throw;
    }
    return *this;
  }

  size_t size() const { return sz; }

  void push_back(const T& value) {
    Node* new_node = NodeTraits::allocate(alloc, 1);
    try {
      NodeTraits::construct(alloc, &new_node->value, value);
      BaseNode* end_node = fake_node.prev;
      end_node->next = static_cast<BaseNode*>(new_node);
      new_node->next = &fake_node;
      new_node->prev = end_node;
      fake_node.prev = static_cast<BaseNode*>(new_node);
      ++sz;
      start.set_iter_node(fake_node.next);
      finish.set_iter_node(&fake_node);
    } catch (...) {
      NodeTraits::deallocate(alloc, new_node, 1);
      throw;
    }
  }

  void pop_back() {
    BaseNode* end_node = fake_node.prev;
    BaseNode* prev_end_node = end_node->prev;
    prev_end_node->next = &fake_node;
    fake_node.prev = prev_end_node;
    NodeTraits::destroy(alloc, &static_cast<Node*>(end_node)->value);
    NodeTraits::deallocate(alloc, static_cast<Node*>(end_node), 1);
    --sz;
    start.set_iter_node(fake_node.next);
    finish.set_iter_node(&fake_node);
  }

  void push_front(const T& value) {
    Node* new_node = NodeTraits::allocate(alloc, 1);
    try {
      NodeTraits::construct(alloc, &new_node->value, value);
      BaseNode* begin_node = fake_node.next;
      begin_node->prev = static_cast<BaseNode*>(new_node);
      fake_node.next = static_cast<BaseNode*>(new_node);
      new_node->next = begin_node;
      new_node->prev = &fake_node;
      ++sz;
      start.set_iter_node(fake_node.next);
      finish.set_iter_node(&fake_node);
    } catch (...) {
      NodeTraits::deallocate(alloc, new_node, 1);
      throw;
    }
  }

  void pop_front() {
    BaseNode* begin_node = fake_node.next;
    BaseNode* next_begin_node = begin_node->next;
    next_begin_node->prev = &fake_node;
    fake_node.next = next_begin_node;
    NodeTraits::destroy(alloc, &static_cast<Node*>(begin_node)->value);
    NodeTraits::deallocate(alloc, static_cast<Node*>(begin_node), 1);
    --sz;
    start.set_iter_node(fake_node.next);
    finish.set_iter_node(&fake_node);
  }

  iterator insert(const_iterator pos, const T& value) {
    Node* new_node = NodeTraits::allocate(alloc, 1);
    try {
      NodeTraits::construct(alloc, &new_node->value, value);
      BaseNode* prev_node = pos.get_iter_node()->prev;
      prev_node->next = static_cast<BaseNode*>(new_node);
      new_node->prev = prev_node;
      pos.get_iter_node()->prev = static_cast<BaseNode*>(new_node);
      new_node->next = pos.get_iter_node();
      ++sz;
      start.set_iter_node(fake_node.next);
      finish.set_iter_node(&fake_node);
      return iterator(new_node);
    } catch (...) {
      NodeTraits::deallocate(alloc, new_node, 1);
      throw;
    }

  }

  iterator erase(const_iterator pos) {
    BaseNode* pos_node = pos.get_iter_node();
    BaseNode* prev_pos = pos_node->prev;
    BaseNode* next_pos = pos_node->next;
    prev_pos->next = next_pos;
    next_pos->prev = prev_pos;
    NodeTraits::destroy(alloc, &static_cast<Node*>(pos_node)->value);
    NodeTraits::deallocate(alloc, static_cast<Node*>(pos_node), 1);
    --sz;
    start.set_iter_node(fake_node.next);
    finish.set_iter_node(&fake_node);
    return iterator(next_pos);
  }

  Alloc get_allocator() const { return Alloc(alloc); }

  ~List() { Destroy(); }

};

template <typename T, typename Alloc>
template <bool IsConst>
class List<T, Alloc>::base_iterator {
 public:
  using iterator_category = std::bidirectional_iterator_tag;
  using pointer = std::conditional_t<IsConst, const T*, T*>;
  using reference = std::conditional_t<IsConst, const T&, T&>;
  using value_type = T;
  using difference_type = std::ptrdiff_t;

  explicit base_iterator(BaseNode* other_iter_node)
      : iter_node(other_iter_node) {}

  base_iterator<IsConst>& operator++() {
    iter_node = iter_node->next;
    return *this;
  }

  base_iterator<IsConst>& operator--() {
    iter_node = iter_node->prev;
    return *this;
  }

  base_iterator<IsConst> operator++(int) {
    base_iterator<IsConst> it = *this;
    iter_node = iter_node->next;
    return it;
  }

  base_iterator<IsConst> operator--(int) {
    base_iterator<IsConst> it = *this;
    iter_node = iter_node->prev;
    return it;
  }

  bool operator==(const base_iterator<IsConst>& it) const {
    return iter_node == it.iter_node;
  }

  bool operator!=(const base_iterator<IsConst>& it) const {
    return !(*this == it);
  }

  reference operator*() const { return static_cast<Node*>(iter_node)->value; }

  pointer operator->() const { return &(static_cast<Node*>(iter_node)->value); }

  operator base_iterator<true>() const {
    base_iterator<true> New(iter_node);
    return New;
  }

  BaseNode* get_iter_node() const { return iter_node; }

  void set_iter_node(BaseNode* new_iter_node) {
    iter_node = new_iter_node;
  }

 private:
  BaseNode* iter_node;
};