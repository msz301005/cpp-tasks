#include <vector>
#include <algorithm>
#include <iostream>
template <typename Type>
class Deque {
 private:
  static const std::ptrdiff_t kBlocks = 16;
  std::vector<Type*> blocks_;
  template <bool IsConst>
  class base_iterator;
  base_iterator<false> start;
  base_iterator<false> finish;
  void IfEmptyBlocks(Type value) {
    Type* block = reinterpret_cast<Type*>(new char[kBlocks * sizeof(Type)]);
    new(block) Type(value);
    blocks_.push_back(block);
    start = base_iterator<false>(this, 0, 0);
    finish = base_iterator<false>(this, 0, 1);
  }

  void allocate_start_block(Type value) {
    Type* block = reinterpret_cast<Type*>(new char[kBlocks * sizeof(Type)]);
    new(block + kBlocks - 1) Type(value);
    blocks_[start.get_block() - 1] = block;
    --start;
  }

  void allocate_blocks() {
    std::vector<Type*> copy_blocks(2 * blocks_.size(), nullptr);
    std::copy(blocks_.begin(), blocks_.end(),
              copy_blocks.begin() + blocks_.size());
    blocks_ = copy_blocks;
    /*for(size_t i = 0, j = blocks_.size() / 2; i < copy_blocks.size(); ++i, ++j) {
      //new(blocks_[j]) Type*(copy_blocks[i]);
      //blocks_[j] = copy_blocks[i];
    }*/
  }

  std::vector<Type*> allocate(size_t cap_new) {
    std::vector<Type*> blocks_new;
    try {
      for (size_t i = 0; i < (cap_new + kBlocks - 1) / kBlocks; ++i) {
        blocks_new.push_back(
            reinterpret_cast<Type*>(new char[kBlocks * sizeof(Type)]));
      }
    }
    catch (...) {
      for (auto& block : blocks_new) {
        delete[] reinterpret_cast<char*>(block);
      }
      throw;
    }
    return blocks_new;
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
  const_iterator cbegin() { return start; };
  const_iterator cend() { return finish; }
  reverse_iterator rbegin() { return std::reverse_iterator(finish); };
  reverse_iterator rend() { return std::reverse_iterator(start); };
  const_reverse_iterator rbegin() const {
    return std::reverse_iterator(start);
  };
  const_reverse_iterator rend() const { return std::reverse_iterator(finish); };
  const_reverse_iterator crbegin() { return std::reverse_iterator(start); };
  const_reverse_iterator crend() { return std::reverse_iterator(finish); };

  Deque() : start(this, 0, 0), finish(this, 0, 0) {}

  explicit Deque(int number) : Deque() {
    try {
      for (int number_copy = number; number_copy > 0; number_copy -= kBlocks) {
        blocks_.push_back(
            reinterpret_cast<Type*>(new char[kBlocks * sizeof(Type)]));
        if (std::is_default_constructible_v<Type>) {
          for (size_t i = 0; i < kBlocks; ++i) {
            new(blocks_.back() + i) Type();
          }
        }
        start = base_iterator<false>(this, 0, 0);
        finish = start + number;
      }
    } catch (...) {
      /*for (auto& block : blocks_) {
        for (size_t i = 0; i < kBlocks; ++i) {
          (block + i)->~Type();
        }
        delete[] reinterpret_cast<char*>(block);
      }
      blocks_.clear();*/
      blocks_.pop_back();
      throw;
    }
  }

  Deque(int number, const Type& value) : Deque() {
    for (int i = 0; i < number; ++i) {
      push_back(value);
    }
  }

  Deque(const Deque<Type>& other) : Deque() {
    for (size_t i = 0; i < other.blocks_.size(); ++i) {
      Type* block = reinterpret_cast<Type*>(new char[kBlocks * sizeof(Type)]);
      for (size_t j = 0; j < kBlocks; ++j) {
        new(block + j) Type(other.blocks_[i][j]);
      }
      blocks_.push_back(block);
      if (i == static_cast<size_t>(other.start.get_block())) {
        start = base_iterator<false>(this, i, other.start.get_index());
      }
      if (i == static_cast<size_t>(other.finish.get_block())) {
        finish = base_iterator<false>(this, i, other.finish.get_index());
      }
    }
  }

  Deque<Type>& operator=(const Deque<Type>& other) {
    blocks_.clear();
    for (size_t i = 0; i < other.blocks_.size(); ++i) {
      Type* block = reinterpret_cast<Type*>(new char[kBlocks * sizeof(Type)]);
      for (size_t j = 0; j < kBlocks; ++j) {
        new(block + j) Type(other.blocks_[i][j]);
      }
      blocks_.push_back(block);
      if (i == static_cast<size_t>(other.start.get_block())) {
        start = base_iterator<false>(this, i, other.start.get_index());
      }
      if (i == static_cast<size_t>(other.finish.get_block())) {
        finish = base_iterator<false>(this, i, other.finish.get_index());
      }
    }
    return *this;
  }

  void push_back(Type value) {
    if (blocks_.empty()) {
      IfEmptyBlocks(value);
    } else if (finish.get_index() != kBlocks - 1) {
      new(blocks_[finish.get_block()] + finish.get_index()) Type(value);
      ++finish;
    } else {
      new(blocks_[finish.get_block()] + finish.get_index()) Type(value);
      blocks_.push_back(
          reinterpret_cast<Type*>(new char[kBlocks * sizeof(Type)]));
      ++finish;
    }
  }

  void push_front(Type value) {
    if (blocks_.empty()) {
      IfEmptyBlocks(value);
    } else if (start.get_index() != 0) {
      --start;
      new(blocks_[start.get_block()] + start.get_index()) Type(value);
    } else if (start.get_block() != 0) {
      allocate_start_block(value);
    } else {
      allocate_blocks();
      start = base_iterator<false>(this, blocks_.size() / 2, start.get_index());
      finish =
          base_iterator<false>(this, finish.get_block() + blocks_.size() / 2,
                               finish.get_index());
      allocate_start_block(value);
    }
  }

  void pop_back() {
    --finish;
    finish->~Type();
  }

  void pop_front() {
    start->~Type();
    ++start;
  }

  void insert(iterator it, const Type& value) {

    if (blocks_.empty()) {
      IfEmptyBlocks(value);
      return;
    }
    Type first_elem = value;
    std::vector<Type> copy;
    try {
      while (it <= finish) {
        if (static_cast<size_t>(it.get_block()) >= blocks_.size()) {
          Type* block =
              reinterpret_cast<Type*>(new char[kBlocks * sizeof(Type)]);
          new(block) Type(first_elem);
          blocks_.push_back(block);
          continue;
        }
        copy.assign(it, it + (kBlocks - it.get_index() - 1));
        *it = first_elem;
        if (finish > it + (kBlocks - it.get_index() - 1)) {
          first_elem = blocks_[it.get_block()][kBlocks - 1];
        }
        std::copy(copy.begin(), copy.end(), it + 1);
        it += kBlocks - it.get_index();
      }
      ++finish;
    } catch (...) {
      copy.clear();
      return;
      // throw;
    }

  }

  void erase(iterator it) {
    while (it < finish) {
      std::vector<Type> copy(it + 1, it + (kBlocks - it.get_index()));
      std::copy(copy.begin(), copy.end(), it);
      if (it.get_block() == finish.get_block()) {
        it = finish;
      } else {
        it += kBlocks - it.get_index();
      }
      if (it < finish) {
        Type elem = blocks_[it.get_block()][0];
        --it;
        *it = elem;
        ++it;
      }
    }
    --finish;
  }

  size_t size() const { return finish - start; }

  Type& operator[](size_t index) { return *(start + index); }

  const Type& operator[](size_t index) const { return *(start + index); }

  Type& at(size_t index) {
    if (index >= size()) {
      throw std::out_of_range("");
    }
    return *(start + index);
  }

  const Type& at(size_t index) const {
    if (index >= size()) {
      throw std::out_of_range("");
    }
    return *(start + index);
  }

  static std::ptrdiff_t get_kBlock() { return kBlocks; }

  ~Deque() {
    for (auto& block : blocks_) {
      for (size_t i = 0; i < kBlocks; ++i) {
        if (std::is_destructible_v<Type>) {
          (block + i)->~Type();
        }
      }
      delete[] reinterpret_cast<char*>(block);
    }
    blocks_.resize(0);
  }
};

template <typename Type>
template <bool IsConst>
class Deque<Type>::base_iterator {
 public:
  using iterator_category = std::random_access_iterator_tag;
  using pointer = std::conditional_t<IsConst, const Type*, Type*>;
  using reference = std::conditional_t<IsConst, const Type&, Type&>;
  using value_type = Type;
  using difference_type = std::ptrdiff_t;
  base_iterator(Deque<Type>* pointer, difference_type block, difference_type index)
      : pointer_deque(pointer)
      , block_(block)
      , index_(index) {};

  // base_iterator<IsConst>& operator=(const base_iterator<IsConst>&) = default;

  base_iterator<IsConst>& operator+=(difference_type delta) {
    while (index_ + delta >= kBlock) {
      ++block_;
      delta -= kBlock - index_;
      index_ = 0;
    }
    while (index_ + delta < 0) {
      --block_;
      delta += index_ + 1;
      index_ = kBlock - 1;
    }
    index_ += delta;
    return *this;
  }

  base_iterator<IsConst>& operator-=(difference_type delta) {
    *this += -delta;
    return *this;
  }

  base_iterator<IsConst> operator+(difference_type delta) const {
    base_iterator<IsConst> res = *this;
    res += delta;
    return res;
  }

  base_iterator<IsConst> operator-(difference_type delta) const {
    base_iterator<IsConst> res = *this;
    res -= delta;
    return res;
  }

  difference_type operator-(const base_iterator<IsConst>& it) const {
    if (block_ == it.block_) {
      return index_ - it.index_;
    }
    difference_type res =
        (std::abs(block_ - it.block_) - 1) * kBlock + kBlock - it.index_ +
        index_;
    return (block_ < it.block_ ? -res : res);
  }

  base_iterator<IsConst>& operator++() { return *this += 1; }

  base_iterator<IsConst>& operator--() { return *this -= 1; }

  base_iterator<IsConst> operator++(int) {
    base_iterator<IsConst> it = *this;
    *this += 1;
    return it;
  }

  base_iterator<IsConst> operator--(int) {
    base_iterator<IsConst> it = *this;
    *this -= 1;
    return it;
  }

  bool operator==(const base_iterator<IsConst>& it) const {
    return *this - it == 0;
  }

  bool operator!=(const base_iterator<IsConst>& it) const {
    return !(*this == it);
  }

  bool operator<(const base_iterator<IsConst>& it) const {
    return it - *this > 0;
  }

  bool operator<=(const base_iterator<IsConst>& it) const {
    return !(it < *this);
  }

  bool operator>(const base_iterator<IsConst>& it) const { return it < *this; }

  bool operator>=(const base_iterator<IsConst>& it) const {
    return !(*this < it);
  }

  reference operator*() const { return pointer_deque->blocks_[block_][index_]; }

  pointer operator->() const { return &pointer_deque->blocks_[block_][index_]; }

  operator base_iterator<true>() const {
    base_iterator<true> New(pointer_deque, block_, index_);
    return New;
  }

  difference_type get_block() const { return block_; }

  difference_type get_index() const { return index_; }

 private:
  Deque<Type>* pointer_deque;
  difference_type block_;
  difference_type index_;
  static const std::ptrdiff_t kBlock = Deque<Type>::kBlocks;
};