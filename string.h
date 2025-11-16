#include <iostream>
#include <cstring>
#include <algorithm>
class String {
 private:
  char* arr = nullptr;
  size_t sz = 0;
  size_t cap = 0;
  explicit String(size_t count) : arr(new char[count + 1]),
                                  sz(count),
                                  cap(count) {
    arr[count] = '\0';
  }

  void swap(String& other) {
    std::swap(arr, other.arr);
    std::swap(sz, other.sz);
    std::swap(cap, other.cap);
  }

  void relocate(size_t delta) {
    String tmp;
    tmp.arr = new char[sz + delta + 1];
    tmp.arr[sz + delta] = '\0';
    tmp.sz = sz;
    tmp.cap = cap + delta;
    swap(tmp);
    std::copy(tmp.arr, tmp.arr + tmp.sz, arr);
  }
 public:

  String() = default;

  String(size_t count, char ch) : String(count) {
    std::fill(arr, arr + count, ch);
  }

  String(const char* str) : String(strlen(str)) {
    std::copy(str, str + sz, arr);
  }

  String(const String& other)
      : arr(new char[other.cap + 1]), sz(other.sz), cap(other.cap) {
    arr[other.cap] = '\0';
    std::copy(other.arr, other.arr + other.sz, arr);
  }

  String& operator=(String other) {
    swap(other);
    return *this;
  }

  char& operator[](size_t ind) { return arr[ind]; }

  const char& operator[](size_t ind) const { return arr[ind]; }

  size_t size() const {
    return sz;
  }

  size_t capacity() const {
    return cap;
  }

  size_t length() const {
    return sz;
  }

  char& front() {
    return arr[0];
  }

  const char& front() const {
    return arr[0];
  }

  char& back() {
    return arr[sz - 1];
  }

  const char& back() const {
    return arr[sz - 1];
  }

  char* data() {
    return arr;
  }

  const char* data() const {
    return arr;
  }

  void pop_back() {
    --sz;
  }

  void push_back(char ch) {
    if (sz == cap) {
      relocate(sz);
    }
    arr[sz++] = ch;
  }

  bool empty() const {
    return !sz;
  }

  void clear() {
    sz = 0;
  }

  size_t find(const String& sub) const {
    for (size_t i = 0; i + sub.sz - 1 < sz; i++) {
      bool flag = true;
      for (size_t j = 0; j < sub.sz; j++) {
        if (arr[i + j] != sub.arr[j]) {
          flag = false;
          break;
        }
      }
      if (flag) {
        return i;
      }
    }
    return sz;
  }

  size_t rfind(const String& sub) const {
    size_t pos = sz;
    for (size_t i = 0; i + sub.sz - 1 < sz; i++) {
      bool flag = true;
      for (size_t j = 0; j < sub.sz; j++) {
        if (arr[i + j] != sub.arr[j]) {
          flag = false;
          break;
        }
      }
      if (flag) {
        pos = i;
      }
    }
    return pos;
  }

  String substr(size_t start, size_t count) const {
    String result;
    result.relocate(count);
    for (auto i = start; i < std::min(start + count, sz); ++i) {
      result += arr[i];
    }
    return result;
  }

  void shrink_to_fit() {
    String tmp;
    tmp.arr = new char[sz + 1];
    tmp.arr[sz] = '\0';
    tmp.sz = tmp.cap = sz;
    swap(tmp);
    std::copy(tmp.arr, tmp.arr + tmp.sz, arr);
  }

  String& operator+=(const String& other) {
    if (cap < sz + other.sz) {
      relocate(other.sz - cap + sz);
    }
    std::copy(other.arr, other.arr + other.sz, arr + sz);
    sz += other.sz;
    return *this;
  }

  String& operator+=(char ch) {
    push_back(ch);
    return *this;
  }

  ~String() {
    delete[] arr;
  }
};

bool operator==(const String& str1, const String& str2) {
  return strcmp(str1.data(), str2.data()) == 0;
}

bool operator!=(const String& str1, const String& str2) {
  return strcmp(str1.data(), str2.data()) != 0;
}

bool operator<(const String& str1, const String& str2) {
  return strcmp(str1.data(), str2.data()) < 0;
}

bool operator>(const String& str1, const String& str2) {
  return str2 < str1;
}

bool operator<=(const String& str1, const String& str2) {
  return !(str2 < str1);
}

bool operator>=(const String& str1, const String& str2) {
  return !(str1 < str2);
}

String operator+(const String& str1, const String& str2) {
  String result = str1;
  result += str2;
  return result;
}

String operator+(const String& str1, char ch) {
  String result = str1;
  result += ch;
  return result;
}

String operator+(char ch, const String& str1) {
  String result(1, ch);
  result += str1;
  return result;
}

std::ostream& operator<<(std::ostream& oss, const String& str) {
  for (size_t i = 0; i < str.size(); ++i) {
    oss << str[i];
  }
  return oss;
}

std::istream& operator>>(std::istream& iss, String& str) {
  str.clear();
  char ch = iss.get();
  while (ch != ' ' && ch != '\n' && !iss.eof()) {
    str.push_back(ch);
    ch = iss.get();
  }
  return iss;
}