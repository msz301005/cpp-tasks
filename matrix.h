#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <array>

class BigInteger;
BigInteger operator+(const BigInteger&, const BigInteger&);
BigInteger operator-(const BigInteger&, const BigInteger&);
BigInteger operator*(const BigInteger&, long long);
BigInteger operator*(long long, const BigInteger&);
BigInteger operator/(const BigInteger&, const BigInteger&);
bool operator==(const BigInteger&, const BigInteger&);
bool operator!=(const BigInteger&, const BigInteger&);
bool operator<(const BigInteger&, const BigInteger&);
bool operator>(const BigInteger&, const BigInteger&);
bool operator<=(const BigInteger&, const BigInteger&);
bool operator>=(const BigInteger&, const BigInteger&);
std::istream& operator>>(std::istream&, BigInteger&);
std::ostream& operator<<(std::ostream&, const BigInteger&);

enum class Sign {
  Negative = -1, Null, Positive
};

Sign abs(const Sign& other) {
  return static_cast<Sign>(std::abs(static_cast<int>(other)));
}

Sign reverseSign(const Sign& other) {
  return static_cast<Sign>(-static_cast<int>(other));
}

Sign operator*(const Sign& first, const Sign& second) {
  return static_cast<Sign>(static_cast<int>(first) * static_cast<int>(second));
}

class BigInteger {
 public:

  BigInteger() = default;

  BigInteger(const BigInteger&) = default;

  BigInteger& operator=(const BigInteger& other) = default;

  BigInteger(int val) {
    BigInteger copy = static_cast<BigInteger>(std::to_string(val));
    *this = copy;
  }

  BigInteger(const std::string& other) {
    popZero();
    auto sz = static_cast<int>(other.size());
    if (std::count(other.begin(), other.end(), '0') == sz
        || (other[0] == '-' &&
            std::count(other.begin(), other.end(), '0') == sz - 1)) {
      sign_ = Sign::Null;
      num = {0};
    } else {
      if (other[0] == '-') {
        sign_ = Sign::Negative;
      }
      for (int i = sz; i > 0; i -= LenBlock) {
        if (i < LenBlock) {
          if (i == 1 && other[0] == '-') {
            continue;
          }
          int fix = (other[0] == '-' ? 1 : 0);
          num.push_back(std::abs(std::stoll(other.substr(fix, i - fix))));
        } else {
          num.push_back(
              std::abs(std::stoll(other.substr(i - LenBlock, LenBlock))));
        }
      }
    }
  }

  explicit operator bool() const { return *this != 0; }

  std::string toString() const;

  Sign getSign() const {
    return sign_;
  }

  void setSign(const Sign& val) {
    sign_ = val;
  }

  size_t size() const {
    return num.size();
  }

  bool getFlagReverse() const {
    return flag_reverse;
  }

  bool empty() const {
    return !num.size();
  }

  long long back() const {
    return num.back();
  }

  long long& operator[](size_t ind) {
    return num[ind];
  }

  const long long& operator[](size_t ind) const {
    return num[ind];
  }

  BigInteger operator-() const {
    BigInteger copy = *this;
    copy.sign_ = reverseSign(copy.sign_);
    return copy;
  }

  BigInteger& operator++() {
    *this += 1;
    return *this;
  }

  BigInteger operator++(int) {
    BigInteger copy = *this;
    *this += 1;
    return copy;
  }

  BigInteger& operator--() {
    *this -= 1;
    return *this;
  }

  BigInteger operator--(int) {
    BigInteger copy = *this;
    *this -= 1;
    return copy;
  }

  BigInteger& operator-=(const BigInteger& other) {
    *this += -other;
    return *this;
  }

  BigInteger& operator+=(const BigInteger& other) {
    if (sign_ == Sign::Null) {
      *this = other;
      return *this;
    }
    if (other.sign_ == Sign::Null) {
      return *this;
    }
    if (sign_ == Sign::Negative && other.sign_ == Sign::Positive) {
      sign_ = Sign::Positive;
      if (other < *this) {
        leftSubtract(*this, other);
        sign_ = Sign::Negative;
      } else {
        rightSubtract(*this, other);
      }
      return *this;
    }
    if (sign_ == Sign::Positive && other.sign_ == Sign::Negative) {
      if (*this < -other) {
        rightSubtract(*this, -other);
        sign_ = Sign::Negative;
      } else {
        leftSubtract(*this, -other);
      }
      return *this;
    }
    int carry = 0;
    for (size_t i = 0; i < std::max(num.size(), other.size()) || carry; ++i) {
      if (i == num.size()) {
        num.push_back(0);
      }
      num[i] += carry + (i < other.size() ? other[i] : 0);
      carry = num[i] >= Base;
      if (carry) {
        num[i] -= Base;
      }
    }
    return *this;
  }

  BigInteger& operator*=(long long other) {
    if (sign_ == Sign::Null) {
      return *this;
    }
    long long carry = 0;
    for (size_t i = 0; i < num.size() || carry; ++i) {
      if (i == num.size()) {
        num.push_back(0);
      }
      long long cur = carry + num[i] * other;
      num[i] = cur % Base;
      carry = cur / Base;
    }
    popZero();
    return *this;
  }

  BigInteger& operator*=(const BigInteger& other) {
    if (sign_ == Sign::Null) {
      return *this;
    }
    BigInteger res;
    res.num.resize(num.size() + other.size());
    res.sign_ = sign_ * other.sign_;
    for (size_t i = 0; i < num.size(); ++i) {
      long long carry = 0;
      for (size_t j = 0; j < other.size() || carry; ++j) {
        long long cur =
            res[i + j] + num[i] * 1ll * (j < other.size() ? other[j] : 0) +
            carry;
        res[i + j] = cur % Base;
        carry = cur / Base;
      }
    }
    res.popZero();
    *this = res;
    return *this;
  }

  BigInteger& operator/=(const BigInteger& other) {
    if (sign_ == Sign::Null) {
      num = {0};
      sign_ = Sign::Null;
      return *this;
    }
    if (other == 1) {
      return *this;
    }
    BigInteger div;
    BigInteger result;
    BigInteger other_copy = other;
    div.flag_reverse = result.flag_reverse = true;
    div.push_back(num.back());
    result.sign_ = sign_ * other.sign_;
    other_copy.sign_ = sign_ = Sign::Positive;
    for (int i = static_cast<int>(num.size()) - 2; i >= 0; --i) {
      if (div == 0) {
        result.push_back(0);
      } else if (div >= other_copy) {
        long long left = -1;
        long long right = Base + 1;
        while (right - left > 1) {
          long long mid = (left + right) / 2;
          (other_copy * mid <= div ? left : right) = mid;
        }
        div -= other_copy * left;
        div.popFrontZero();
        result.push_back(left);
      }
      div.push_back(num[i]);
      div.popFrontZero();
    }
    if (div == 0) {
      result.push_back(0);
    } else if (div >= other_copy) {
      long long left = -1;
      long long right = Base + 1;
      while (right - left > 1) {
        long long mid = (left + right) / 2;
        (other_copy * mid <= div ? left : right) = mid;
      }
      div -= other_copy * left;
      div.popFrontZero();
      result.push_back(left);
    }
    reverse(result.num.begin(), result.num.end());
    result.flag_reverse = false;
    *this = result;
    return *this;
  }

  BigInteger& operator%=(const BigInteger& other) {
    BigInteger div = (*this / other);
    div *= other;
    *this -= div;
    return *this;
  }
 private:
  static const long long Base;
  static const long long LenBlock;
  std::vector<long long> num;
  bool flag_reverse = false;
  Sign sign_ = Sign::Positive;

  static void leftSubtract(BigInteger& first, const BigInteger& second) {
    int carry = 0;
    for (size_t i = 0; i < second.size() || carry; ++i) {
      size_t j = (first.flag_reverse ? first.size() - i - 1 : i);
      first[j] -= carry + (i < second.size() ? second[i] : 0);
      carry = first[j] < 0;
      if (carry) {
        first[j] += Base;
      }
    }
    first.popZero();
  }

  static void rightSubtract(BigInteger& second, const BigInteger& first) {
    BigInteger second_copy(first);
    int carry = 0;
    int second_sz = static_cast<int>(second.size());
    for (size_t i = 0; i < second.size() || carry; ++i) {
      size_t j = (second.flag_reverse ? second_sz - i - 1 : i);
      second_copy[i] -= carry + (i < second.size() ? second[j] : 0);
      carry = second_copy[i] < 0;
      if (carry) {
        second_copy[i] += Base;
      }
    }
    second_copy.popZero();
    second = second_copy;
  }

  void popZero() {
    while (num.size() > 1 && num.back() == 0) {
      num.pop_back();
    }
  }

  void popFrontZero() {
    reverse(num.begin(), num.end());
    popZero();
    reverse(num.begin(), num.end());
  }

  void push_back(long long val) {
    num.push_back(val);
  }

};

const long long BigInteger::Base = 1'000'000'000;
const long long BigInteger::LenBlock = 9;

std::string BigInteger::toString() const {
  std::string result;
  std::stringstream str;
  str << *this;
  str >> result;
  return result;
}

BigInteger operator ""_bi(unsigned long long x) {
  return BigInteger(x);
}

BigInteger operator+(const BigInteger& num1, const BigInteger& num2) {
  BigInteger result = num1;
  result += num2;
  return result;
}

BigInteger operator-(const BigInteger& num1, const BigInteger& num2) {
  BigInteger result = num1;
  result -= num2;
  return result;
}

BigInteger operator*(const BigInteger& num1, const BigInteger& num2) {
  BigInteger result = num1;
  result *= num2;
  return result;
}

BigInteger operator*(const BigInteger& num1, long long num2) {
  BigInteger result = num1;
  result *= num2;
  return result;
}

BigInteger operator*(long long num2, const BigInteger& num1) {
  BigInteger result = num1;
  result *= num2;
  return result;
}

BigInteger operator/(const BigInteger& num1, const BigInteger& num2) {
  BigInteger result = num1;
  result /= num2;
  return result;
}

BigInteger operator%(const BigInteger& num1, const BigInteger& num2) {
  BigInteger result = num1;
  result %= num2;
  return result;
}

bool operator==(const BigInteger& first, const BigInteger& second) {
  bool isNull = true;
  for (size_t i = 0; i < first.size(); i++) {
    if (first[i] != second[i]) {
      return false;
    }
    if (first[i] || second[i]) {
      isNull = false;
    }
  }
  if (isNull) {
    return true;
  }
  return !(first.getSign() != first.getSign() || first.size() != second.size());
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
  return !(first == second);
}

bool operator<(const BigInteger& first, const BigInteger& second) {
  bool flag = !first.getFlagReverse() && !second.getFlagReverse();
  if (flag) {
    if (first.getSign() != second.getSign()) {
      return first.getSign() < second.getSign();
    }
  }
  if (first.size() != second.size() && flag) {
    return (first.size() < second.size()) ^ (first.getSign() == Sign::Negative);
  } else if (first.size() != second.size()) {
    return (first.size() < second.size());
  }

  int first_sz = static_cast<int>(first.size());
  int second_sz = static_cast<int>(second.size());
  for (int i = 0; i < first_sz; i++) {
    int j_first = (!first.getFlagReverse() ? first_sz - i - 1 : i);
    int j_second = (!second.getFlagReverse() ? second_sz - i - 1 : i);
    if (first[j_first] != second[j_second]) {
      if (flag) {
        return (first[j_first] < second[j_second]) ^
               (first.getSign() == Sign::Negative);
      }
      return (first[j_first] < second[j_second]);
    }
  }
  return false;
}

bool operator>(const BigInteger& first, const BigInteger& second) {
  return second < first;
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
  return !(second < first);
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
  return !(first < second);
}

std::istream& operator>>(std::istream& iss, BigInteger& num) {
  std::string str;
  iss >> str;
  BigInteger num1 = static_cast<BigInteger>(str);
  num = num1;
  return iss;
}

std::ostream& operator<<(std::ostream& oss, const BigInteger& num) {
  if (num == 0) {
    oss << '0';
    return oss;
  }
  if (num.getSign() == Sign::Negative) {
    oss << '-';
  }
  oss << (num.empty() ? 0 : num.back());
  auto sz = static_cast<int>(num.size());
  for (int i = sz - 2; i >= 0; --i) {
    int size = static_cast<int>(std::to_string(num[i]).size());
    for (int j = 0; j < 9 - size; j++) {
      oss << '0';
    }
    oss << num[i];
  }
  return oss;
}

BigInteger Gcd(BigInteger first, BigInteger second) {
  if (first.getSign() == Sign::Negative) {
    first.setSign(Sign::Positive);
  }
  if (second.getSign() == Sign::Negative) {
    second.setSign(Sign::Positive);
  }
  BigInteger* first_ptr = &first;
  BigInteger* second_ptr = &second;
  while (*second_ptr) {
    *first_ptr %= *second_ptr;
    std::swap(first_ptr, second_ptr);
  }
  return *first_ptr;
}

class Rational {
 public:

  const BigInteger& getNum() const { return num; }

  const BigInteger& getDen() const { return den; }

  Rational() = default;

  Rational(const Rational&) = default;

  Rational& operator=(const Rational&) = default;

  Rational(const BigInteger& other) : num(other), den(1) {}

  Rational(int val) : num(val), den(1) {};

  Rational& operator=(int val) {
    num = val;
    den = 1;
    return *this;
  }

  Rational& operator=(const BigInteger& other) {
    num = other;
    den = 1;
    return *this;
  }

  std::string toString() {
    reduction();
    std::string res = num.toString();
    if (num != 0 && den != 1) {
      res += '/';
      res += den.toString();
    }
    return res;
  }

  std::string asDecimal(size_t precision = 0) const {
    std::string res;
    if (num.getSign() == Sign::Negative) {
      res += '-';
    }
    BigInteger int_part = num / den;
    res += (int_part.getSign() == Sign::Positive ? int_part.toString()
                                                 : (-int_part).toString());
    if (precision > 0) {
      BigInteger factor = 1;
      for (size_t i = 0; i < precision; ++i) {
        factor *= 10;
      }
      auto val = num * factor / den;
      std::string part = (val.getSign() == Sign::Positive ? val.toString()
                                                          : (-val).toString());
      if (int_part != 0) {
        part = part.substr(res.size());
      }
      res += '.';
      for (int i = 0;
           i < static_cast<int>(precision) - static_cast<int>(part.size());
           ++i) {
        res += '0';
      }
      res += part;
    }
    return res;
  }

  explicit operator double() const {
    std::stringstream str;
    double res;
    str << asDecimal(300);
    str >> res;
    return res;
  }

  Rational operator-() const {
    Rational copy = *this;
    copy.setSign(reverseSign(num.getSign()));
    return copy;
  }

  Rational& operator+=(const Rational& other) {
    BigInteger gcd = Gcd(den, other.den);
    BigInteger lcm = den * other.den / gcd;
    num *= lcm / den;
    num += other.num * lcm / other.den;
    den = lcm;
    reduction();
    return *this;
  }

  Rational& operator-=(const Rational& other) {
    BigInteger gcd = Gcd(den, other.den);
    BigInteger lcm = den * other.den / gcd;
    num *= lcm / den;
    num -= other.num * lcm / other.den;
    den = lcm;
    reduction();
    return *this;
  }

  Rational& operator*=(const Rational& other) {
    num *= other.num;
    den *= other.den;
    reduction();
    return *this;
  }

  Rational& operator/=(const Rational& other) {
    Sign sign = num.getSign() * other.num.getSign();
    num *= other.den;
    den *= other.num;
    setSign(sign);
    reduction();
    return *this;
  }

 private:
  void reduction() {
    if (num != 0) {
      BigInteger gcd = Gcd(num, den);
      num /= gcd;
      den /= gcd;
    }
  }

  void setSign(const Sign& val) {
    num.setSign(val);
    den.setSign(Sign::Positive);
  }

  BigInteger num;
  BigInteger den;
};

bool operator==(const Rational& first, const Rational& second) {
  return first.getNum() * second.getDen() == first.getDen() * second.getNum();
}

bool operator!=(const Rational& first, const Rational& second) {
  return !(first == second);
}

bool operator<(const Rational& first, const Rational& second) {
  return first.getNum() * second.getDen() < first.getDen() * second.getNum();
}

bool operator>(const Rational& left, const Rational& right) {
  return right < left;
}

bool operator<=(const Rational& left, const Rational& right) {
  return !(right < left);
}

bool operator>=(const Rational& left, const Rational& right) {
  return !(left < right);
}

Rational operator+(const Rational& rat1, const Rational& rat2) {
  Rational result = rat1;
  result += rat2;
  return result;
}

Rational operator-(const Rational& rat1, const Rational& rat2) {
  Rational result = rat1;
  result -= rat2;
  return result;
}

Rational operator*(const Rational& rat1, const Rational& rat2) {
  Rational result = rat1;
  result *= rat2;
  return result;
}

Rational operator/(const Rational& rat1, const Rational& rat2) {
  Rational result = rat1;
  result /= rat2;
  return result;
}

std::istream& operator>>(std::istream& in, Rational& first) {
  std::string str;
  in >> str;
  first = Rational(str);
  return in;
}

std::ostream& operator<<(std::ostream& out, const Rational& first) {
  out << first.getNum() << '/' << first.getDen();
  return out;
}

// Next comes Residue and Matrix code

template <size_t N, size_t l = 1, size_t r = N>
struct ValueSqrt {
  const static size_t mid = (l + r) / 2;
  const static size_t value = ValueSqrt<N,
                                        (mid * mid >= N ? l : mid + 1),
                                        (mid * mid >= N ? mid : r)>::value;
};

template <size_t N, size_t value_>
struct ValueSqrt<N, value_, value_> {
  static const size_t value = value_;
};

template <size_t N>
const static size_t Sqrt = ValueSqrt<N>::value;

template <size_t N, size_t M>
struct IsPrime {
  const static bool value = (N % M != 0) && IsPrime<N, M - 1>::value;
};

template <size_t N>
struct IsPrime<N, 1> {
  const static bool value = true;
};

template <size_t N>
struct IsPrime<N, 0> {
  const static bool value = true;
};

template <size_t N>
const static bool
    is_prime = IsPrime<N, Sqrt<N> - (Sqrt<N> * Sqrt<N> > N)>::value;

template <size_t N>
class Residue {
 public:
  Residue() : remainder(0) {}
  Residue(int n, int N_ = N) : remainder((N_ + (n % N_)) % N_) {}
  explicit operator int() { return remainder; }
  size_t rem() const { return remainder; }

  Residue<N>& operator+=(const Residue<N>& other) {
    *this = Residue<N>((remainder + other.remainder) % N);
    return *this;
  }

  Residue<N>& operator-=(const Residue<N>& other) {
    *this = Residue<N>(
        static_cast<int>(remainder) - static_cast<int>(other.remainder));
    return *this;
  }

  Residue<N>& operator*=(const Residue<N>& other) {
    *this = Residue<N>((remainder * other.remainder) % N);
    return *this;
  }

  Residue<N>& operator/=(const Residue<N>& other) {
    static_assert(is_prime<N>);
    for (size_t i = 0; i < N; ++i) {
      if ((other.remainder * i) % N == remainder) {
        *this = Residue<N>(i);
        return *this;
      }
    }
    *this = Residue<N>(0);
    return *this;
  }

 private:
  size_t remainder;
};

template <size_t N>
Residue<N> operator+(const Residue<N>& first, const Residue<N>& second) {
  Residue<N> result = first;
  result += second;
  return result;
}

template <size_t N>
Residue<N> operator-(const Residue<N>& first, const Residue<N>& second) {
  Residue<N> result = first;
  result -= second;
  return result;
}

template <size_t N>
Residue<N> operator*(const Residue<N>& first, const Residue<N>& second) {
  Residue<N> result = first;
  result *= second;
  return result;
}

template <size_t N>
Residue<N> operator/(const Residue<N>& first, const Residue<N>& second) {
  Residue<N> result = first;
  result /= second;
  return result;
}

template <size_t N>
bool operator==(const Residue<N>& first, const Residue<N>& second) {
  return first.rem() == second.rem();
}

template <size_t N>
std::istream& operator>>(std::istream& in, Residue<N>& first) {
  int val;
  in >> val;
  first = Residue<N>(val);
  return in;
}

template <size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& first) {
  out << first.rem();
  return out;
}

template <size_t M, size_t N, typename Field=Rational>
class Matrix {
  std::vector<std::vector<Field>> matrix;

  std::pair<size_t, size_t> nonZeroElem(size_t pos_i, size_t pos_j) const {
    for (size_t j = pos_j; j < N; ++j) {
      for (size_t i = pos_i; i < M; ++i) {
        if (matrix[i][j] != Field(0)) {
          return {i, j};
        }
      }
    }
    return {M + 1, N + 1};
  }

  std::pair<Matrix<M, N, Field>, Field> forwardGaussianMove(Matrix<M,
                                                                   N,
                                                                   Field>& unity_matrix) const {
    auto copy_matrix = *this;
    Field coefficient = 1;
    size_t i = 0, j = 0;
    for (; i < M && j < N; ++i, ++j) {
      auto non_zero_elem = copy_matrix.nonZeroElem(i, j);
      if (non_zero_elem.first == M + 1 && non_zero_elem.second == N + 1) {
        return {copy_matrix, coefficient};
      }
      if (i != non_zero_elem.first) {
        copy_matrix.swapRows(i, non_zero_elem.first);
        unity_matrix.swapRows(i, non_zero_elem.first);
        coefficient *= Field(-1);
      }
      j = non_zero_elem.second;
      auto factor = copy_matrix[i][j];
      coefficient *= factor;
      copy_matrix.divisionRow(i, factor);
      unity_matrix.divisionRow(i, factor);
      for (size_t down_i = i + 1; down_i < M; ++down_i) {
        if (copy_matrix[down_i][j] == Field(0)) {
          continue;
        }
        factor = copy_matrix[down_i][j];
        for (size_t j1 = 0; j1 < N; ++j1) {
          copy_matrix[down_i][j1] -= factor * copy_matrix[i][j1];
          unity_matrix[down_i][j1] -= factor * unity_matrix[i][j1];
        }
      }
    }
    return {copy_matrix, coefficient};
  }

  std::pair<Matrix<M, N, Field>, Field> forwardGaussianMove() const {
    Matrix<M, N, Field> trash;
    return forwardGaussianMove(trash);
  }

  Matrix<M, N, Field> reverseGaussianMove(Matrix<M,
                                                 N,
                                                 Field>& unity_matrix) const {
    auto copy_matrix = *this;
    int i = static_cast<int>(M) - 1;
    int j = static_cast<int>(N) - 1;
    for (; i >= 0 && copy_matrix[i][j] != Field(1); --i);
    for (; i >= 0 && j >= 0; --i, --j) {
      for (; j >= 0 && copy_matrix[i][j] != Field(1); --j);
      for (int up_i = i - 1; up_i >= 0; --up_i) {
        if (copy_matrix[up_i][j] == Field(0)) {
          continue;
        }
        auto factor = copy_matrix[up_i][j];
        for (int j1 = 0; j1 < static_cast<int>(N); ++j1) {
          copy_matrix[up_i][j1] -= factor * copy_matrix[i][j1];
          unity_matrix[up_i][j1] -= factor * unity_matrix[i][j1];
        }
      }
    }
    return copy_matrix;
  }

  void divisionRow(size_t number, const Field& value) {
    if (value == Field(0)) {
      return;
    }
    for (auto& elem : matrix[number]) {
      elem /= value;
    }
  }

  void swapRows(size_t first, size_t second) {
    std::swap(matrix[first], matrix[second]);
  }

  static std::vector<Field> nullCol() {
    return std::vector<Field>(N, Field(0));
  }

  static Field detTriangularMatrix(const Matrix<M, N, Field>& matrix) {
    Field result = 1;
    for (size_t diag = 0; diag < M; ++diag) {
      result *= matrix[diag][diag];
    }
    return result;
  }

 public:
  Matrix() { matrix.resize(M, std::vector<Field>(N, Field(0))); }

  Matrix(std::initializer_list<std::initializer_list<Field>> matrix_) {
    for (auto& row : matrix_) {
      std::vector<Field> v;
      for (auto& elem : row) {
        v.push_back(elem);
      }
      matrix.push_back(v);
    }
  }

  static Matrix<M, N, Field> unityMatrix() {
    static_assert(M == N);
    Matrix<M, N, Field> unity_matrix;
    for (size_t i = 0; i < M; ++i) {
      unity_matrix[i][i] = Field(1);
    }
    return unity_matrix;
  }

  Matrix<N, M, Field> transposed() const {
    Matrix<N, M, Field> transposed_matrix;
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        transposed_matrix[j][i] = matrix[i][j];
      }
    }
    return transposed_matrix;
  }

  Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& other) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        matrix[i][j] += other.matrix[i][j];
      }
    }
    return *this;
  }

  Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& other) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        matrix[i][j] -= other.matrix[i][j];
      }
    }
    return *this;
  }

  Matrix<M, N, Field>& operator*=(const Field& value) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        matrix[i][j] *= value;
      }
    }
    return *this;
  }

  template <size_t K>
  Matrix<M, K, Field>& operator*=(const Matrix<N, K, Field> other) {
    static_assert(N == K);
    *this = *this * other;
    return *this;
  }

  std::vector<Field>& operator[](size_t i) { return matrix[i]; }

  const std::vector<Field>& operator[](size_t i) const { return matrix[i]; }

  Field det() const {
    static_assert(M == N);
    auto [copy_matrix, coefficient] = forwardGaussianMove();
    return detTriangularMatrix(copy_matrix) * coefficient;
  }

  size_t rank() const {
    Matrix<M, N, Field> triangular_matrix = forwardGaussianMove().first;
    size_t result = 0;
    auto null_col = nullCol();
    for (size_t i = 0; i < M; ++i) {
      result += (triangular_matrix[i] != null_col);
    }
    return result;
  }

  Field trace() const {
    static_assert(N == M);
    Field result = 0;
    for (size_t diag = 0; diag < M; ++diag) {
      result += matrix[diag][diag];
    }
    return result;
  }

  std::array<Field, N> getRow(size_t i) const {
    std::array<Field, N> row;
    for (size_t j = 0; j < N; ++j) {
      row[j] = matrix[i][j];
    }
    return row;
  }

  std::array<Field, M> getColumn(size_t j) const {
    std::array<Field, M> col;
    for (size_t i = 0; i < M; ++i) {
      col[i] = matrix[i][j];
    }
    return col;
  }

  void invert() {
    static_assert(M == N);
    Matrix<M, N, Field> unity_matrix = unityMatrix();
    *this = forwardGaussianMove(unity_matrix).first;
    *this = reverseGaussianMove(unity_matrix);
    std::cerr << *this << std::endl;
    *this = unity_matrix;
  }

  Matrix<M, N, Field> inverted() {
    static_assert(M == N);
    Matrix<M, N, Field> invert_matrix = *this;
    invert_matrix.invert();
    return invert_matrix;
  }

};

template <size_t M, size_t N, typename Field>
bool operator==(const Matrix<M, N, Field>& first,
                const Matrix<M, N, Field>& second) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (first[i][j] != second[i][j]) {
        return false;
      }
    }
  }
  return true;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field>& first,
                              const Matrix<M, N, Field>& second) {
  Matrix<M, N, Field> result = first;
  result += second;
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& first,
                              const Matrix<M, N, Field>& second) {
  Matrix<M, N, Field> result = first;
  result -= second;
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Matrix<M,
                                           N,
                                           Field>& first, const Field& value) {
  Matrix<M, N, Field> result = first;
  result *= value;
  return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Field& value, const Matrix<M,
                                                               N,
                                                               Field>& first) {
  Matrix<M, N, Field> result = first;
  result *= value;
  return result;
}

template <size_t M, size_t K, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Matrix<M, K, Field>& first,
                              const Matrix<K, N, Field>& second) {
  Matrix<M, N, Field> matrix;
  for (size_t i_m = 0; i_m < M; ++i_m) {
    for (size_t i_n = 0; i_n < N; ++i_n) {
      Field sum = 0;
      for (size_t i_k = 0; i_k < K; ++i_k) {
        sum += first[i_m][i_k] * second[i_k][i_n];
      }
      matrix[i_m][i_n] = sum;
    }
  }
  return matrix;
}

template <size_t N, typename Field=Rational>
using SquareMatrix = Matrix<N, N, Field>;

template <size_t M, size_t N, typename Field>
std::ostream& operator<<(std::ostream& out, const Matrix<M, N, Field>& matrix) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      out << matrix[i][j] << ' ';
    }
    out << '\n';
  }
  return out;
}

template <size_t M, size_t N, typename Field>
std::istream& operator>>(std::istream& in, Matrix<M, N, Field>& matrix) {
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      in >> matrix[i][j];
    }
  }
  return in;
}
