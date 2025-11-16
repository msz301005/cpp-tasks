#include <bits/stdc++.h>

namespace {
bool isEqual(long double x, long double y) {
  long double eps = 1e-5;
  return std::abs(x - y) < eps;
}
}

struct Point;
class Line;

class Shape {
 public:
  virtual long double perimeter() const = 0;
  virtual long double area() const = 0;
  virtual ~Shape() = default;

  bool operator==(const Shape& another) const;

  virtual bool isCongruentTo(const Shape& another) const = 0;

  virtual bool isSimilarTo(const Shape& another) const = 0;

  virtual bool containsPoint(const Point& point) const = 0;

  virtual void rotate(const Point& center, long double angle) = 0;

  virtual void reflect(const Point& center) = 0;

  virtual void reflect(const Line& axis) = 0;

  virtual void scale(const Point& center, long double coefficient) = 0;
};

struct Point {
  long double x;
  long double y;
  Point() = default;
  Point(long double x, long double y) : x(x), y(y) {}

  virtual bool operator==(const Point& other) const {
    return isEqual(x, other.x) && isEqual(y, other.y);
  }

  virtual bool operator!=(const Point& other) const {
    return !(*this == other);
  }

  Point& operator+=(const Point& other) {
    x += other.x;
    y += other.y;
    return *this;
  }

  Point& operator/=(long double other) {
    x /= other;
    y /= other;
    return *this;
  }

  Point& operator-=(const Point& other) {
    x -= other.x;
    y -= other.y;
    return *this;
  }

  void rotate(const Point& center, long double cos, long double sin) {
    auto copy_x = x;
    x = center.x + (x - center.x) * cos - (y - center.y) * sin;
    y = center.y + (copy_x - center.x) * sin + (y - center.y) * cos;
  }

  void rotate(const Point& center, long double angle) {
    angle *= M_PI / 180;
    auto copy_x = x;
    x = center.x + (x - center.x) * std::cos(angle) -
        (y - center.y) * std::sin(angle);
    y = center.y + (copy_x - center.x) * std::sin(angle) +
        (y - center.y) * std::cos(angle);
  }

  void reflect(const Point& center);

  void reflect(const Line& axis);

  void scale(const Point& center, long double coefficient);

};

Point operator+(const Point& first, const Point& second) {
  Point result = first;
  result += second;
  return result;
}

Point operator-(const Point& first, const Point& second) {
  Point result = first;
  result -= second;
  return result;
}

Point operator/(const Point& first, long double second) {
  Point result = first;
  result /= second;
  return result;
}

Point operator*(long double first, const Point& second) {
  Point result = second;
  result.x *= first;
  result.y *= first;
  return result;
}

class Vector : public Point {
 public:
  using Point::Point;
  Vector() = default;
  Vector(const Point& p1, const Point& p2) {
    x = p2.x - p1.x;
    y = p2.y - p1.y;
  }

  long double len() const { return std::sqrt(x * x + y * y); }

  long double operator*(const Vector& other) {
    return x * other.y - other.x * y;
  }
};

class Line {
 public:
  Line(long double k, long double b) : A_(-k), B_(1.0), C_(-b) {}
  Line(const Point& p, long double k) : A_(-k), B_(1.0), C_(k * p.x - p.y) {}
  Line(const Point& p1, const Point& p2)
      : A_(p2.y - p1.y)
      , B_(p1.x - p2.x)
      , C_(p1.y * p2.x - p1.x * p2.y) {}
  bool operator==(const Line& other) const {
    return isEqual(A_ * other.B_, other.A_ * B_) &&
           isEqual(A_ * other.C_, other.A_ * C_);
  }

  bool operator!=(const Line& other) const {
    return !(*this == other);
  }

  long double GetA() const { return A_; }

  long double GetB() const { return B_; }

  long double GetC() const { return C_; }

  void SetA(long double A) { A_ = A; }

  void SetB(long double B) { B_ = B; }

  void SetC(long double C) { C_ = C; }

 private:
  long double A_, B_, C_;
};

class Polygon : public Shape {
 public:

  Polygon(const std::vector<Point>& ver) : ver_(ver) {}

  template <class... Args>
  Polygon(Args... args) : ver_({args...}) {}

  size_t verticesCount() const { return ver_.size(); }

  const std::vector<Point>& getVertices() const { return ver_; }

  bool isConvex() const {
    bool flag_positive = false;
    bool flag_negative = false;
    for (size_t i = 0; i < ver_.size() - 2; ++i) {
      Vector v1(ver_[i], ver_[(i + 1) % ver_.size()]);
      Vector v2(ver_[(i + 1) % ver_.size()], ver_[(i + 2) % ver_.size()]);
      (v1 * v2 >= 0 ? flag_positive : flag_negative) = true;
    }
    return flag_positive ^ flag_negative;
  }

  long double area() const override {
    long double res = 0;
    int size = static_cast<int>(ver_.size());
    for (int i = 1; i < size - 1; ++i) {
      res += Vector(ver_[0], ver_[i]) * Vector(ver_[0], ver_[i + 1]);
    }
    return std::abs(res / 2.0);
  }

  long double perimeter() const override {
    long double res = 0;
    int size = static_cast<int>(ver_.size());
    for (int i = 0; i < size - 1; ++i) {
      res += Vector(ver_[i], ver_[i + 1]).len();
    }
    return res;
  }

//  bool operator==(const Shape& another) const override;
  bool IsEqual(const Shape& another) const;
  bool isCongruentTo(const Shape& another) const override;

  bool isSimilarTo(const Shape& another) const override;

  bool containsPoint(const Point& point) const override;

  void rotate(const Point& center, long double angle) override;

  void reflect(const Point& center) override;

  void reflect(const Line& axis) override;

  void scale(const Point& center, long double coefficient) override;

 protected:
  std::vector<Point> ver_;
};

class Ellipse : public Shape {
 public:
  friend bool IsEqual(const Ellipse&, const Ellipse&);

  Ellipse(Point f1, Point f2, long double s)
      : f1_(std::move(f1))
      , f2_(std::move(f2))
      , center_(Point((f1_.x + f2_.x) / 2, (f1_.y + f2_.y) / 2))
      , a_(s / 2)
      , c_(Vector(center_, f1_).len())
      , e_(c_ / a_) { b_ = std::sqrt(a_ * a_ - c_ * c_); }

  std::pair<Point, Point> focuses() { return {f1_, f2_}; }

  std::pair<Line, Line> directrices() const {
    Vector v1(center_, f1_);
    Vector v2(center_, f2_);
    Point p_dir1(center_.x + v1.x / v1.len() * a_ / e_,
                 center_.y + v1.y / v1.len() * a_ / e_);
    Point p_dir2(center_.x + v2.x / v2.len() * a_ / e_,
                 center_.y + v2.y / v2.len() * a_ / e_);
    Line f1_f2(f1_, f2_);
    Line dir1 = f1_f2;
    Line dir2 = f1_f2;
    dir1.SetA(-f1_f2.GetB());
    dir2.SetA(-f1_f2.GetB());
    dir1.SetB(f1_f2.GetA());
    dir2.SetB(f1_f2.GetA());
    dir1.SetC(dir1.GetB() * p_dir1.x - dir1.GetA() + p_dir1.y);
    dir2.SetC(dir2.GetB() * p_dir2.x - dir2.GetA() + p_dir2.y);
    return {dir1, dir2};
  }

  long double eccentricity() const { return e_; }

  Point center() const { return center_; }

  long double area() const override { return M_PI * a_ * b_; }

  long double perimeter() const override {
    return M_PI * (a_ + b_) * (1 + (3 * std::pow((a_ - b_) / (a_ + b_), 2)) /
                                   (10 + std::sqrt(4 - (3 * std::pow(
                                       (a_ - b_) / (a_ + b_), 2)))));
  }

  // bool operator==(const Shape& another) const override;
  bool IsEqual(const Shape& another) const;
  bool isCongruentTo(const Shape& another) const override;

  bool isSimilarTo(const Shape& another) const override;

  bool containsPoint(const Point& point) const override;

  void rotate(const Point& center, long double angle) override;

  void reflect(const Point& center) override;

  void reflect(const Line& axis) override;

  void scale(const Point& center, long double coefficient) override;

 protected:
  Point f1_, f2_, center_;
  long double a_, b_, c_, e_;
};

class Circle : public Ellipse {
 public:
  Circle(const Point& center, long double radius)
      :
      Ellipse(center, center, 2 * radius) {}

  long double radius() const { return a_; }

  long double area() const override { return M_PI * a_ * a_; }

  long double perimeter() const override { return 2 * M_PI * a_; }

//  bool operator==(const Shape& another) const override;
  bool IsEqual(const Shape& another) const;
  bool isCongruentTo(const Shape& another) const override;

  bool isSimilarTo(const Shape& another) const override;

  bool containsPoint(const Point& point) const override;

  void rotate(const Point& center, long double angle) override;

  void reflect(const Point& center) override;

  void reflect(const Line& axis) override;

  void scale(const Point& center, long double coefficient) override;
};

class Rectangle : public Polygon {
 public:
  Rectangle(const std::vector<Point>& ver)
      : Polygon(ver)
      , p1_(ver[0])
      , p2_(ver[1])
      , p3_(ver[2])
      , p4_(ver[3]) {}

  Rectangle(Point p1, Point p2, Point p3, Point p4)
      : Polygon(p1, p2, p3, p4)
      , p1_(std::move(p1))
      , p2_(std::move(p2))
      , p3_(std::move(p3))
      , p4_(std::move(p4)) {}

  Rectangle(const Point& p1, const Point& p3, long double k)
      : p1_(p1)
      , p2_(p1)
      , p3_(p3)
      , p4_(p3) {
    if (k < 1) {
      k = 1 / k;
    }
    long double cos_2 = (1 - k * k) / (1 + k * k);
    long double sin_2 = 2 * k / (1 + k * k);
    Point center = (p1 + p3) / 2;
    p2_.rotate(center, -cos_2, -sin_2);
    p4_.rotate(center, -cos_2, -sin_2);
    ver_ = {p1_, p2_, p3_, p4_};
  }

  Point center() const { return (p1_ + p3_) / 2; }

  std::pair<Line, Line> diagonals() const {
    return {Line(p1_, p3_), Line(p2_, p4_)};
  }

  long double area() const override {
    return Vector(p1_, p2_).len() * Vector(p3_, p4_).len();
  }

  long double perimeter() const override {
    return 2 * (Vector(p1_, p2_).len() + Vector(p3_, p4_).len());
  }

//  bool operator==(const Shape& another) const override;
  bool IsEqual(const Shape& another) const;
  bool isCongruentTo(const Shape& another) const override;

  bool isSimilarTo(const Shape& another) const override;

  bool containsPoint(const Point& point) const override;

  void rotate(const Point& center, long double angle) override;

  void reflect(const Point& center) override;

  void reflect(const Line& axis) override;

  void scale(const Point& center, long double coefficient) override;

 protected:
  Point p1_, p2_, p3_, p4_;
};

class Square : public Rectangle {
 public:
  using Rectangle::Rectangle;
  Square(const Point&, const Point&, double) = delete;
  Square(const Point& p1, const Point& p2) : Rectangle(p1, p2, 1.0) {}

  Circle circumscribedCircle() const {
    auto center_ = center();
    return {center_, Vector(p1_, center_).len()};
  }

  Circle inscribedCircle() const {
    auto center_ = center();
    return {center_, Vector(p1_, p2_).len() / 2};
  }

  long double area() const override {
    return std::pow(Vector(p1_, p2_).len(), 2);
  }

  long double perimeter() const override { return 4 * Vector(p1_, p2_).len(); }

//  bool operator==(const Shape& another) const override;
  bool IsEqual(const Shape& another) const;

  bool isCongruentTo(const Shape& another) const override;

  bool isSimilarTo(const Shape& another) const override;

  bool containsPoint(const Point& point) const override;

  void rotate(const Point& center, long double angle) override;

  void reflect(const Point& center) override;

  void reflect(const Line& axis) override;

  void scale(const Point& center, long double coefficient) override;

};

class Triangle : public Polygon {
 public:
  Triangle(const std::vector<Point>& ver)
      : Polygon(ver)
      , p1_(ver[0])
      , p2_(ver[1])
      , p3_(ver[2]) {}

  Triangle(Point p1, Point p2, Point p3)
      : Polygon(p1, p2, p3)
      , p1_(std::move(p1))
      , p2_(std::move(p2))
      , p3_(std::move(p3)) {}

  Circle circumscribedCircle() const {
    long double det = 4 * ((p2_.x - p1_.x) * (p3_.y - p1_.y) -
                           (p3_.x - p1_.x) * (p2_.y - p1_.y));
    long double
        d2 = p2_.x * p2_.x + p2_.y * p2_.y - p1_.x * p1_.x - p1_.y * p1_.y;
    long double
        d3 = p3_.x * p3_.x + p3_.y * p3_.y - p1_.x * p1_.x - p1_.y * p1_.y;
    Point center;
    center.x = (2 * d2 * (p3_.y - p1_.y) - 2 * d3 * (p2_.y - p1_.y)) / det;
    center.y = (2 * d3 * (p2_.x - p1_.x) - 2 * d2 * (p3_.x - p1_.x)) / det;
    long double radius = std::sqrt(
        std::pow(p1_.x - center.x, 2) + std::pow(p1_.y - center.y, 2));
    return {center, radius};
  }

  Circle inscribedCircle() const {
    long double r1r2 =
        std::sqrt(std::pow(p1_.x - p2_.x, 2) + std::pow(p1_.y - p2_.y, 2));
    long double r2r3 =
        std::sqrt(std::pow(p2_.x - p3_.x, 2) + std::pow(p2_.y - p3_.y, 2));
    long double r3r1 =
        std::sqrt(std::pow(p3_.x - p1_.x, 2) + std::pow(p3_.y - p1_.y, 2));
    Point center;
    center.x =
        (r1r2 * p3_.x + r2r3 * p1_.x + r3r1 * p2_.x) / (r1r2 + r2r3 + r3r1);
    center.y =
        (r1r2 * p3_.y + r2r3 * p1_.y + r3r1 * p2_.y) / (r1r2 + r2r3 + r3r1);
    Line side(p1_, p2_);
    long double radius = std::abs(
        side.GetA() * center.x + side.GetB() * center.y + side.GetC()) /
                         std::sqrt(side.GetA() * side.GetA() +
                                   side.GetB() * side.GetB());
    return {center, radius};
  }

  Point centroid() const {
    Vector r1(p1_, p2_);
    Vector r2(p1_, p3_);
    return p1_ + (r1 + r2) / 3;
  }

  Point orthocenter() const {
    Line p1p3(p1_, p3_);
    Line p2p3(p2_, p3_);
    long double A1 = p1p3.GetA();
    long double B1 = p1p3.GetB();
    long double A2 = p2p3.GetA();
    long double B2 = p2p3.GetB();
    long double det = -B1 * A2 + B2 * A1;
    Point center;
    center.x = ((p2_.x * B1 - p2_.y * A1) * (-A2) -
                (p1_.x * B2 - p1_.y * A2) * (-A1)) / det;
    center.y =
        ((p1_.x * B2 - p1_.y * A2) * B1 - (p2_.x * B1 - p2_.y * A1) * B2) / det;
    return center;
  }

  Line EulerLine() const {
    return {circumscribedCircle().center(), orthocenter()};
  }

  Circle ninePointsCircle() const {
    return Triangle((p1_ + p2_) / 2, (p1_ + p3_) / 2,
                    (p2_ + p3_) / 2).circumscribedCircle();
  }

  long double area() const override {
    return std::abs(Vector(p1_, p2_) * Vector(p1_, p3_)) / 2.0;
  }

  long double perimeter() const override {
    return Vector(p1_, p2_).len() + Vector(p1_, p3_).len() +
           Vector(p2_, p3_).len();
  }

//  bool operator==(const Shape& another) const override;
  bool IsEqual(const Shape& another) const;

  bool isCongruentTo(const Shape& another) const override;

  bool isSimilarTo(const Shape& another) const override;

  bool containsPoint(const Point& point) const override;

  void rotate(const Point& center, long double angle) override;

  void reflect(const Point& center) override;

  void reflect(const Line& axis) override;

  void scale(const Point& center, long double coefficient) override;

 protected:
  Point p1_, p2_, p3_;
};

void Point::reflect(const Line& axis) {
  long double coefficient = -(axis.GetA() * x + axis.GetB() * y + axis.GetC()) /
                            (axis.GetA() * axis.GetA() +
                             axis.GetB() * axis.GetB());
  x += 2 * axis.GetA() * coefficient;
  y += 2 * axis.GetB() * coefficient;
}

void Point::reflect(const Point& center) {
  x = 2 * center.x - x;
  y = 2 * center.y - y;
}

void Point::scale(const Point& center, long double coefficient) {
  x = coefficient * (x - center.x) + center.x;
  y = coefficient * (y - center.y) + center.y;
}

bool IsEqual(const std::vector<Point>& first, const std::vector<Point>& second) {
  //std::cerr << "first_polygon\n";
  /*for(auto& elem : first) {
    std::cerr << elem.x << ' ' << elem.y << std::endl;
  }*/
  /*std::cerr << "second_polygon\n";
  for(auto& elem : second) {
    std::cerr << elem.x << ' ' << elem.y << std::endl;
  }*/
  int size_first = static_cast<int>(first.size());
  size_t size_second = second.size();
  if (first.size() != size_second) {
    //std::cerr << "res: 0\n";
    return false;
  }
  int pos = 0;
  bool flag1 = false;
  bool flag2 = true;
  for (size_t i = 0; i < size_second; ++i) {
    if (first[0] == second[i]) {
      flag1 = true;
      pos = static_cast<int>(i);
      break;
    }
  }
  if (!flag1) {
    //std::cerr << "res: 0\n";
    return false;
  }
  int copy_pos = pos + 1;
  for (size_t i = 1; i < first.size();
       ++i, copy_pos = (copy_pos + 1) % size_first) {
    if (first[i] != second[copy_pos]) {
      flag1 = false;
      break;
    }
  }
  copy_pos = pos - 1;
  for (size_t i = 1; i < first.size();
       ++i, copy_pos = (copy_pos <= 0 ? size_first - 1 : copy_pos - 1)) {
    if (first[i] != second[copy_pos]) {
      flag2 = false;
      break;
    }
  }
  //std::cerr << "res: " << (flag1 || flag2) << std::endl;
  return flag1 || flag2;
}

bool IsEqual(const std::vector<long double>& first, const std::vector<long double>& second,
             long double coefficient = 1.0) {
  int size_first = static_cast<int>(first.size());
  size_t size_second = second.size();
  if (first.size() != size_second) {
    return false;
  }
  std::vector<int> array_pos;
  bool flag1 = false;
  bool flag2 = true;
  for (size_t i = 0; i < size_second; ++i) {
    if (isEqual(coefficient * first[0], second[i])) {
      flag1 = true;
      array_pos.push_back(static_cast<int>(i));
    }
  }
  if (!flag1) {
    return false;
  }
  for (auto pos : array_pos) {
    flag1 = flag2 = true;
    int copy_pos = pos + 1;
    for (size_t i = 1; i < first.size();
         ++i, copy_pos = (copy_pos + 1) % size_first) {
      if (!isEqual(coefficient * first[i], second[copy_pos])) {
        flag1 = false;
        break;
      }
    }
    copy_pos = pos - 1;
    for (size_t i = 1; i < first.size();
         ++i, copy_pos = (copy_pos <= 0 ? size_first - 1 : copy_pos - 1)) {
      if (!isEqual(coefficient * first[i], second[copy_pos])) {
        flag2 = false;
        break;
      }
    }
    if (flag1 || flag2) {
      return true;
    }
  }
  return false;
}

bool IsEqual(const Ellipse& first, const Ellipse& second) {
  bool focus = std::make_pair(first.f1_, first.f2_) ==
               std::make_pair(second.f1_, second.f2_)
               ||
               std::make_pair(first.f2_, first.f1_) ==
               std::make_pair(second.f1_, second.f2_);
  return focus && first.a_ == second.a_ && first.b_ == second.b_;
}

bool Polygon::IsEqual(const Shape& other) const {
  const auto* ellipse = dynamic_cast<const Ellipse*>(&other);
  if (ellipse != nullptr) {
    return false;
  }
  const auto* polygon = dynamic_cast<const Polygon*>(&other);
  /*std::cerr << "real second polygon\n";
  for(auto& elem : polygon->ver_) {
    std::cerr << elem.x << ' ' << elem.y << std::endl;
  }*/
  return ::IsEqual(ver_, polygon->ver_);
}

bool Ellipse::IsEqual(const Shape& other) const {
  const auto* polygon = dynamic_cast<const Polygon*>(&other);
  if (polygon != nullptr) {
    return false;
  }
  const auto* ellipse = dynamic_cast<const Ellipse*>(&other);
  return ::IsEqual(*this, *ellipse);
}

bool Circle::IsEqual(const Shape& other) const {
  return this->Ellipse::IsEqual(other);
}

bool Rectangle::IsEqual(const Shape& other) const {
  return this->Polygon::IsEqual(other);
}

bool Square::IsEqual(const Shape& other) const {
  return this->Polygon::IsEqual(other);
}

bool Triangle::IsEqual(const Shape& other) const {
  return this->Polygon::IsEqual(other);
}

bool Shape::operator==(const Shape& another) const {
  const auto* polygon = dynamic_cast<const Polygon*>(&another);
  if (polygon != nullptr) {
    /*std::cerr << "real first polygon\n";
    for(auto& elem : polygon->getVertices()) {
      std::cerr << elem.x << ' ' << elem.y << std::endl;
    }*/
    return polygon->IsEqual(*this);
  }
  const auto* ellipse = dynamic_cast<const Ellipse*>(&another);
  if (ellipse != nullptr) {
    return ellipse->IsEqual(*this);
  }
  return false;
}

bool Polygon::isCongruentTo(const Shape& another) const {
  // std::cerr << "polygon::isCongruentTo\n";
  const auto* ellipse = dynamic_cast<const Ellipse*>(&another);
  if (ellipse != nullptr) {
    return false;
  }
  /*std::cerr << "polygon1\n";
  for (auto& elem : ver_) {
    std::cerr << elem.x << ' ' << elem.y << std::endl;
  }*/
  const auto* polygon = dynamic_cast<const Polygon*>(&another);
  /*std::cerr << "polygon2\n";
  for (auto& elem : polygon->ver_) {
    std::cerr << elem.x << ' ' << elem.y << std::endl;
  }*/
  std::vector<long double> side1, side2;
  for (size_t i = 1; i < ver_.size(); ++i) {
    side1.push_back(Vector(ver_[i - 1], ver_[i]).len());
  }
  side1.push_back(Vector(ver_.back(), ver_[0]).len());
  for (size_t i = 1; i < polygon->ver_.size(); ++i) {
    side2.push_back(Vector(polygon->ver_[i - 1], polygon->ver_[i]).len());
  }
  side2.push_back(Vector(polygon->ver_.back(), polygon->ver_[0]).len());

  /*std::cerr << "side1\n";
  for (auto& elem : side1) {
    std::cerr << elem << std::endl;
  }
  std::cerr << "side2\n";
  for (auto& elem : side2) {
    std::cerr << elem << std::endl;
  }
  std::cerr << "res: " << static_cast<bool>(::IsEqual(side1, side2))
            << std::endl;*/
  return ::IsEqual(side1, side2);
}

bool Ellipse::isCongruentTo(const Shape& another) const {
  const auto* polygon = dynamic_cast<const Polygon*>(&another);
  if (polygon != nullptr) {
    return false;
  }
  const auto* ellipse = dynamic_cast<const Ellipse*>(&another);
  return a_ == ellipse->a_ && b_ == ellipse->b_;
}

bool Circle::isCongruentTo(const Shape& another) const {
  return this->Ellipse::isCongruentTo(another);
}

bool Rectangle::isCongruentTo(const Shape& another) const {
  return this->Polygon::isCongruentTo(another);
}

bool Square::isCongruentTo(const Shape& another) const {
  return this->Polygon::isCongruentTo(another);
}

bool Triangle::isCongruentTo(const Shape& another) const {
  return this->Polygon::isCongruentTo(another);
}

bool Polygon::isSimilarTo(const Shape& another) const {
  const auto* ellipse = dynamic_cast<const Ellipse*>(&another);
  if (ellipse != nullptr) {
    return false;
  }
  // std::cerr << "polygon1\n";
  /*for(auto& elem : ver_) {
    std::cerr << elem.x << ' ' << elem.y << std::endl;
  }*/
  const auto* polygon = dynamic_cast<const Polygon*>(&another);
  // std::cerr << "polygon2\n";
  /*for(auto& elem : polygon->ver_) {
    std::cerr << elem.x << ' ' << elem.y << std::endl;
  }*/
  std::vector<long double> side1, side2;
  for (size_t i = 1; i < ver_.size(); ++i) {
    side1.push_back(Vector(ver_[i - 1], ver_[i]).len());
  }
  side1.push_back(Vector(ver_.back(), ver_[0]).len());
  for (size_t i = 1; i < polygon->ver_.size(); ++i) {
    side2.push_back(Vector(polygon->ver_[i - 1], polygon->ver_[i]).len());
  }
  side2.push_back(Vector(polygon->ver_.back(), polygon->ver_[0]).len());
  long double perimetr1 = 0;
  long double perimetr2 = 0;

  for (auto elem : side1) {
    perimetr1 += elem;
  }
  for (auto elem : side2) {
    perimetr2 += elem;
  }
  return ::IsEqual(side1, side2, perimetr2 / perimetr1);
}

bool Ellipse::isSimilarTo(const Shape& another) const {
  const auto* polygon = dynamic_cast<const Polygon*>(&another);
  if (polygon != nullptr) {
    return false;
  }
  const auto* ellipse = dynamic_cast<const Ellipse*>(&another);
  return a_ * ellipse->b_ == b_ * ellipse->a_;
}

bool Circle::isSimilarTo(const Shape& another) const {
  return this->Ellipse::isSimilarTo(another);
}

bool Rectangle::isSimilarTo(const Shape& another) const {
  return this->Polygon::isSimilarTo(another);
}

bool Square::isSimilarTo(const Shape& another) const {
  return this->Polygon::isSimilarTo(another);
}

bool Triangle::isSimilarTo(const Shape& another) const {
  return this->Polygon::isSimilarTo(another);
}

bool Polygon::containsPoint(const Point& point) const {
  bool flag_positive = false;
  bool flag_negative = false;
  for (size_t i = 0; i < ver_.size(); ++i) {
    (Vector(point, ver_[i]) * Vector(point, ver_[(i + 1) % ver_.size()]) >= 0
     ? flag_positive : flag_negative) = true;
  }
  return flag_negative ^ flag_positive;
}

bool Ellipse::containsPoint(const Point& point) const {
  return std::pow(point.x - center_.x, 2) / (a_ * a_) +
         std::pow(point.y - center_.y, 2) / (b_ * b_) <= 1;
}

bool Circle::containsPoint(const Point& point) const {
  return this->Ellipse::containsPoint(point);
}

bool Rectangle::containsPoint(const Point& point) const {
  return this->Polygon::containsPoint(point);
}

bool Square::containsPoint(const Point& point) const {
  return this->Polygon::containsPoint(point);
}

bool Triangle::containsPoint(const Point& point) const {
  return this->Polygon::containsPoint(point);
}

void Polygon::rotate(const Point& center, long double angle) {
  //std::cerr << "polygon::rotate\n";
  /*std::cerr << "angle: " << angle << std::endl;
  std::cerr << "center: " << center.x << ' ' << center.y << std::endl;
  std::cerr << "area: " << std::fixed << std::setprecision(5) << area()
            << std::endl;*/
  for (auto& point : ver_) {
    // std::cerr << point.x << ' ' << point.y << ' ';
    point.rotate(center, angle);
    // std::cerr << point.x << ' ' << point.y << std::endl;
  }
  // std::cerr << "new area: " << std::fixed << std::setprecision(5) << area()
  // << std::endl;
}

void Ellipse::rotate(const Point& center, long double angle) {
  std::cerr << "Ellipse::rotate\n";
  f1_.rotate(center, angle);
  f2_.rotate(center, angle);
  center_.rotate(center, angle);
}

void Circle::rotate(const Point& center, long double angle) {
  std::cerr << "Circle::rotate\n";
  Ellipse::rotate(center, angle);
}

void Rectangle::rotate(const Point& center, long double angle) {
  std::cerr << "Rectangle::rotate\n";
  Polygon::rotate(center, angle);
}

void Square::rotate(const Point& center, long double angle) {
  std::cerr << "Square::rotate\n";
  Polygon::rotate(center, angle);
}

void Triangle::rotate(const Point& center, long double angle) {
  std::cerr << "Triangle::rotate\n";
  Polygon::rotate(center, angle);
}

void Polygon::reflect(const Point& center) {
  std::cerr << "polygon::reflect\n";
  for (auto& point : ver_) {
    point.reflect(center);
  }
}

void Ellipse::reflect(const Point& center) {
  std::cerr << "Ellipse::reflect\n";
  f1_.reflect(center);
  f2_.reflect(center);
  center_.reflect(center);
}

void Circle::reflect(const Point& center) {
  std::cerr << "Circle::reflect\n";

  Ellipse::reflect(center);
}

void Rectangle::reflect(const Point& center) {
  std::cerr << "Rectangle::reflect\n";

  Polygon::reflect(center);
}

void Square::reflect(const Point& center) {
  std::cerr << "Square::reflect\n";

  Polygon::reflect(center);
}

void Triangle::reflect(const Point& center) {
  std::cerr << "Triangle::reflect\n";

  Polygon::reflect(center);
}

void Polygon::reflect(const Line& axis) {
  std::cerr << "Polygon::reflect\n";
  for (auto& point : ver_) {
    point.reflect(axis);
  }
}

void Ellipse::reflect(const Line& axis) {
  f1_.reflect(axis);
  f2_.reflect(axis);
  center_.reflect(axis);
}

void Circle::reflect(const Line& axis) {
  Ellipse::reflect(axis);
}

void Rectangle::reflect(const Line& axis) {
  Polygon::reflect(axis);
}

void Square::reflect(const Line& axis) {
  Polygon::reflect(axis);
}

void Triangle::reflect(const Line& axis) {
  Polygon::reflect(axis);
}

void Polygon::scale(const Point& center, long double coefficient) {
  for (auto& elem : ver_) {
    elem.scale(center, coefficient);
  }
}

void Ellipse::scale(const Point& center, long double coefficient) {
  f1_.scale(center, coefficient);
  f2_.scale(center, coefficient);
  *this = {f1_, f2_, 2 * a_ * coefficient};
}

void Circle::scale(const Point& center, long double coefficient) {
  Ellipse::scale(center, coefficient);
}

void Rectangle::scale(const Point& center, long double coefficient) {
  Polygon::scale(center, coefficient);
}

void Square::scale(const Point& center, long double coefficient) {
  Polygon::scale(center, coefficient);
}

void Triangle::scale(const Point& center, long double coefficient) {
  Polygon::scale(center, coefficient);
}