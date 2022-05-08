#include "vector.h"

Vector::Vector() : Vector(0) {}

Vector::Vector(int size, int default_value) : values_(size, default_value) {}

Vector::Vector(const std::initializer_list<long double>& list)
    : values_(list) {}

int Vector::GetSize() const {
  return values_.size();
}

void Vector::Shift(int amount) {
  int n = GetSize();
  Vector result(n);
  for (int i = std::max(0, -amount); i < n && i + amount < n; i++) {
    result[i + amount] = values_[i];
  }
  *this = result;
}

void Vector::AddValue(long double value) {
  values_.push_back(value);
}

void Vector::PopValue() {
  values_.pop_back();
}

long double& Vector::operator[](int index) {
  return values_[index];
}

const long double& Vector::operator[](int index) const {
  return values_[index];
}

Vector Vector::operator+(const Vector& other) const {
  Vector result = *this;
  result += other;
  return result;
}

Vector& Vector::operator+=(const Vector& other) {
  Add(other);
  return *this;
}

Vector Vector::operator-(const Vector& other) const {
  return *this + (other * -1);
}

Vector Vector::operator*(long double value) const {
  Vector result = *this;
  result *= value;
  return result;
}

Vector& Vector::operator*=(long double value) {
  for (int i = 0; i < GetSize(); i++) {
    values_[i] *= value;
  }
  return *this;
}

/**
 *
 * Adds [left; right) segment of a shifted vector, multiplied by a given
 * coefficient
 * @param other Vector to add
 * @param multiplier Coefficient of *other* vector
 * @param shift Shift of *other* vector, shift < 0: left, shift > 0: right
 * @param left Begin of the *other* vector, inclusive
 * @param right End of the *other* vector, exclusive
 */
void Vector::Add(const Vector& other, long double multiplier, int shift,
                 int left, int right) {
  int begin = std::max(left, -shift);

  for (int index = begin;
       index < other.GetSize() && index + shift < GetSize() && index != right;
       index++) {
    values_[index + shift] += other[index] * multiplier;
  }
}

void Vector::Multiply(long double multiplier, int left, int right) {
  for (int index = left; index < values_.size() && index != right; index++) {
    values_[index] *= multiplier;
  }
}

long double Vector::Dot(const Vector& other) const {
  long double result = 0;
  for (int i = 0; i < GetSize(); i++) {
    result += values_[i] * other[i];
  }
  return result;
}

long double Vector::Length() const {
  return std::sqrt(Dot(*this));
}

std::ostream& operator<<(std::ostream& out, const Vector& vector) {
  out << "[";
  for (int i = 0; i < vector.GetSize() - 1; i++) {
    out << vector[i] << ", ";
  }
  if (vector.GetSize() > 0) {
    out << vector[vector.GetSize() - 1];
  }
  out << "]";
  return out;
}
