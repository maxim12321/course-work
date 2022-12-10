#include "vector.h"

Vector::Vector() : Vector(0) {}

Vector::Vector(int size, int default_value) : values_(size, default_value) {}

Vector::Vector(const std::initializer_list<long double>& list)
  : values_(list) {}

int Vector::GetSize() const {
  return values_.size();
}


long double& Vector::operator[](int index) {
  return values_[index];
}

const long double& Vector::operator[](int index) const {
  return values_[index];
}

std::ostream& operator<<(std::ostream& out, const Vector& vector) {
  out << "[";
  for (int i = 0; i < vector.GetSize() - 1; i++) {
    out << vector[i] << " ";
  }
  if (vector.GetSize() > 0) {
    out << vector[vector.GetSize() - 1];
  }
  out << "]";
  return out;
}
