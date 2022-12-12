#pragma once

#include <cmath>
#include <cstdint>
#include <ostream>
#include <vector>

class Vector {
 public:
  Vector();
  explicit Vector(int size, int default_value = 0);
  Vector(const std::initializer_list<long double>& list);

  int GetSize() const;

  void Shift(int amount);

  void AddValue(long double value = 0);
  void PopValue();

  long double& operator[](int index);
  const long double& operator[](int index) const;

  Vector operator+(const Vector& other) const;
  Vector& operator+=(const Vector& other);

  Vector operator-(const Vector& other) const;

  Vector operator*(long double value) const;
  Vector& operator*=(long double value);

  void Add(const Vector& other, long double multiplier = 1, int shift = 0,
           int left = 0, int right = -1);

  void Multiply(long double multiplier, int left = 0, int right = -1);

  long double Dot(const Vector& other) const;
  long double Length() const;

  long double* Data();
  const long double* Data() const;

 private:
  std::vector<long double> values_;
};

std::ostream& operator<<(std::ostream& out, const Vector& vector);
