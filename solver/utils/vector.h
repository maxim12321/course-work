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

  long double& operator[](int index);
  const long double& operator[](int index) const;

 private:
  std::vector<long double> values_;
};

std::ostream& operator<<(std::ostream& out, const Vector& vector);
