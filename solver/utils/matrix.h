#pragma once

#include <random>
#include <vector>

class Vector;

class Matrix {
  friend class Vector;

public:
  Matrix();
  Matrix(int n, int m);

  int GetRowCount() const;
  int GetColumnCount() const;

  Vector GetRow(int i);
  Vector GetColumn(int i);

  double &operator()(int i, int j);
  const double &operator()(int i, int j) const;

private:
  int rows_;
  int columns_;

  std::vector<double> matrix_;
};

// outputs the matrix in transposed form
std::ostream &operator<<(std::ostream &out, const Matrix &matrix);

class Vector {
  friend class Matrix;

private:
  Vector(Matrix &base, int i, int j, int step, int size);

public:
  int GetSize() const;

  double &operator[](int index);
  const double &operator[](int index) const;

private:
  Matrix &base_;
  int i_;
  int j_;
  int step_;
  int size_;
};
