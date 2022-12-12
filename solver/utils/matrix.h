#pragma once

#include <random>
#include <vector>

#include "vector.h"

class Matrix {
 public:
  Matrix();
  Matrix(int n, int m);
  Matrix(const std::initializer_list<Vector>& list);

  static Matrix Diagonal(int size, int diagonal_value = 1);
  static Matrix Diagonal(const Vector& vector);

  int GetRowCount() const;
  int GetColumnCount() const;

  void SwapRows(int row1, int row2);

  Matrix Transposed() const;
  void Transpose();

  Vector MultiplyAsColumn(const Vector& other) const;
  Vector MultiplyTransposedAsColumn(const Vector& other) const;

  void AddRow(const Vector& row);

  Vector& operator[](int index);
  const Vector& operator[](int index) const;

  Matrix operator*(const Matrix& other) const;

  Matrix operator+(const Matrix& other) const;
  Matrix& operator+=(const Matrix& other);

  void Store(std::ostream& out, int row_begin = 0, int row_end = -1) const;
  static Matrix Load(std::istream& in);

 private:
  int rows_;
  int columns_;

  std::vector<Vector> row_values_;
};

std::ostream& operator<<(std::ostream& out, const Matrix& matrix);
