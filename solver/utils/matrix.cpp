#include <iostream>
#include "matrix.h"

Matrix::Matrix() : Matrix(0, 0) {}

Matrix::Matrix(int n, int m)
  : rows_(n), columns_(m), row_values_(n, Vector(m)) {}

Matrix::Matrix(const std::initializer_list<Vector>& list) {
  rows_ = list.size();
  if (rows_ == 0) {
    columns_ = 0;
  }
  else {
    columns_ = list.begin()->GetSize();
  }
  for (const auto& row : list) {
    row_values_.push_back(row);
  }
}

Matrix Matrix::Diagonal(int size, int diagonal_value) {
  Matrix result(size, size);
  for (int i = 0; i < size; i++) {
    result[i][i] = diagonal_value;
  }
  return result;
}

Matrix Matrix::Diagonal(const Vector& vector) {
  int n = vector.GetSize();
  Matrix result(n, n);
  for (int i = 0; i < n; i++) {
    result[i][i] = vector[i];
  }
  return result;
}

int Matrix::GetRowCount() const {
  return rows_;
}

int Matrix::GetColumnCount() const {
  return columns_;
}

void Matrix::SwapRows(int row1, int row2) {
  std::swap(row_values_[row1], row_values_[row2]);
}

Matrix Matrix::Transposed() const {
  Matrix result(columns_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < columns_; j++) {
      result[j][i] = row_values_[i][j];
    }
  }
  return result;
}

void Matrix::Transpose() {
  *this = Transposed();
}

Vector Matrix::MultiplyAsColumn(const Vector& other) const {
  Vector result(rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.GetSize(); j++) {
      result[i] += row_values_[i][j] * other[j];
    }
  }
  return result;
}

Vector Matrix::MultiplyTransposedAsColumn(const Vector& other) const {
  Vector result(columns_);
  for (int i = 0; i < columns_; i++) {
    for (int j = 0; j < other.GetSize(); j++) {
      result[i] += row_values_[j][i] * other[j];
    }
  }
  return result;
}

Vector& Matrix::operator[](int index) {
  return row_values_[index];
}

const Vector& Matrix::operator[](int index) const {
  return row_values_[index];
}

Matrix Matrix::operator*(const Matrix& other) const {
  Matrix result(rows_, other.columns_);
  for (int k = 0; k < columns_; k++) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < other.columns_; j++) {
        result[i][j] += row_values_[i][k] * other[k][j];
      }
    }
  }
  return result;
}

Matrix Matrix::operator+(const Matrix& other) const {
  Matrix result = *this;
  for (int i = 0; i < rows_; i++) {
    result[i] += other[i];
  }
  return result;
}

Matrix& Matrix::operator+=(const Matrix& other) {
  *this = *this + other;
  return *this;
}

void Matrix::AddRow(const Vector& row) {
  if (columns_ == 0) {
    columns_ = row.GetSize();
  }

  row_values_.push_back(row);
  rows_++;
}

void Matrix::Store(std::ostream& out, int row_begin, int row_end) const {
  row_end = (row_end >= 0 ? std::min(row_end, rows_) : rows_);

  out << row_end - row_begin << " " << columns_ << "\n";
  for (int i = row_begin; i < row_end; ++i) {
    for (int j = 0; j < columns_; ++j) {
      out << row_values_[i][j] << " ";
    }
    out << "\n";
  }
}

Matrix Matrix::Load(std::istream& in) {
  std::string temp;
  int rows, columns;

  in >> temp;
  rows = std::stoi(temp);

  in >> temp;
  columns = std::stoi(temp);

  Matrix result(rows, columns);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      in >> temp;
      result[i][j] = std::stold(temp);
    }
  }
  return result;
}

std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
  out << matrix.GetRowCount() << ' ' << matrix.GetColumnCount() << '\n';
  out << "[";
  for (int i = 0; i < matrix.GetRowCount() - 1; i++) {
    out << matrix[i] << "\n ";
  }
  if (matrix.GetRowCount() > 0) {
    out << matrix[matrix.GetRowCount() - 1];
  }
  out << "]";
  return out;
}
