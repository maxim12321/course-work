#include "matrix.h"

#include <ostream>

Matrix::Matrix() : Matrix(0, 0) {}

Matrix::Matrix(int n, int m) : rows_(n), columns_(m), matrix_(n * m, 0) {}

int Matrix::GetRowCount() const { return rows_; }

int Matrix::GetColumnCount() const { return columns_; }

Vector Matrix::GetRow(int i) { return GetVector(i, 0, 1, columns_); }

Vector Matrix::GetColumn(int i) { return GetVector(0, i, columns_, rows_); }

double &Matrix::operator()(int i, int j) { return matrix_[i * columns_ + j]; }

const double &Matrix::operator()(int i, int j) const {
  return matrix_[i * columns_ + j];
}

Vector Matrix::GetVector(int i, int j, int step, int size) {
  return Vector(*this, i, j, step, size);
}

// outputs the matrix in transposed form
std::ostream &operator<<(std::ostream &out, const Matrix &matrix) {
  auto n = matrix.GetRowCount();
  auto m = matrix.GetColumnCount();

  auto print_column = [&](int j) {
    out << "[";
    for (int i = 0; i < n - 1; i++) {
      out << matrix(i, j) << " ";
    }
    if (n > 0) {
      out << matrix(n - 1, j);
    }
    out << "]";
  };

  out << m << ' ' << n << '\n';
  out << "[";
  for (int j = 0; j < m - 1; j++) {
    print_column(j);
    out << "\n ";
  }
  if (m > 0) {
    print_column(m - 1);
  }
  out << "]";
  return out;
}

Vector::Vector(Matrix &base, int i, int j, int step, int size)
    : base_(base), i_(i), j_(j), step_(step), size_(size) {}

int Vector::GetSize() const { return size_; }

double &Vector::operator[](int index) {
  return base_.matrix_[i_ * base_.columns_ + j_ + index * step_];
}

const double &Vector::operator[](int index) const {
  return base_.matrix_[i_ * base_.columns_ + j_ + index * step_];
}

Vector::Iterator Vector::begin() {
  Iterator it;
  it.base_ = this;
  it.step_ = step_;
  it.vector_index_ = 0;
  it.matrix_column_size_ = base_.GetColumnCount();
  it.matrix_start_index_ = i_ * it.matrix_column_size_ + j_;

  return it;
}

Vector::Iterator Vector::end() {
  Iterator it;
  it.base_ = this;
  it.step_ = step_;
  it.vector_index_ = size_;
  it.matrix_column_size_ = base_.GetColumnCount();
  it.matrix_start_index_ = i_ * it.matrix_column_size_ + j_;

  return it;
}

Vector::Iterator::Iterator(const Vector::Iterator &it)
  : base_(it.base_), step_(it.step_), vector_index_(it.vector_index_),
    matrix_start_index_(it.matrix_start_index_), matrix_column_size_(it.matrix_column_size_) {
}

Vector::Iterator &Vector::Iterator::operator++() {
  ++vector_index_;
  return *this;
}

bool Vector::Iterator::operator==(Vector::Iterator other) const {
  return vector_index_ == other.vector_index_;
}

bool Vector::Iterator::operator!=(Vector::Iterator other) const {
  return vector_index_ != other.vector_index_;
}

double &Vector::Iterator::operator*() {
  return base_->operator[](vector_index_);
}

const double &Vector::Iterator::operator*() const {
  return base_->operator[](vector_index_);
}

std::pair<int, int> Vector::Iterator::Index() {
  int index = matrix_start_index_ + vector_index_ * step_;
  return {index / matrix_column_size_, index % matrix_column_size_};
}
