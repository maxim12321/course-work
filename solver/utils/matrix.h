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

protected:
  Vector GetVector(int i, int j, int step, int size);

  int rows_;
  int columns_;

private:
  std::vector<double> matrix_;
};

// outputs the matrix in transposed form
std::ostream &operator<<(std::ostream &out, const Matrix &matrix);

class Vector {
  friend class Matrix;

public:
  class Iterator {
    friend class Vector;

  public:
    Iterator() = default;
    Iterator(const Iterator &it);

    Iterator &operator++();

    friend Iterator operator+(Iterator iter, int offset) {
      Vector::Iterator new_iter(iter);
      new_iter.vector_index_ += offset;
      return new_iter;
    }

    friend Iterator operator-(Iterator iter, int offset) {
      Vector::Iterator new_iter(iter);
      new_iter.vector_index_ -= offset;
      return new_iter;
    }

    bool operator==(Iterator other) const;
    bool operator!=(Iterator other) const;

    double &operator*();
    const double &operator*() const;

    std::pair<int, int> Index();
  private:
    Vector *base_;
    int step_;
    int vector_index_;
    int matrix_start_index_;
    int matrix_column_size_;
  };

private:
  Vector(Matrix &base, int i, int j, int step, int size);

public:
  int GetSize() const;

  double &operator[](int index);
  const double &operator[](int index) const;

  Iterator begin();
  Iterator end();

private:
  Matrix &base_;
  int i_;
  int j_;
  int step_;
  int size_;
};
