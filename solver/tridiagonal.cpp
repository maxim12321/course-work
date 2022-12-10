#include "tridiagonal.h"

#include <iostream>

const long double kEpsilon = 0.0000000001;

TridiagonalMatrix::TridiagonalMatrix(int N)
    : N_(N), matrix_(N, std::vector<double>(3, 0)) {}

int TridiagonalMatrix::GetN() const {
  return N_;
}

void TridiagonalMatrix::SwapRows(int row1, int row2) {
  std::swap(matrix_[row1], matrix_[row2]);
}

void TridiagonalMatrix::AddUpRows(int source, int dest, double source_mult) {
  for (int i = 0; i < 3; ++i) {
    matrix_[dest][i] += matrix_[source][i] * source_mult;
  }
}

void TridiagonalMatrix::LeftShiftRow(int row) {
  matrix_[row][0] = matrix_[row][1];
  matrix_[row][1] = matrix_[row][2];
  matrix_[row][2] = 0;
}

std::vector<double> &TridiagonalMatrix::operator[](int index) {
  return matrix_[index];
}
const std::vector<double> &TridiagonalMatrix::operator[](int index) const {
  return matrix_[index];
}

bool TridiagonalAlgorithm(TridiagonalMatrix &matrix, Vector &right) {
  bool result = Gauss(matrix, right) && ReverseGauss(matrix, right);
  return result;
}

bool Gauss(TridiagonalMatrix &matrix, Vector &right) {
  int n = matrix.GetN();

  for (int i = 0; i < n - 1; i++) {
    if (std::abs(matrix[i][0]) < std::abs(matrix[i + 1][0])) {
      matrix.SwapRows(i, i + 1);
      std::swap(right[i], right[i + 1]);
    }
    if (std::abs(matrix[i][0]) < kEpsilon) {
      return false;
    }

    long double coefficient = matrix[i + 1][0] / matrix[i][0];
    matrix.AddUpRows(i, i + 1, -coefficient);
    right[i + 1] -= right[i] * coefficient;

    matrix.LeftShiftRow(i + 1);
  }
  return true;
}

bool ReverseGauss(const TridiagonalMatrix &matrix, Vector &right) {
  int n = matrix.GetN();

  for (int i = n - 1; i >= 0; i--) {
    if (std::abs(matrix[i][0]) < kEpsilon) {
      return false;
    }
    right[i] /= matrix[i][0];

    if (i > 0) {
      right[i - 1] -= matrix[i - 1][1] * right[i];
    }
    if (i > 1) {
      right[i - 2] -= matrix[i - 2][2] * right[i];
    }
  }
  return true;
}
