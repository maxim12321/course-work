#include "tridiagonal.h"

#include <iostream>

const long double kEpsilon = 0.0000000001;

bool TridiagonalAlgorithm(Matrix& matrix, Vector& right) {
//    std::cout << matrix << "\n" << right << "\n";
    bool result = Gauss(matrix, right) && ReverseGauss(matrix, right);
//    std::cout << right << "\n" << std::endl;
    return result;
}

bool Gauss(Matrix& matrix, Vector& right) {
  int n = matrix.GetRowCount();

  for (int i = 0; i < n - 1; i++) {
    if (std::abs(matrix[i][0]) < std::abs(matrix[i + 1][0])) {
      matrix.SwapRows(i, i + 1);
      std::swap(right[i], right[i + 1]);
    }
    if (std::abs(matrix[i][0]) < kEpsilon) {
      return false;
    }

    long double coefficient = matrix[i + 1][0] / matrix[i][0];
    matrix[i + 1].Add(matrix[i], -coefficient);
    right[i + 1] -= right[i] * coefficient;

    matrix[i + 1].Shift(-1);
  }
  return true;
}

bool ReverseGauss(const Matrix& matrix, Vector& right) {
  int n = matrix.GetRowCount();

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
