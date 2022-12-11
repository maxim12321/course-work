#pragma once

#include <vector>

#include "utils/matrix.h"

class TridiagonalMatrix {
public:
  TridiagonalMatrix(int N);

  int GetN() const;

  void SwapRows(int row1, int row2);

  void AddUpRows(int source, int dest, double source_mult);

  void LeftShiftRow(int row);

  std::vector<double> &operator[](int index);
  const std::vector<double> &operator[](int index) const;

private:
  int N_;
  std::vector<std::vector<double>> matrix_;
};

bool TridiagonalAlgorithm(TridiagonalMatrix &matrix, Vector &right);

bool Gauss(TridiagonalMatrix &matrix, Vector &right);

bool ReverseGauss(const TridiagonalMatrix &matrix, Vector &right);
