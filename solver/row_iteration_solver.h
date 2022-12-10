#pragma once

#include "properties_wrapper.h"
#include "tridiagonal.h"
#include "utils/matrix.h"

class RowIterationSolver : public PropertiesWrapper {
public:
  explicit RowIterationSolver(PropertiesManager *properties);

  Matrix CalculateNextIteration(const Matrix &prev_iter,
                                const Matrix &prev_layer);

private:
  void BuildTridiagonal(const Matrix &prev_iter, const Matrix &prev_layer, Vector &row, int k);

  void BottomRow(const Matrix &prev_iter, const Matrix &prev_layer,
                 Vector &row);
  void MiddleRow(const Matrix &prev_iter, const Matrix &prev_layer,
                 Vector &row, int k);
  void TopRow(const Matrix &prev_iter, const Matrix &prev_layer,
              Vector &row);

private:
  int N_;
  int M_;

  TridiagonalMatrix tridiagonal_;
  Matrix next_;
};
