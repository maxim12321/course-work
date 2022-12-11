#pragma once

#include "../properties_wrapper.h"
#include "../utils/matrix.h"

class ColumnIterationSolver : public PropertiesWrapper {
public:
  explicit ColumnIterationSolver(PropertiesManager *properties);

  // previous = T^{s}, semi_previous = T^{s+1/2} (from RowSolver)
  Matrix CalculateNextIteration(const Matrix &prev_iter,
                                const Matrix &semi_prev_iter);

private:
  void BuildTridiagonal(const Matrix &prev_iter, const Matrix &semi_prev_iter,
                        Vector &column, int i);
  void LeftColumn(const Matrix &prev_iter, const Matrix &semi_prev_iter,
                  Vector &column);
  void MiddleColumn(const Matrix &prev_iter, const Matrix &semi_prev_iter,
                    Vector &column, int i);
  void RightColumn(const Matrix &prev_iter, const Matrix &semi_prev_iter,
                   Vector &column);

private:
  int N_;
  int M_;

  Matrix tridiagonal_;
  Matrix next_;
};
