#pragma once

#include "properties_manager.h"
#include "utils/matrix.h"

class ColumnIterationSolver {
public:
    explicit ColumnIterationSolver(PropertiesManager* properties);

    // prev_iter = T^{s}, semi_previ_iter = T^{s+1/2} (from RowSolver)
    Matrix CalculateNextIteration(const Matrix& prev_iter,
                                  const Matrix& semi_prev_iter);

private:
    void LeftColumn(const Matrix& prev_iter,
                    const Matrix& semi_prev_iter,
                    Vector& column);
    void MiddleColumn(const Matrix& prev_iter,
                      const Matrix& semi_prev_iter,
                      Vector& column, int i);
    void RightColumn(const Matrix& prev_iter,
                     const Matrix& semi_prev_iter,
                     Vector& column);

    PropertiesManager* properties_;
    int N_;
    int M_;
};

