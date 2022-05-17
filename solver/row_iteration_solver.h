#pragma once

#include "properties_manager.h"
#include "utils/matrix.h"

class RowIterationSolver {
public:
    explicit RowIterationSolver(PropertiesManager* properties);

    Matrix CalculateNextIteration(const Matrix& prev_iter, const Matrix& prev_layer);

private:
    void BottomRow(const Matrix& prev_iter, const Matrix& prev_layer, Vector& row);
    void MiddleRow(const Matrix& prev_iter, const Matrix& prev_layer, Vector& row, int k);
    void TopRow(const Matrix& prev_iter, const Matrix& prev_layer, Vector& row);

private:
    PropertiesManager* properties_;

    int N_;
    int M_;
};

