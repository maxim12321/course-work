#pragma once

#include "properties_manager.h"
#include "utils/matrix.h"

class RowIterationSolver {
public:
    explicit RowIterationSolver(PropertiesManager* properties);

    Matrix CalculateNextLayer(const Matrix& previous);

private:
    void BottomRow(const Matrix& previous, Vector& row);
    void MiddleRow(const Matrix& previous, Vector& row);
    void TopRow(const Matrix& previous, Vector& row);

private:
    PropertiesManager* properties_;
};

