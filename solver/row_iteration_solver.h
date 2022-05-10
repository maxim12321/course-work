#pragma once

#include "properties_manager.h"
#include "utils/matrix.h"

class RowIterationSolver {
public:
    explicit RowIterationSolver(PropertiesManager* properties);

    Matrix CalculateNextLayer(const Matrix& previous);

private:
    PropertiesManager* properties_;
};

