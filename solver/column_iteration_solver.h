#pragma once

#include "properties_manager.h"
#include "utils/matrix.h"

class ColumnIterationSolver {
public:
    explicit ColumnIterationSolver(PropertiesManager* properties);

    // previous = T^{s}, semi_previous = T^{s+1/2} (from RowSolver)
    Matrix CalculateNextLayer(const Matrix& previous,
                              const Matrix& semi_previous);

private:
    PropertiesManager* properties_;
};

