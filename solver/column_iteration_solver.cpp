#include "column_iteration_solver.h"

ColumnIterationSolver::ColumnIterationSolver(PropertiesManager* properties)
    : properties_(properties) {}

Matrix ColumnIterationSolver::CalculateNextIteration(const Matrix& previous,
                                                     const Matrix& semi_previous) {
    // TODO: calculate T^{s+1} from T^s and T^{s+1/2}
    return previous + semi_previous;
}
