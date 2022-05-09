#include "row_iteration_solver.h"

RowIterationSolver::RowIterationSolver(PropertiesManager* properties)
    : properties_(properties) {}

Matrix RowIterationSolver::CalculateNextLayer(const Matrix& previous) {
    // TODO: calculate T^{s+1/2} from T^s
    return previous;
}
