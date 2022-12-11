#include "explicit_solver.h"

ExplicitSolver::ExplicitSolver(int p_rank,
                               int p_size,
                               PropertiesManager* properties,
                               SolverBase::Callback callback)
    : SolverBase(p_rank, p_size, properties, std::move(callback)),
      props_(properties) {
}

void ExplicitSolver::Solve() {
}
