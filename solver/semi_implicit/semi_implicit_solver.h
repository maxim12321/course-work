#pragma once

#include <functional>

#include "../properties_manager.h"
#include "../solver_base.h"
#include "../utils/matrix.h"
#include "column_iteration_solver.h"
#include "row_iteration_solver.h"

class SemiImplicitSolver : public SolverBase {
 public:
  SemiImplicitSolver(int p_rank, int p_size, PropertiesManager* properties, Callback callback);

  void Solve() override;

 private:
  // Calculate new values for *temperature_* based on current values,
  // invoke *on_layer_ready_* callback
  void CalculateNextLayer();

  void SchedulerRoutine();
  void WorkerRoutine();

  bool HasConverged(const Matrix& current, const Matrix& next);

 private:
  RowIterationSolver row_solver_;
  ColumnIterationSolver column_solver_;
};
