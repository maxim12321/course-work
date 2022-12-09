#pragma once

#include <functional>

#include "column_iteration_solver.h"
#include "properties_manager.h"
#include "row_iteration_solver.h"
#include "utils/matrix.h"

class Solver {
public:
  using Callback = std::function<void(const Matrix &)>;

  Solver(int p_rank, int p_size, PropertiesManager *properties, Callback callback);

  void Start();

private:
  // Initialize fields with values from PropertiesManager...
  void Initialize();

  // Calculate new values for *temperature_* based on current values,
  // invoke *on_layer_ready_* callback
  void CalculateNextLayer();

  void SchedulerRoutine();
  void WorkerRoutine();

  bool HasConverged(const Matrix &current, const Matrix &next);

private:
  const int kSchedulerRank = 0;

  int p_rank_;
  int p_size_;

  PropertiesManager *properties_;
  Callback on_layer_ready_;

  int nx_;
  int nz_;

  Matrix current_temp_;
  Matrix previous_temp_;

  RowIterationSolver row_solver_;
  ColumnIterationSolver column_solver_;
};
