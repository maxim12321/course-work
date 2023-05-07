#pragma once

#include <functional>

#include "tridiagonal.h"
#include "../properties/properties_manager.h"
#include "../solver_base.h"
#include "../utils/matrix.h"

using DerivateCoeffs = std::tuple<double, double, double>;

class SemiImplicitSolver : public SolverBase {
 public:
  SemiImplicitSolver(int p_rank, int p_size, Properties* properties, Callback callback);

  void Solve() override;

 private:
  // Calculate new values for *temperature_* based on current values,
  // invoke *on_layer_ready_* callback
  void CalculateNextLayer();

  void SchedulerRoutine();
  void WorkerRoutine();

  void InitDerivativeApproxs();

  bool HasConverged(const ComputationGrid &current, const ComputationGrid &next);

  void RunRowIters(const ComputationGrid &prev_layer, const ComputationGrid &prev_iter);
  void BuildRowIterSystem(const ComputationGrid &prev_layer, const ComputationGrid &prev_iter, Vector &row);

  void RunColumnIters(const ComputationGrid &prev_iter, const ComputationGrid &semi_next_iter);
  void BuildColumnIterSystem(const ComputationGrid &prev_iter, const ComputationGrid &semi_next_iter, Vector &column);

  DerivateCoeffs XForwardDerivateCoeffs(int x, int z, double from_temp = 0.l, double to_temp = 0.l);
  DerivateCoeffs XBackwardDerivateCoeffs(int x, int z, double from_temp = 0.l, double to_temp = 0.l);
  DerivateCoeffs ZForwardDerivateCoeffs(int x, int z, double from_temp = 0.l, double to_temp = 0.l);
  DerivateCoeffs ZBackwardDerivateCoeffs(int x, int z, double from_temp = 0.l, double to_temp = 0.l);

 private:
  std::vector<std::vector<std::pair<int, int>>> x_derivates_ids_grid_;
  std::vector<std::function<DerivateCoeffs(int, int, double, double)>> x_derivates_;

  std::vector<std::vector<std::pair<int, int>>> z_derivates_ids_grid_;
  std::vector<std::function<DerivateCoeffs(int, int, double, double)>> z_derivates_;

  ComputationGrid semi_next_temp_;
  ComputationGrid next_temp_;

  TridiagonalMatrix row_system_;
  TridiagonalMatrix column_system_;
};
