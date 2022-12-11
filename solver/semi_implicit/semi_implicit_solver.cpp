#include "semi_implicit_solver.h"

#include <iostream>

#include <chrono>
#include <utility>

SemiImplicitSolver::SemiImplicitSolver(int p_rank, int p_size, PropertiesManager* properties, Callback callback)
    : SolverBase(p_rank, p_size, properties, std::move(callback)),
      nx_(properties->GetGridWidth()),
      nz_(properties->GetGridHeight()),
      row_solver_(properties), column_solver_(properties) {
  Initialize();
}

void SemiImplicitSolver::Solve() {
  auto start = std::chrono::steady_clock::now();

  for (size_t i = 0; i < properties_->GetTimeLayers(); ++i) {
    CalculateNextLayer();
    if (i % 100 == 0) {
      // Print some debug info for long calculations
    }
  }

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
}

void SemiImplicitSolver::Initialize() {
  current_temp_ = properties_->InitializeGrids(nx_, nz_);
  previous_temp_ = current_temp_;
  on_layer_ready_(current_temp_.Transposed());
}

void SemiImplicitSolver::CalculateNextLayer() {
  for (int iteration = 0; iteration < properties_->GetMaxIterations();
       ++iteration) {
    Matrix semi_next =
        row_solver_.CalculateNextIteration(current_temp_, previous_temp_);
    Matrix next_temp =
        column_solver_.CalculateNextIteration(current_temp_, semi_next);

    if (HasConverged(current_temp_, next_temp)) {
      current_temp_ = next_temp;
      break;
    }
    current_temp_ = next_temp;
  }

  previous_temp_ = current_temp_;
  on_layer_ready_(current_temp_.Transposed());
}

bool SemiImplicitSolver::HasConverged(const Matrix& current, const Matrix& next) {
  for (int i = 0; i < current_temp_.GetRowCount(); ++i) {
    for (int j = 0; j < current_temp_.GetColumnCount(); ++j) {
      double delta = std::fabs(current[i][j] - next[i][j]);

      if (delta > properties_->GetEpsilon1() * current[i][j] +
          properties_->GetEpsilon2()) {
        return false;
      }
    }
  }
  return true;
}

void SemiImplicitSolver::SchedulerRoutine() {
}

void SemiImplicitSolver::WorkerRoutine() {
}
