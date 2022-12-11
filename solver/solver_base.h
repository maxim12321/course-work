#pragma once

#include <functional>

#include "properties_manager.h"
#include "utils/matrix.h"

class SolverBase {
 public:
  using Callback = std::function<void(const Matrix&)>;

  SolverBase(int p_rank, int p_size, PropertiesManager* properties, Callback callback)
      : p_rank_(p_rank),
        p_size_(p_size),
        nx_(properties->GetGridWidth()),
        nz_(properties->GetGridHeight()),
        properties_(properties),
        on_layer_ready_(std::move(callback)) {
    Initialize();
  }

  virtual void Solve() = 0;

 private:
  // Initialize fields with values from PropertiesManager...
  void Initialize() {
    current_temp_ = properties_->InitializeGrids(nx_, nz_);
    previous_temp_ = current_temp_;
    on_layer_ready_(current_temp_.Transposed());
  }

 protected:
  const int kSchedulerRank = 0;

  int p_rank_;
  int p_size_;

  int nx_;
  int nz_;

  Matrix current_temp_;
  Matrix previous_temp_;

  PropertiesManager* properties_;
  Callback on_layer_ready_;
};
