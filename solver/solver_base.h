#pragma once

#include <functional>

#include "properties_manager.h"
#include "utils/matrix.h"

class SolverBase {
 public:
  using Callback = std::function<void(const Matrix&)>;

  SolverBase(int p_rank, int p_size, PropertiesManager* properties, Callback callback)
      : p_rank_(p_rank), p_size_(p_size),
        properties_(properties),
        on_layer_ready_(std::move(callback)) {}

  virtual void Solve() = 0;

 protected:
  const int kSchedulerRank = 0;

  int p_rank_;
  int p_size_;

  PropertiesManager* properties_;
  Callback on_layer_ready_;
};
