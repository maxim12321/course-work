#pragma once

#include "../solver_base.h"
#include "../properties_wrapper.h"

class ExplicitSolver : public SolverBase {
 public:
  ExplicitSolver(int p_rank, int p_size, PropertiesManager* properties, Callback callback);

  void Solve() override;

 private:
  PropertiesWrapper props_;
};
