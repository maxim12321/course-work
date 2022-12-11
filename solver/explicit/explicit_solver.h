#pragma once

#include "../solver_base.h"
#include "../properties_wrapper.h"

class ExplicitSolver : public SolverBase, public PropertiesWrapper {
 public:
  ExplicitSolver(int p_rank, int p_size, PropertiesManager* properties, Callback callback);

  void Solve() override;

 private:
  void CalculateNextLayer();

  long double GetNodeValue(int i, int k);

  long double T(int i, int k);
};
