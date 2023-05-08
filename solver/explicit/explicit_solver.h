#pragma once

#include "../solver_base.h"
#include "../properties_wrapper.h"
#include "node_edge_info.h"

class ExplicitSolver : public SolverBase, public PropertiesWrapper {
 public:
  ExplicitSolver(int p_rank,
                 int p_size,
                 PropertiesManager* properties,
                 const std::string& result_file_name);

  void Solve() override;

 private:
  void PrepareNodeEdges();

  void CalculateNextLayer();

  long double GetNodeValue(int i, int k);

  long double GetNodeValue(NodeEdgeInfo* node_edge_info);

  long double T(int i, int k);

 private:
  int rows_per_process_;
  int row_begin_;
  int row_end_;

  std::vector<std::vector<NodeEdgeInfo>> nodes;

  std::ofstream output_;
};
