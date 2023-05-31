#pragma once

#include <mpi.h>
#include "../solver_base.h"
#include "../properties_wrapper.h"
#include "node_edge_info.h"
#include "../utils/result_saver.h"

class ExplicitSolver : public SolverBase, public PropertiesWrapper {
 public:
  ExplicitSolver(int p_rank,
                 int p_size,
                 PropertiesManager* properties,
                 std::string result_file_name,
                 int num_threads = 1);

  void Solve() override;

 private:
  const int kLoggingSkip = 10;

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

  std::string result_file_name_;

  int num_threads_;

  std::vector<MPI_Request> read_requests_;
  std::vector<MPI_Request> write_requests_;
};
