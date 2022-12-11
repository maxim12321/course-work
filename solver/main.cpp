#include <iostream>
#include <memory>
#include <mpi.h>

#include "semi_implicit/semi_implicit_solver.h"

const int kMasterRank = 0;
const std::string kSemiImplicitSolverType = "semi-implicit";
const std::string kExplicitSolverType = "explicit";

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int process_count;
  int current_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &process_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);

  if (current_rank == kMasterRank && argc < 3) {
    std::cerr << "The wrong number of arguments. Required 2, but got "
              << argc - 1 << "." << std::endl;
    std::cerr << "Usage: [input-file-path] [output-file-path] [solver-type, default: "
              << kSemiImplicitSolverType << "]"
              << std::endl;
    return 1;
  }

  std::string input_filepath = argv[1];
  std::string output_filepath = argv[2];
  std::string solver_type = (argc == 3 ? kSemiImplicitSolverType : argv[3]);

  if (solver_type != kSemiImplicitSolverType && solver_type != kExplicitSolverType) {
    std::cerr << "Invalid solver type! Must be one of [" << kSemiImplicitSolverType << ", "
              << kExplicitSolverType << "]" << std::endl;
    return 1;
  }

  auto properties = PropertiesManagerFabric::Create(input_filepath);

  std::ofstream out(output_filepath, std::ios::out);
  if (out.fail()) {
    std::cerr << "An error occurred while opening the input file" << std::endl;
    return 2;
  }

  auto callback = [&out](const Matrix& matrix) {
    out << matrix;
    out << "\n\n";
  };

  std::shared_ptr<SolverBase> solver;

  if (solver_type == kSemiImplicitSolverType) {
    solver = std::make_shared<SemiImplicitSolver>(
        current_rank,
        process_count,
        &properties,
        callback
    );
  } else if (solver_type == kExplicitSolverType) {
    // TODO: solver = std::make_shared<ExplicitSolver>(...);
  }

  solver->Solve();

  out.close();

  MPI_Finalize();
  return 0;
}
