#include <iostream>
#include <memory>
#include <mpi.h>

#include "explicit/explicit_solver.h"
#include "semi_implicit/semi_implicit_solver.h"
#include "properties/properties_manager.h"

const int kMasterRank = 0;
int main(int argc, char** argv) {
  // MPI_Init(&argc, &argv);

  // int process_count;
  // int current_rank;

  // MPI_Comm_size(MPI_COMM_WORLD, &process_count);
  // MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);

  // if (current_rank == kMasterRank && argc < 3) {
  //   std::cerr << "The wrong number of arguments. Required 2, but got "
  //             << argc - 1 << "." << std::endl;
  //   std::cerr << "Usage: [input-file-path] [output-file-path] [solver-type, default: "
  //             << kSemiImplicitSolverType << "]"
  //             << std::endl;
  //   return 1;
  // }
  if (argc < 3) {
    std::cerr << "The wrong number of arguments. Required 2, but got "
              << argc - 1 << "." << std::endl;
    std::cerr << "Usage: [input-file-path] [output-file-path] [solver-type, default: "
              << kSemiImplicitSolverName << "]"
              << std::endl;
    return 1;
  }

  std::string input_filepath = argv[1];
  std::string output_filepath = argv[2];
  std::string solver_name = (argc == 3 ? kSemiImplicitSolverName : argv[3]);

  if (kSolverNameToType.count(solver_name) == 0) {
    std::string valid_solver_types;
    for (const auto &solvers : kSolverNameToType) {
        valid_solver_types += solvers.first + " ";
    }
    std::cerr << "Invalid solver type! Must be one of [ " << valid_solver_types <<  "]" << std::endl;
    return 1;
  }

  auto properties = PropertiesLoader::Load(input_filepath);

  std::ofstream out(output_filepath, std::ios::out);
  if (out.fail()) {
    std::cerr << "An error occurred while opening the input file" << std::endl;
    return 2;
  }

  auto callback = [&out](const Matrix& matrix) {
    out << matrix;
    out << "\n\n";
  };

  SolverType solver_type = kSolverNameToType.at(solver_name);
  std::shared_ptr<SolverBase> solver;

  switch (solver_type) {
  case SolverType::SemiImplicit:
    solver = std::make_shared<SemiImplicitSolver>(
        0,
        0,
        &properties,
        callback
    );
    break;
  case SolverType::Explicit:
    // solver = std::make_shared<ExplicitSolver>(
    //   0,
    //   0,
    //   &properties,
    //   callback
    // );
    break;
  }

  solver->Solve();

  out.close();

  // MPI_Finalize();
  return 0;
}
