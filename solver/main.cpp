#include <iostream>
#include <memory>
#include <mpi.h>

#include "explicit/explicit_solver.h"
#include "semi_implicit/semi_implicit_solver.h"

const int kMasterRank = 0;
const std::string kSemiImplicitSolverType = "semi-implicit";
const std::string kExplicitSolverType = "explicit";
const std::string kProcessResultsFilePrefix = "results_";

std::string GetResultsFileName(int process_rank) {
  return kProcessResultsFilePrefix + std::to_string(process_rank) + ".txt";
}

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

  std::ofstream process_results(GetResultsFileName(current_rank), std::ios::out);
  auto callback = [&process_results](const Matrix& matrix) {
    matrix.Store(process_results);
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
    solver = std::make_shared<ExplicitSolver>(
        current_rank,
        process_count,
        &properties,
        GetResultsFileName(current_rank)
    );
  }

  solver->Solve();

  // Wait for all processes to finish
  MPI_Barrier(MPI_COMM_WORLD);

  if (current_rank == kMasterRank) {
    std::ofstream out(output_filepath, std::ios::out);
    if (out.fail()) {
      std::cerr << "An error occurred while opening the input file" << std::endl;
      return 2;
    }

    int time_layers = properties.GetTimeLayers();

    std::vector<std::ifstream> results;
    results.reserve(process_count);
    for (int i = 0; i < process_count; ++i) {
      results.emplace_back(GetResultsFileName(i), std::ios::in);

      int curr_layers;
      results.back() >> curr_layers;
      time_layers = std::min(time_layers, curr_layers);
    }

    out << time_layers << "\n";

    for (int time = 0; time < time_layers; ++time) {
      Matrix result(0, 0);

      for (int i = 0; i < process_count; ++i) {
        Matrix rows = Matrix::Load(results[i]);
        for (int row = 0; row < rows.GetRowCount(); ++row) {
          result.AddRow(rows[row]);
        }
      }

      out << result.Transposed() << "\n\n";
    }

    out.close();
  }

  MPI_Finalize();
  return 0;
}
