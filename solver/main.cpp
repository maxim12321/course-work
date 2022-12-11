#include <cassert>
#include <iostream>

#include "semi_implicit/semi_implicit_solver.h"
#include <mpi.h>

const int kMasterRank = 0;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int process_count;
  int current_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &process_count);
  MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);

  if (current_rank == kMasterRank && argc != 3) {
    std::cerr << "The wrong number of arguments. Required 2, but got "
              << argc - 1 << "." << std::endl;
    std::cerr << "Usage: [input-file-path] [output-file-path]" << std::endl;
    return 1;
  }

  std::string input_filepath = argv[1];
  std::string output_filepath = argv[2];
  auto properties = PropertiesManagerFabric::Create(input_filepath);

  std::ofstream out(output_filepath, std::ios::out);
  if (out.fail()) {
    std::cerr << "An error occurred while opening the input file" << std::endl;
    assert(false);
  }

  SemiImplicitSolver solver(current_rank, process_count, &properties, [&](const Matrix &matrix) {
    out << matrix;
    out << "\n\n";
  });
  solver.Solve();

  out.close();

  MPI_Finalize();
  return 0;
}
