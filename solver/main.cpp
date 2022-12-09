#include <cassert>
#include <iostream>

#include "solver.h"
#include <mpi.h>

#define ROOT 0

int main(int argc, char **argv) {
  int prs = 0;
  int p = 0;

  if (p == ROOT && argc != 3) {
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

  Solver solver(p, prs, &properties, [&](const Matrix &matrix) {
    out << matrix;
    out << "\n\n";
  });
  solver.Start();

  out.close();

  return 0;
}
