#include "explicit_solver.h"

#include <chrono>
#include <iostream>
#include <mpi.h>

ExplicitSolver::ExplicitSolver(int p_rank,
                               int p_size,
                               PropertiesManager* properties,
                               const std::string& result_file_name)
    : SolverBase(p_rank, p_size, properties, [](auto) {}),
      PropertiesWrapper(properties),
      rows_per_process_((nx_ + 2 + p_size - 1) / p_size),
      row_begin_(p_rank * rows_per_process_),
      row_end_((p_rank + 1) * rows_per_process_),
      output_(result_file_name, std::ios::out) {
  std::cout << "[#" << p_rank << "]: " << row_begin_ << "-" << row_end_ << std::endl;
}

void ExplicitSolver::Solve() {
  output_ << (properties_->GetTimeLayers() / kLoggingSkip) + 1 << "\n";

  PrepareNodeEdges();

  auto start = std::chrono::steady_clock::now();

  current_temp_.Store(output_, row_begin_, row_end_);

  for (size_t i = 1; i <= properties_->GetTimeLayers(); ++i) {
    CalculateNextLayer();
    if (i % kLoggingSkip == 0) {
      current_temp_.Store(output_, row_begin_, row_end_);
    }
  }

  output_.close();

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
}

void ExplicitSolver::CalculateNextLayer() {
  UpdateTemperatureLambda(previous_temp_);

  if (p_rank_ > 0) {
    MPI_Status status;
    MPI_Recv(previous_temp_[row_begin_ - 1].Data(), nz_ + 2, MPI_LONG_DOUBLE, p_rank_ - 1,
             100, MPI_COMM_WORLD, &status);

    MPI_Send(previous_temp_[row_begin_].Data(), nz_ + 2, MPI_LONG_DOUBLE, p_rank_ - 1,
             101, MPI_COMM_WORLD);
  }
  if (p_rank_ + 1 < p_size_) {
    MPI_Send(previous_temp_[row_end_ - 1].Data(), nz_ + 2, MPI_LONG_DOUBLE, p_rank_ + 1,
             100, MPI_COMM_WORLD);

    MPI_Status status;
    MPI_Recv(previous_temp_[row_end_].Data(), nz_ + 2, MPI_LONG_DOUBLE, p_rank_ + 1,
             101, MPI_COMM_WORLD, &status);
  }

  int row_end = std::min(row_end_, nx_ + 2);

// #pragma omp parallel for num_threads(4)
  for (int i = row_begin_; i < row_end; ++i) {
    for (int k = 0; k <= nz_ + 1; ++k) {
      // current_temp_[i][k] = GetNodeValue(i, k);
      current_temp_[i][k] = GetNodeValue(&nodes[i][k]);
    }
  }

  previous_temp_ = current_temp_;
}

void ExplicitSolver::PrepareNodeEdges() {
  const int N = nx_;
  const int M = nz_;

  nodes.resize(N + 2, std::vector<NodeEdgeInfo>(M + 2));

  for (int i = 0; i <= N + 1; i++) {
    for (int k = 0; k <= M + 1; k++) {
      NodeEdgeInfo& node = nodes[i][k];
      node.i = i;
      node.k = k;

      if (i == 0) {
        node.width = 0.5 * dx(1);
      } else if (i == N + 1) {
        node.width = 0.5 * dx(N);
      } else {
        node.width = dx(i);
      }

      if (k == 0) {
        node.height = 0.5 * dz(1);
      } else if (k == M + 1) {
        node.height = 0.5 * dz(M);
      } else {
        node.height = dz(k);
      }

      if (i == 1 || i == N) {
        node.width *= 0.75;
      }
      if (k == 1 || k == M) {
        node.height *= 0.75;
      }

      if (i == 0) {
        node.left_edge = EdgeType::kAir;
        node.right_edge = (k == 0 || k == M + 1) ? EdgeType::kMixed : EdgeType::kMaterial;
        node.bottom_edge = (k == 0) ? EdgeType::kAir : EdgeType::kMixed;
        node.top_edge = (k == M + 1) ? EdgeType::kAir : EdgeType::kMixed;
      } else if (i == N + 1) {
        node.left_edge = (k == 0 || k == M + 1) ? EdgeType::kMixed : EdgeType::kMaterial;
        node.right_edge = EdgeType::kAir;
        node.bottom_edge = (k == 0) ? EdgeType::kAir : EdgeType::kMixed;
        node.top_edge = (k == M + 1) ? EdgeType::kAir : EdgeType::kMixed;
      } else if (k == 0) {
        node.left_edge = EdgeType::kMixed;
        node.right_edge = EdgeType::kMixed;
        node.top_edge = EdgeType::kMaterial;
        node.bottom_edge = EdgeType::kAir;
      } else if (k == M + 1) {
        node.left_edge = EdgeType::kMixed;
        node.right_edge = EdgeType::kMixed;
        node.top_edge = EdgeType::kAir;
        node.bottom_edge = EdgeType::kMaterial;
      } else {
        node.left_edge = EdgeType::kMaterial;
        node.right_edge = EdgeType::kMaterial;
        node.top_edge = EdgeType::kMaterial;
        node.bottom_edge = EdgeType::kMaterial;
      }
    }
  }
}

long double ExplicitSolver::GetNodeValue(NodeEdgeInfo* node) {
  int i = node->i;
  int k = node->k;

  long double right_side = manager_->GetHeatX(i, k) + manager_->GetHeatZ(i, k);

  auto get_air = [&](double alpha, long double size) -> long double {
    return -alpha * (T(i, k) - t_out) / size;
  };

  auto get_material = [&](int dx, int dk) -> long double {
    long double lbd = lambda(i + 0.5 * dx, k + 0.5 * dk);
    if (dx < 0 || dk < 0) {
      lbd *= -1;
    }

    long double size, next_size;
    if (dx == 0) {
      size = node->height;
      next_size = nodes[i][k + dk].height;
    } else {
      size = node->width;
      next_size = nodes[i + dx][k].width;
    }

    return lbd * (T(i + dx, k + dk) - T(i, k)) / size / (0.5 * size + 0.5 * next_size);
  };

  auto get_for_edge =
      [&](EdgeType type, double alpha, int dx, int dk, long double size) -> long double {
        switch (type) {
          case EdgeType::kAir:
            return get_air(alpha, size);
          case EdgeType::kMaterial:
            return get_material(dx, dk);
          case EdgeType::kMixed:
            return 0.5 * get_air(alpha, size) + 0.5 * get_material(dx, dk);
          case EdgeType::kNone:
            return 0;
        }
      };

  // left
  right_side -= get_for_edge(node->left_edge, alpha1, -1, 0, node->width);
  // right
  right_side += get_for_edge(node->right_edge, alpha2, 1, 0, node->width);
  // bottom
  right_side -= get_for_edge(node->bottom_edge, alpha3, 0, -1, node->height);
  // top
  right_side += get_for_edge(node->top_edge, alpha4, 0, 1, node->height);

  return T(i, k) + dt / rho(i, k) / c(i, k) * right_side;
}

long double ExplicitSolver::T(int i, int k) {
  return previous_temp_[i][k];
}
