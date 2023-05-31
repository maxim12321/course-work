#include "explicit_solver.h"

#include <chrono>
#include <iostream>
#include <utility>
#include <mpi.h>

ExplicitSolver::ExplicitSolver(int p_rank,
                               int p_size,
                               PropertiesManager* properties,
                               std::string result_file_name,
                               int num_threads)
    : SolverBase(p_rank, p_size, properties, [](auto) {}),
      PropertiesWrapper(properties),
      rows_per_process_((nx_ + 2 + p_size - 1) / p_size),
      row_begin_(p_rank * rows_per_process_),
      row_end_((p_rank + 1) * rows_per_process_),
      result_file_name_(std::move(result_file_name)),
      num_threads_(num_threads),
      read_requests_(2, MPI_REQUEST_NULL){
  if (p_rank == 0) {
    std::cout << p_size << " processes, " << num_threads << " threads\n";
    std::cout << properties->GetTimeLayers() << " time layers\n";
    std::cout << "Grid: " << properties->GetGridWidth() << "x" << properties->GetGridHeight() << "\n";
  }
}

void ExplicitSolver::Solve() {
  ResultSaver result_saver(result_file_name_);

  result_saver.GetStream() << (properties_->GetTimeLayers() / kLoggingSkip) + 1 << "\n";

  PrepareNodeEdges();

  auto start = std::chrono::steady_clock::now();

  result_saver.Save(current_temp_, row_begin_, row_end_);

  for (size_t i = 1; i <= properties_->GetTimeLayers(); ++i) {
    CalculateNextLayer();
    if (i % kLoggingSkip == 0) {
      result_saver.Save(current_temp_, row_begin_, row_end_);
    }
  }


  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

  MPI_Waitall(write_requests_.size(), write_requests_.data(), MPI_STATUSES_IGNORE);
}

void ExplicitSolver::CalculateNextLayer() {
  UpdateTemperatureLambda(previous_temp_);

  if (p_rank_ > 0) {
    MPI_Request request;
    MPI_Isend(previous_temp_[row_begin_].Data(), nz_ + 2, MPI_LONG_DOUBLE, p_rank_ - 1,
              101, MPI_COMM_WORLD, &request);
    write_requests_.push_back(request);

    MPI_Irecv(previous_temp_[row_begin_ - 1].Data(), nz_ + 2, MPI_LONG_DOUBLE, p_rank_ - 1,
              100, MPI_COMM_WORLD, &read_requests_[0]);
  }
  if (p_rank_ + 1 < p_size_) {
    MPI_Request request;
    MPI_Isend(previous_temp_[row_end_ - 1].Data(), nz_ + 2, MPI_LONG_DOUBLE, p_rank_ + 1,
              100, MPI_COMM_WORLD, &request);
    write_requests_.push_back(request);

    MPI_Irecv(previous_temp_[row_end_].Data(), nz_ + 2, MPI_LONG_DOUBLE, p_rank_ + 1,
              101, MPI_COMM_WORLD, &read_requests_[1]);
  }

  int row_end = std::min(row_end_, nx_ + 2);

#pragma omp parallel for num_threads(num_threads_)
  for (int i = row_begin_ + 1; i < row_end - 1; ++i) {
    for (int k = 0; k <= nz_ + 1; ++k) {
      current_temp_[i][k] = GetNodeValue(&nodes[i][k]);
    }
  }

  MPI_Waitall(read_requests_.size(), read_requests_.data(), MPI_STATUSES_IGNORE);

  for (int k = 0; k <= nz_ + 1; ++k) {
    current_temp_[row_begin_][k] = GetNodeValue(&nodes[row_begin_][k]);
  }
  for (int k = 0; k <= nz_ + 1; ++k) {
    current_temp_[row_end - 1][k] = GetNodeValue(&nodes[row_end - 1][k]);
  }

  std::swap(previous_temp_, current_temp_);
}

void ExplicitSolver::PrepareNodeEdges() {
  const int N = nx_;
  const int M = nz_;

  const int B = properties_->GetBackingStartI();

  nodes.resize(N + 2);

  for (int i = std::max(0, row_begin_ - 1); i <= std::min(row_end_, N + 1); i++) {
    nodes[i].resize(M + 2);

    for (int k = 0; k <= M + 1; k++) {
      NodeEdgeInfo& node = nodes[i][k];
      node.i = i;
      node.k = k;
      node.width = dx(i);
      node.height = dz(k);

      if (i == 0 || i == N + 1) {
        node.width /= 2;
      }

      if (k == 0 || k == B + 1) {
        node.height /= 2;
      }

      if (i == 1 || i == N) {
        node.width *= 0.75;
      }
      if (k == 1 || k == B) {
        node.height *= 0.75;
      }

      if (i == 0) {
        node.left_edge = EdgeType::kAir;
        node.right_edge = (k == 0 || k == B + 1) ? EdgeType::kMixed : EdgeType::kMaterial;
        node.bottom_edge = (k == 0) ? EdgeType::kAir : EdgeType::kMixed;
        node.top_edge = (k == B + 1) ? EdgeType::kAir : EdgeType::kMixed;
      } else if (i == N + 1) {
        node.left_edge = (k == 0 || k == B + 1) ? EdgeType::kMixed : EdgeType::kMaterial;
        node.right_edge = EdgeType::kAir;
        node.bottom_edge = (k == 0) ? EdgeType::kAir : EdgeType::kMixed;
        node.top_edge = (k == B + 1) ? EdgeType::kAir : EdgeType::kMixed;
      } else if (k == 0) {
        node.left_edge = EdgeType::kMixed;
        node.right_edge = EdgeType::kMixed;
        node.top_edge = EdgeType::kMaterial;
        node.bottom_edge = EdgeType::kAir;
      } else if (k == B + 1) {
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

      if (i > properties_->GetToolStartI() && i <= properties_->GetToolFinishI()) {
        if (k == B + 1) {
          node.top_edge = EdgeType::kMaterial;
        } else if (k == M + 1) {
          node.left_edge = EdgeType::kMixed;
          node.right_edge = EdgeType::kMixed;
          node.top_edge = EdgeType::kAir;
        }
      } else {
        if (k >= B + 2) {
          node.left_edge = EdgeType::kNone;
          node.right_edge = EdgeType::kNone;
          node.top_edge = EdgeType::kNone;
          node.bottom_edge = EdgeType::kNone;
        }
      }

      if (i == properties_->GetToolStartI() + 1 && k > B + 1 && k <= M + 1) {
        node.left_edge = EdgeType::kAir;
        if (k < M + 1) {
          node.top_edge = EdgeType::kMixed;
        }
        if (k > B + 2) {
          node.bottom_edge = EdgeType::kMixed;
        }
      }

      if (i == properties_->GetToolFinishI() && k > B + 1 && k <= M + 1) {
        node.right_edge = EdgeType::kAir;
        if (k < M + 1) {
          node.top_edge = EdgeType::kMixed;
        }
        if (k > B + 2) {
          node.bottom_edge = EdgeType::kMixed;
        }
      }
    }
  }
}

long double ExplicitSolver::GetNodeValue(NodeEdgeInfo* node) {
  int i = node->i;
  int k = node->k;

  if (node->top_edge == EdgeType::kNone) {
    return T(i, k);
  }

  long double right_side = manager_->GetHeatX(i, k) + manager_->GetHeatZ(i, k);

  auto get_air = [&](double alpha, long double size) -> long double {
    return -alpha * (T(i, k) - t_out) / size;
  };

  auto get_material = [&](int dx, int dk) -> long double {
    long double lbd = lambda(i + 0.5 * dx, k + 0.5 * dk);

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
          case EdgeType::kAir:return get_air(alpha, size);
          case EdgeType::kMaterial:return get_material(dx, dk);
          case EdgeType::kMixed:return 0.5 * get_air(alpha, size) + 0.5 * get_material(dx, dk);
          case EdgeType::kNone:return 0;
        }
      };

  double left_alpha = alpha1;
  double right_alpha = alpha2;
  if (k > properties_->GetBackingStartI()) {
    left_alpha = properties_->GetAlpha4Tool();
    right_alpha = properties_->GetAlpha4Tool();
  }

  // left
  right_side += get_for_edge(node->left_edge, left_alpha, -1, 0, node->width);
  // right
  right_side += get_for_edge(node->right_edge, right_alpha, 1, 0, node->width);
  // bottom
  right_side += get_for_edge(node->bottom_edge, alpha3, 0, -1, node->height);
  // top
  right_side += get_for_edge(node->top_edge, alpha4, 0, 1, node->height);

  return T(i, k) + dt / rho(i, k) / c(i, k) * right_side;
}

long double ExplicitSolver::T(int i, int k) {
  return previous_temp_[i][k];
}
