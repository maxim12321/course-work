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
  auto start = std::chrono::steady_clock::now();

  current_temp_.Store(output_, row_begin_, row_end_);

  for (size_t i = 0; i < properties_->GetTimeLayers(); ++i) {
    CalculateNextLayer();
    if (i % 100 == 0) {
      // Print some debug info for long calculations
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
      current_temp_[i][k] = GetNodeValue(i, k);
    }
  }

  previous_temp_ = current_temp_;
  current_temp_.Store(output_, row_begin_, row_end_);
}

long double ExplicitSolver::GetNodeValue(int i, int k) {
  const int N = nx_;
  const int M = nz_;

  // GetToolWaveHeight() - 0.25 dz(M)
  const double height_diff = manager_->GetToolHeight() - manager_->GetToolPenetration();

  // i=0, 0<=k<=M+1
  auto lambda_x_0k = [&](int k) -> long double {
    // TODO: maybe k+0.5, as in lambda_x_Nk
    int lambda_z = std::min(M, k);
    return lambda(0.25, lambda_z) * (T(1, k) - T(0, k)) / (0.5 * dx(1));
  };

  // i=N+1, 0<=k<=M+1
  auto lambda_x_Nk = [&](int k) -> long double {
    double lambda_z = k + 0.5;
    if (k == 0) {
      lambda_z = 0;
    } else if (k == M + 1) {
      lambda_z = M;
    }
    return lambda(N - 0.25, lambda_z) * (T(N + 1, k) - T(N, k)) / (0.5 * dx(N));
  };

  // i in tool, k=M+1
  auto lambda_x_iM = [&](int i) -> long double {
    return lambda(i + 0.5, M) * (T(i + 1, M + 1) - T(i, M + 1)) / (dxb(i + 1));
  };

  // i in tool, k=M+1
  auto lambda_xb_iM = [&](int i) -> long double {
    return lambda(i - 0.5, M) * (T(i, M + 1) - T(i - 1, M + 1)) / (dxb(i));
  };

  // 0<=i<=N+1, k=0
  auto lambda_z_i0 = [&](int i) -> long double {
    int lambda_x = std::min(N, i);
    return lambda(lambda_x, 0.25) * (T(i, 1) - T(i, 0)) / (0.5 * dz(1));
  };

  // 0<=i<=M+1, k=M+1
  auto lambda_z_iM = [&](int i) -> long double {
    int lambda_x = std::min(N, i);
    return lambda(lambda_x, M - 0.25) * (T(i, M + 1) - T(i, M)) / (0.5 * dz(M));
  };

  // 1<=i<=N, k=0
  auto lambda_xx_i0 = [&](int i) -> long double {
    return (
        lambda(i + 0.5, 0) * (T(i + 1, 0) - T(i, 0)) / (0.5 * dxb(i + 1)) -
            lambda(i - 0.5, 0) * (T(i, 0) - T(i - 1, 0)) / (0.5 * dxb(i))
    ) / dx(i);
  };

  // 1<=i<=N, k=M+1
  auto lambda_xx_iM = [&](int i) -> long double {
    return (
        lambda(i + 0.5, M) * (T(i + 1, M + 1) - T(i, M + 1)) / dxb(i + 1) -
            lambda(i - 0.5, M) * (T(i, M + 1) - T(i - 1, M + 1)) / dxb(i)
    ) / dx(i);
  };

  auto lambda_xx_ik = [&](int i, int k) -> long double {
    return (
        lambda(i + 0.5, k) * (T(i + 1, k) - T(i, k)) / (0.5 * dxb(i + 1)) -
            lambda(i - 0.5, k) * (T(i, k) - T(i - 1, k)) / (0.5 * dxb(i))
    ) / dx(i);
  };

  // i=0 or i=N+1, 1<=k<=M
  auto lambda_zz_0k = [&](int i, int k) -> long double {
    int lambda_x = std::min(N, i);
    /// lambda(x, k+1) - lambda(x, k) in doc
    return (
        lambda(lambda_x, k + 0.5) * (T(i, k + 1) - T(i, k)) / dzb(k + 1) -
            lambda(lambda_x, k - 0.5) * (T(i, k) - T(i, k - 1)) / dzb(k)
    ) / dz(k);
  };

  auto lambda_zz_tool = [&](int i) -> long double {
    return (
        lambda(i, M) * (manager_->GetToolInitTemperature() - T(i, M + 1)) / height_diff -
            lambda(i, M - 0.25) * (T(i, M + 1) - T(i, M)) / (0.5 * dz(M))
    ) / manager_->GetToolWaveHeight();
  };

  auto lambda_zz_ik = [&](int i, int k) -> long double {
    return (
        lambda(i, k + 0.5) * (T(i, k + 1) - T(i, k)) / dzb(k + 1) -
            lambda(i, k - 0.5) * (T(i, k) - T(i, k - 1)) / dzb(k)
    ) / dz(i);
  };

  long double right_side = 0;

  if (i == 0 && k == 0) {
    right_side = lambda_x_0k(k) / (0.25 * dx(1)) +
        lambda_z_i0(0) / (0.25 * dz(1)) -
        alpha1 * (T(0, 0) - t_out) / (0.25 * dx(1)) -
        alpha3 * (T(0, 0) - t_out) / (0.25 * dz(1));
  }

  if (i == 0 && k >= 1 && k <= M) {
    right_side = lambda_x_0k(k) / (0.25 * dx(1)) +
        lambda_zz_0k(0, k) -
        alpha1 * (T(0, k) - t_out) / (0.25 * dx(1));
  }

  if (i == 0 && k == M + 1) {
    right_side = lambda_x_0k(k) / (0.25 * dx(1)) -
        lambda_z_iM(0) / (0.25 * dz(M)) -
        alpha1 * (T(0, M + 1) - t_out) / (0.25 * dx(1)) -
        alpha4 * (T(0, M + 1) - t_out) / (0.25 * dz(M));
  }

  if (i == N + 1 && k == 0) {
    right_side = -lambda_x_Nk(k) / (0.25 * dx(N)) +
        lambda_z_i0(N + 1) / (0.25 * dz(1)) -
        alpha2 * (T(N + 1, 0) - t_out) / (0.25 * dx(N)) -
        alpha3 * (T(N + 1, 0) - t_out) / (0.25 * dz(1));
  }

  if (i == N + 1 && k >= 1 && k <= M) {
    right_side = -lambda_x_Nk(k) / (0.25 * dx(N)) +
        lambda_zz_0k(N + 1, k) -
        alpha2 * (T(N + 1, k) - t_out) / (0.25 * dx(N));
  }

  if (i == N + 1 && k == M + 1) {
    right_side = -lambda_x_Nk(M + 1) / (0.25 * dx(N)) -
        lambda_z_iM(N + 1) / (0.25 * dz(M)) -
        alpha2 * (T(N + 1, M + 1) - t_out) / (0.25 * dx(N)) -
        alpha3 * (T(N + 1, M + 1) - t_out) / (0.25 * dz(M));
  }

  if (i >= 1 && i <= N && k == 0) {
    right_side = lambda_xx_i0(i) +
        lambda_z_i0(i) / (0.25 * dz(1)) -
        alpha3 * (T(i, 0) - t_out) / (0.25 * dz(1));
  }

  if (k == M + 1 &&
      ((i >= 1 && i < manager_->GetToolStartI()) ||
          (i > manager_->GetToolFinishI() + 1 && i <= M)
      )) {
    right_side = lambda_xx_iM(i) -
        lambda_z_iM(i) / (0.25 * dz(M)) -
        alpha4 * (T(i, M + 1) - t_out) / (0.25 * dz(M));
  }

  if (i == manager_->GetToolStartI() && k == M + 1) {
    right_side = lambda_xx_iM(i) -
        lambda_z_iM(i) / (0.25 * dz(M)) +
        manager_->GetHeatX(i, M) -
        alpha4 * (T(i, M + 1) - t_out) / (0.25 * dz(M));
  }

  if (i == manager_->GetToolStartI() + 1 && k == M + 1) {
    right_side = lambda_x_iM(i) / dx(i) -
        0.25 * dz(M) / manager_->GetToolWaveHeight() * lambda_xb_iM(i) / dx(i) +
        lambda_zz_tool(i) +  // точно плюс???
        0.25 * dz(M) / manager_->GetToolWaveHeight() * manager_->GetHeatX(i, k) -
        height_diff / manager_->GetToolWaveHeight() * alpha4 * (T(i, M + 1) - t_out) / dx(i);
  }

  if (i > manager_->GetToolStartI() + 1 && i < manager_->GetToolFinishI() && k == M + 1) {
    right_side = lambda_xx_iM(i) + lambda_zz_tool(i);
  }

  if (i == manager_->GetToolFinishI() && k == M + 1) {
    right_side = 0.25 * dz(M) / manager_->GetToolWaveHeight() * lambda_x_iM(i) / dx(i) -
        lambda_xb_iM(i) / dx(i) +
        lambda_zz_tool(i) -  // точно минус???
        0.25 * dz(M) / manager_->GetToolWaveHeight() * manager_->GetHeatX(i, k) -
        height_diff / manager_->GetToolWaveHeight() * alpha4 * (T(i, M + 1) - t_out) / dx(i);
  }

  if (i == manager_->GetToolFinishI() + 1 && k == M + 1) {
    right_side = lambda_xx_iM(i) -
        lambda_z_iM(i) / (0.25 * dz(M)) +
        manager_->GetHeatX(i, M) -
        alpha4 * (T(i, M + 1) - t_out) / (0.25 * dz(M));
  }

  if (i >= 1 && i <= N && k >= 1 && k <= M) {
    right_side = lambda_xx_ik(i, k) +
        lambda_zz_ik(i, k) +
        manager_->GetHeatX(i, k) +
        manager_->GetHeatZ(i, k);
  }

  return T(i, k) + dt / rho(i, k) / c(i, k) * right_side;
}

long double ExplicitSolver::T(int i, int k) {
  return previous_temp_[i][k];
}
