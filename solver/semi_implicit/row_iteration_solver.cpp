#include "row_iteration_solver.h"

#include <cassert>

RowIterationSolver::RowIterationSolver(PropertiesManager *properties)
    : PropertiesWrapper(properties), N_(properties->GetGridWidth()),
      M_(properties->GetGridHeight()), tridiagonal_(N_ + 2),
      next_(N_ + 2, M_ + 2) {}

Matrix RowIterationSolver::CalculateNextIteration(const Matrix &prev_iter,
                                                  const Matrix &prev_layer) {
  UpdateTemperatureLambda(prev_layer);

  for (int k = 0; k < M_ + 2; ++k) {
    auto row = next_.GetColumn(k);
    BuildTridiagonal(prev_iter, prev_layer, row, k);
    assert(TridiagonalAlgorithm(tridiagonal_, row));
  }

  return next_;
}

void RowIterationSolver::BuildTridiagonal(const Matrix &prev_iter,
                                          const Matrix &prev_layer, Vector &row,
                                          int k) {
  if (k == 0) {
    BottomRow(prev_iter, prev_layer, row);
  } else if (k == M_ + 1) {
    TopRow(prev_iter, prev_layer, row);
  } else {
    MiddleRow(prev_iter, prev_layer, row, k);
  }
}

void RowIterationSolver::BottomRow(const Matrix &prev_iter,
                                   const Matrix &prev_layer, Vector &row) {
  auto Lambda = [&](int i) -> long double {
    int lambda_ind = i == N_ + 1 ? N_ : i;
    return lambda(lambda_ind, /*0.25*/ 0.5) / (0.5 * dz(1)) *
           (prev_iter(i, 1) - prev_iter(i, 0));
  };

  long double R;

  // (0, 0) node
  {
    R = c(0, 0) * rho(0, 0) +
        dt / (0.25 * dx(1)) * (lambda(0.5, 0) / (0.5 * dx(1)) + alpha1) +
        dt / (0.25 * dz(1)) * alpha3;

    tridiagonal_[0][0] = 1;
    tridiagonal_[0][1] = dt / (0.25 * dx(1)) * lambda(0.5, 0) / (0.5 * dx(1));
    tridiagonal_[0][1] /= -R;
    tridiagonal_[0][2] = 0;

    row[0] = c(0, 0) * rho(0, 0) * prev_layer(0, 0) +
             dt * t_out / 0.25 * (alpha1 / dx(1) + alpha3 / dz(1)) +
             dt * Lambda(0) / (0.25 * dz(1));
    row[0] /= R;
  }

  // (N + 1, 0) node
  {
    R = c(N_, 0) * rho(N_, 0) +
        dt / (0.25 * dx(N_)) * (lambda(N_ - 0.5, 0) / (0.5 * dx(N_)) + alpha2) +
        dt / (0.25 * dz(1)) * alpha3;

    tridiagonal_[N_ + 1][1] = 1;
    tridiagonal_[N_ + 1][0] =
        dt / (0.25 * dx(N_)) * lambda(N_ - 0.5, 0) / (0.5 * dx(N_));
    tridiagonal_[N_ + 1][0] /= -R;
    tridiagonal_[N_ + 1][2] = 0;

    row[N_ + 1] = c(N_, 0) * rho(N_, 0) * prev_layer(N_ + 1, 0) +
                  dt * t_out / 0.25 * (alpha2 / dx(N_) + alpha3 / dz(1)) +
                  dt * Lambda(N_ + 1) / (0.25 * dz(1));
    row[N_ + 1] /= R;
  }

  // (i, 0) node, i = 1..N
  for (int i = 1; i <= N_; ++i) {
    tridiagonal_[i][0] = dt * lambda(i - 0.5, 0) / (0.5 * dx(i) * dxb(i));
    tridiagonal_[i][2] = dt * lambda(i + 0.5, 0) / (0.5 * dx(i) * dxb(i + 1));
    tridiagonal_[i][1] = c(i, 0) * rho(i, 0) + dt * alpha3 / (0.25 * dz(1));
    tridiagonal_[i][1] += tridiagonal_[i][0] + tridiagonal_[i][2];
    tridiagonal_[i][1] *= -1;

    row[i] = c(i, 0) * rho(i, 0) * prev_layer(i, 0) +
             dt / (0.25 * dz(1)) * (alpha3 * t_out + Lambda(i));
    row[i] *= -1;
  }
}

void RowIterationSolver::MiddleRow(const Matrix &prev_iter,
                                   const Matrix &prev_layer, Vector &row,
                                   int k) {
  // For (0, k) && (N + 1, k) doesn't match with same shit in doc
  auto Lambda_zz = [&](int i) {
    int lambda_x = i == N_ + 1 ? N_ : i;
    double lambda_z = i == 0 || i == N_ + 1 ? k + 1 : k + 0.5;
    // Probably with / 0.5
    return (lambda(lambda_x, lambda_z) *
                (prev_iter(i, k + 1) - prev_iter(i, k)) / dzb(k + 1) -
            lambda(lambda_x, lambda_z - 1) *
                (prev_iter(i, k) - prev_iter(i, k - 1)) / dzb(k)) /
           dz(k);
  };

  long double R;
  // (0, k) node
  {
    R = c(0, k) * rho(0, k) +
        dt / (0.25 * dx(1)) * (alpha1 + lambda(0.5, k) / (0.5 * dx(1)));

    tridiagonal_[0][0] = 1;
    tridiagonal_[0][1] = dt / (0.25 * dx(1)) * lambda(0.5, k) / (0.5 * dx(1));
    tridiagonal_[0][1] /= -R;
    tridiagonal_[0][2] = 0;

    row[0] = c(0, k) * rho(0, k) * prev_layer(0, k) +
             dt * (alpha1 * t_out / (0.25 * dx(1)) + Lambda_zz(0));
    row[0] /= R;
  }

  // (N + 1, k) node
  {
    R = c(N_, k) * rho(N_, k) +
        dt / (0.25 * dx(N_)) * (alpha2 + lambda(N_ - 0.5, k) / (0.5 * dx(N_)));

    tridiagonal_[N_ + 1][1] = 1;
    tridiagonal_[N_ + 1][0] =
        dt / (0.25 * dx(N_)) * lambda(N_ - 0.5, k) / (0.5 * dx(N_));
    tridiagonal_[N_ + 1][0] /= -R;
    tridiagonal_[N_ + 1][2] = 0;

    row[N_ + 1] = c(N_, k) * rho(N_, k) * prev_layer(N_ + 1, k) +
                  dt * (alpha2 * t_out / (0.25 * dx(N_)) + Lambda_zz(N_ + 1));
    row[N_ + 1] /= R;
  }

  // (i, k), i = 1..N, nodes
  for (int i = 1; i <= N_; ++i) {
    // Probably without / 0.5
    tridiagonal_[i][0] = dt * lambda(i - 0.5, k) / (0.5 * dx(i) * dxb(i));
    tridiagonal_[i][2] = dt * lambda(i + 0.5, k) / (0.5 * dx(i) * dxb(i + 1));

    tridiagonal_[i][1] = c(i, k) * rho(i, k);
    tridiagonal_[i][1] += tridiagonal_[i][0] + tridiagonal_[i][2];
    tridiagonal_[i][1] *= -1;

    row[i] = c(i, k) * rho(i, k) * prev_layer(i, k) +
             dt * (Lambda_zz(i) + manager_->GetHeatX(i, k) +
                   manager_->GetHeatZ(i, k));
    row[i] *= -1;
  }
}

void RowIterationSolver::TopRow(const Matrix &prev_iter,
                                const Matrix &prev_layer, Vector &row) {
  auto Lambda_z = [&](int i) -> long double {
    // according to doc
    int lambda_ind = i == N_ + 1 ? N_ : i;
    return lambda(lambda_ind, M_ - /*0.25*/ 0.5) *
           (prev_iter(i, M_ + 1) - prev_iter(i, M_)) / (0.5 * dz(M_));
  };

  long double R;
  // (0, M + 1) node
  {
    R = c(0, M_) * rho(0, M_) +
        dt / (0.25 * dx(1)) * (alpha1 + lambda(0.5, M_) / (0.5 * dx(1))) +
        dt / (0.25 * dz(M_)) * alpha4;

    tridiagonal_[0][0] = 1;
    tridiagonal_[0][1] = dt / (0.25 * dx(1)) * lambda(0.5, M_) / (0.5 * dx(1));
    tridiagonal_[0][1] /= -R;
    tridiagonal_[0][2] = 0;

    row[0] = c(0, M_) * rho(0, M_) * prev_layer(0, M_ + 1) +
             dt / (0.25 * dz(M_)) * (alpha4 * t_out - Lambda_z(0)) +
             dt / (0.25 * dx(1)) * alpha1 * t_out;
    row[0] /= R;
  }

  // (N + 1, M + 1) node
  {
    R = c(N_, M_) * rho(N_, M_) +
        dt / (0.25 * dx(N_)) *
            (alpha2 + lambda(N_ - 0.5, M_) / (0.5 * dx(N_))) +
        dt / (0.25 * dz(M_)) * alpha4;

    tridiagonal_[N_ + 1][1] = 1;
    tridiagonal_[N_ + 1][0] =
        dt / (0.25 * dx(N_)) * lambda(N_ - 0.5, M_) / (0.5 * dx(N_));
    tridiagonal_[N_ + 1][0] /= -R;
    tridiagonal_[N_ + 1][2] = 0;

    row[N_ + 1] = c(N_, M_) * rho(N_, M_) * prev_layer(N_ + 1, M_ + 1) +
                  dt / (0.25 * dz(M_)) * (alpha3 * t_out - Lambda_z(N_ + 1)) +
                  dt / (0.25 * dx(N_)) * alpha2 * t_out;
    row[N_ + 1] /= R;
  }

  // (i, M + 1), i = 1..(Nx * Ib - 1) + (Nx * If + 2)..N, nodes and part of (Nx
  // * Ib, M + 1) and (Nx * If + 1, M + 1) nodes
  for (int i = 1; i <= N_; ++i) {
    if (i == manager_->GetToolStartI() + 1) {
      i = manager_->GetToolFinishI();
      continue;
    }

    tridiagonal_[i][0] = dt * lambda(i - 0.5, M_) / (0.5 * dx(i) * dxb(i));
    tridiagonal_[i][2] = dt * lambda(i + 0.5, M_) / (0.5 * dx(i) * dxb(i + 1));
    tridiagonal_[i][1] = c(i, M_) * rho(i, M_) + dt * alpha4 / (0.25 * dz(M_));
    tridiagonal_[i][1] += tridiagonal_[i][0] + tridiagonal_[i][2];
    tridiagonal_[i][1] *= -1;

    row[i] = c(i, M_) * rho(i, M_) * prev_layer(i, M_ + 1) +
             dt / (0.25 * dz(M_)) * (-Lambda_z(i) + alpha4 * t_out);
    row[i] *= -1;
  }

  // correct (Nx * Ib, M + 1) node
  {
    int i = manager_->GetToolStartI();
    row[i] -= dt * manager_->GetHeatX(i, M_);
  }

  // correct (Nx * If + 1, M + 1) node
  {
    int i = manager_->GetToolFinishI() + 1;
    row[i] -= dt * manager_->GetHeatX(i, M_);
  }

  // (i, M + 1), i = (Nx * Ib + 2)..(Nx * If - 1), nodes and part of (Nx * Ib +
  // 1, M + 1) and (Nx * If, M + 1) nodes
  for (int i = manager_->GetToolStartI() + 1; i <= manager_->GetToolFinishI();
       ++i) {
    tridiagonal_[i][0] = dt * lambda(i - 0.5, M_) / (0.5 * dx(i) * dxb(i));
    tridiagonal_[i][2] = dt * lambda(i + 0.5, M_) / (0.5 * dx(i) * dxb(i + 1));

    tridiagonal_[i][1] =
        c(i, M_) * rho(i, M_) +
        dt / manager_->GetToolWaveHeight() * lambda(i, M_) /
            (manager_->GetToolHeight() - manager_->GetToolPenetration());
    tridiagonal_[i][1] += tridiagonal_[i][0] + tridiagonal_[i][2];
    tridiagonal_[i][1] *= -1;

    row[i] =
        c(i, M_) * rho(i, M_) * prev_layer(i, M_ + 1) +
        dt / manager_->GetToolWaveHeight() *
            (lambda(i, M_) * manager_->GetToolInitTemperature() /
                 (manager_->GetToolHeight() - manager_->GetToolPenetration()) -
             Lambda_z(i));
    row[i] *= -1;
  }

  // correct (Nx * Ib + 1, M + 1) node
  {
    int i = manager_->GetToolStartI() + 1;
    tridiagonal_[i][1] += tridiagonal_[i][0];
    tridiagonal_[i][0] *= 0.25 * dz(M_) / manager_->GetToolWaveHeight();
    tridiagonal_[i][1] -= tridiagonal_[i][0];

    tridiagonal_[i][1] -= dt / dx(i) *
                          (manager_->GetToolWaveHeight() - 0.25 * dz(M_)) /
                          manager_->GetToolWaveHeight() * alpha4;

    row[i] -=
        dt / (dx(i) * manager_->GetToolWaveHeight()) *
        (0.25 * dz(M_) * manager_->GetHeatOutputX() +
         (manager_->GetToolWaveHeight() - 0.25 * dz(M_)) * alpha4 * t_out);
  }

  // correct (Nx * If, M + 1) node
  {
    int i = manager_->GetToolFinishI();
    tridiagonal_[i][1] += tridiagonal_[i][2];
    tridiagonal_[i][2] *= 0.25 * dz(M_) / manager_->GetToolWaveHeight();
    tridiagonal_[i][1] -= tridiagonal_[i][2];

    tridiagonal_[i][1] -= dt / dx(i) *
                          (manager_->GetToolWaveHeight() - 0.25 * dz(M_)) /
                          manager_->GetToolWaveHeight() * alpha4;

    row[i] -=
        dt / (dx(i) * manager_->GetToolWaveHeight()) *
        (0.25 * dz(M_) * manager_->GetHeatOutputX() +
         (manager_->GetToolWaveHeight() - 0.25 * dz(M_)) * alpha4 * t_out);
  }
}
