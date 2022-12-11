#include "column_iteration_solver.h"

#include <cassert>

ColumnIterationSolver::ColumnIterationSolver(PropertiesManager *properties)
    : PropertiesWrapper(properties), N_(properties->GetGridWidth()),
      M_(properties->GetGridHeight()), tridiagonal_(M_ + 2),
      next_(N_ + 2, M_ + 2) {}

Matrix ColumnIterationSolver::CalculateNextIteration(const Matrix &prev_iter,
                                              const Matrix &semi_prev_iter) {
  UpdateTemperatureLambda(prev_iter);

  for (int i = 0; i < N_ + 2; ++i) {
    auto column = next_.GetRow(i);
    BuildTridiagonal(prev_iter, semi_prev_iter, column, i);
    assert(TridiagonalAlgorithm(tridiagonal_, column));
  }

  return next_;
}

void ColumnIterationSolver::BuildTridiagonal(const Matrix &prev_iter,
                                             const Matrix &semi_prev_iter,
                                             Vector &column, int i) {
  if (i == 0) {
    LeftColumn(prev_iter, semi_prev_iter, column);
  } else if (i == N_ + 1) {
    RightColumn(prev_iter, semi_prev_iter, column);
  } else {
    MiddleColumn(prev_iter, semi_prev_iter, column, i);
  }
}

void ColumnIterationSolver::LeftColumn(const Matrix &prev_iter,
                                       const Matrix &semi_prev_iter,
                                       Vector &column) {
  long double R;
  long double Lambda_z;
  // (0, 0) node
  {
    R = c(0, 0) * rho(0, 0) +
        dt / 0.25 *
            (lambda(0, /*0.25*/ 0.5) / (0.5 * dz(1) * dz(1)) + alpha3 / dz(1) +
             alpha1 / dx(1));

    Lambda_z = lambda(0, /*0.25*/ 0.5) / (0.5 * dz(1)) *
               (prev_iter(0, 1) - prev_iter(0, 0));

    tridiagonal_[0][0] = 1;
    tridiagonal_[0][1] =
        -dt / (0.25 * dz(1)) * lambda(0, /*0.25*/ 0.5) / (0.5 * dz(1)) / R;
    tridiagonal_[0][2] = 0;

    column[0] =
        (c(0, 0) * rho(0, 0) * semi_prev_iter(0, 0) -
         dt / (0.25 * dz(1)) * Lambda_z +
         semi_prev_iter(0, 0) * dt / 0.25 * (alpha3 / dz(1) + alpha1 / dx(1))) /
        R;
  }

  // (0, M + 1) node
  {
    R = c(0, M_) * rho(0, M_) +
        dt / 0.25 *
            (lambda(0, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_) * dz(M_)) +
             alpha1 / dx(1) + alpha4 / dz(M_));

    Lambda_z = lambda(0, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_)) *
               (prev_iter(0, M_ + 1) - prev_iter(0, M_));

    tridiagonal_[M_ + 1][1] = 1;
    tridiagonal_[M_ + 1][0] = -dt / (0.25 * dz(M_)) *
                              lambda(0, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_)) / R;
    tridiagonal_[M_ + 1][2] = 0;

    column[M_ + 1] = (c(0, M_) * rho(0, M_) * semi_prev_iter(0, M_ + 1) +
                      dt / (0.25 * dz(M_)) * Lambda_z +
                      semi_prev_iter(0, M_ + 1) * dt / 0.25 *
                          (alpha4 / dz(M_) + alpha1 / dx(1))) /
                     R;
  }

  // (0, i), i = 1..M, nodes
  for (int k = 1; k <= M_; ++k) {
    tridiagonal_[k][0] = dt * lambda(0, k - 0.5) / (dz(k) * dzb(k));
    tridiagonal_[k][2] = dt * lambda(0, k + 0.5) / (dz(k) * dzb(k + 1));

    tridiagonal_[k][1] = c(0, k) * rho(0, k) + tridiagonal_[k][0] +
                         tridiagonal_[k][2] + dt * alpha1 / (0.25 * dx(1));
    tridiagonal_[k][1] *= -1;

    Lambda_z = (lambda(0, k + 0.5) / dzb(k + 1) *
                    (prev_iter(0, k + 1) - prev_iter(0, k)) -
                lambda(0, k - 0.5) / dzb(k) *
                    (prev_iter(0, k) - prev_iter(0, k - 1))) /
               dz(k);

    column[k] = c(0, k) * rho(0, k) * semi_prev_iter(0, k) +
                dt / (0.25 * dx(1)) * alpha1 * semi_prev_iter(0, k) -
                dt * Lambda_z;
    column[k] *= -1;
  }
}

void ColumnIterationSolver::MiddleColumn(const Matrix &prev_iter,
                                         const Matrix &semi_prev_iter,
                                         Vector &column, int i) {
  long double R;
  long double Lambda_z;

  // (i, 0) node
  {
    R = c(i, 0) * rho(i, 0) +
        dt / 0.25 *
            (lambda(i, /*0.25*/ 0.5) / (0.5 * dz(1) * dz(1)) + alpha3 / dz(1));

    Lambda_z = lambda(i, /*0.25*/ 0.5) / (0.5 * dz(1)) *
               (prev_iter(i, 1) - prev_iter(i, 0));

    tridiagonal_[0][0] = 1;
    tridiagonal_[0][1] =
        -dt / (0.25 * dz(1)) * lambda(i, /*0.25*/ 0.5) / (0.5 * dz(1)) / R;
    tridiagonal_[0][2] = 0;

    column[0] = (c(i, 0) * rho(i, 0) * semi_prev_iter(i, 0) -
                 dt / (0.25 * dz(1)) * Lambda_z +
                 dt * alpha3 * semi_prev_iter(i, 0) / (0.25 * dz(1))) /
                R;
  }

  // (N + 1, M + 1) node
  {
    // Common part for all nodes in [1, NxIb - 1]+[NxIf + 2, Nx] and (NxIb) and
    // (NxIf + 1) and [NxIb + 2, NxIf - 1] and (NxIb + 1) and (NxIf)
    Lambda_z = lambda(i, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_)) *
               (prev_iter(i, M_ + 1) - prev_iter(i, M_));
    tridiagonal_[M_ + 1][1] = 1;
    tridiagonal_[M_ + 1][0] =
        -dt * lambda(i, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_));
    tridiagonal_[M_ + 1][2] = 0;

    R = c(i, M_) * rho(i, M_);

    // Correct
    if (i <= manager_->GetToolStartI() || manager_->GetToolFinishI() + 1 <= i) {
      R += dt / (0.25 * dz(M_)) *
           (lambda(i, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_)) + alpha4);

      tridiagonal_[M_ + 1][0] /= 0.25 * dz(M_) * R;

      column[M_ + 1] =
          (c(i, M_) * rho(i, M_) * semi_prev_iter(i, M_ + 1) + /* (-) in doc ?*/
           dt / (0.25 * dz(M_)) * Lambda_z +
           semi_prev_iter(i, M_ + 1) * dt * alpha4 / (0.25 * dz(M_))) /
          R;
    }

    if (i > manager_->GetToolStartI() + 1 && i < manager_->GetToolFinishI()) {
      R += dt / manager_->GetToolWaveHeight() *
           (lambda(i, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_)) +
            lambda(i, M_) /
                (manager_->GetToolHeight() - manager_->GetToolPenetration()));
      tridiagonal_[M_ + 1][0] /= manager_->GetToolWaveHeight() * R;

      column[M_ + 1] = (c(i, M_) * rho(i, M_) * semi_prev_iter(i, M_ + 1) +
                        dt / manager_->GetToolWaveHeight() *
                            (lambda(i, M_) * semi_prev_iter(i, M_ + 1) /
                                 (manager_->GetToolHeight() -
                                  manager_->GetToolPenetration()) +
                             Lambda_z)) /
                       R;
    }

    if (i == manager_->GetToolStartI() + 1 || i == manager_->GetToolFinishI()) {
      R += dt / manager_->GetToolWaveHeight() *
           (lambda(i, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_)) +
            (manager_->GetToolWaveHeight() - 0.25 * dz(M_)) * alpha4 / dx(i) +
            lambda(i, M_) /
                (manager_->GetToolHeight() - manager_->GetToolPenetration()));

      tridiagonal_[M_ + 1][0] /= manager_->GetToolWaveHeight() * R;

      column[M_ + 1] = (c(i, M_) * rho(i, M_) * semi_prev_iter(i, M_ + 1) +
                        dt / manager_->GetToolWaveHeight() *
                            (lambda(i, M_) * semi_prev_iter(i, M_ + 1) /
                                 (manager_->GetToolHeight() -
                                  manager_->GetToolPenetration()) +
                             (manager_->GetToolWaveHeight() - 0.25 * dz(M_)) *
                                 alpha4 * semi_prev_iter(i, M_ + 1) / dx(i) +
                             Lambda_z)) /
                       R;
    }
  }

  // (i, k), k = 1..M_, nodes
  for (int k = 1; k <= M_; ++k) {
    // Probably with / 0.5 (in doc)
    tridiagonal_[k][0] = dt * lambda(i, k - 0.5) / (dz(k) * dzb(k));
    tridiagonal_[k][2] = dt * lambda(i, k + 0.5) / (dz(k) * dzb(k + 1));
    tridiagonal_[k][1] =
        c(i, k) * rho(i, k) + tridiagonal_[k][0] + tridiagonal_[k][2];
    tridiagonal_[k][1] *= -1;

    Lambda_z = (lambda(i, k + 0.5) / dzb(k + 1) *
                    (prev_iter(i, k + 1) - prev_iter(i, k)) -
                lambda(i, k - 0.5) / dzb(k) *
                    (prev_iter(i, k) - prev_iter(i, k - 1))) /
               dz(k);

    column[k] = c(i, k) * rho(i, k) * semi_prev_iter(i, k) - dt * Lambda_z;
    column[k] *= -1;
  }
}

void ColumnIterationSolver::RightColumn(const Matrix &prev_iter,
                                        const Matrix &semi_prev_iter,
                                        Vector &column) {
  long double R;
  long double Lambda_z;

  // (N + 1, 0) node
  {
    R = c(N_, 0) * rho(N_, 0) +
        dt / 0.25 *
            (lambda(N_, /*0.25*/ 0.5) / (0.5 * dz(1) * dz(1)) + alpha3 / dz(1) +
             alpha2 / dx(N_));
    Lambda_z = lambda(N_, /*0.25*/ 0.5) / (0.5 * dz(1)) *
               (prev_iter(N_ + 1, 1) - prev_iter(N_ + 1, 0));
    tridiagonal_[0][0] = 1;
    tridiagonal_[0][1] =
        -dt / (0.25 * dz(1)) * lambda(N_, /*0.25*/ 0.5) / (0.5 * dz(1)) / R;
    tridiagonal_[0][2] = 0;

    column[0] = (c(N_, 0) * rho(N_, 0) * semi_prev_iter(N_ + 1, 0) -
                 dt / (0.25 * dz(1)) * Lambda_z +
                 semi_prev_iter(N_ + 1, 0) * dt / 0.25 *
                     (alpha3 / dz(1) + alpha2 / dx(N_))) /
                R;
  }

  // (N + 1, M + 1) node
  {
    R = c(N_, M_) * rho(N_, M_) +
        dt / 0.25 *
            (lambda(N_, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_) * dz(M_)) +
             alpha2 / dx(N_) + alpha4 / dz(M_));
    Lambda_z = lambda(N_, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_)) *
               (prev_iter(N_ + 1, M_ + 1) - prev_iter(N_ + 1, M_));

    tridiagonal_[M_ + 1][1] = 1;
    tridiagonal_[M_ + 1][0] = -dt / (0.25 * dz(M_)) *
                              lambda(N_, M_ - /*0.25*/ 0.5) / (0.5 * dz(M_)) /
                              R;
    tridiagonal_[M_ + 1][2] = 0;

    column[M_ + 1] = (c(N_, M_) * rho(N_, M_) *
                          semi_prev_iter(N_ + 1, M_ + 1) + /* (-) in doc ?*/
                      dt / (0.25 * dz(M_)) * Lambda_z +
                      semi_prev_iter(N_ + 1, M_ + 1) * dt / 0.25 *
                          (alpha4 / dz(M_) + alpha2 / dx(N_))) /
                     R;
  }

  // (N + 1, i), i = 1..M, nodes
  for (int k = 1; k <= M_; ++k) {
    tridiagonal_[k][0] = dt * lambda(N_, k - 0.5) / (dz(k) * dzb(k));
    tridiagonal_[k][2] = dt * lambda(N_, k + 0.5) / (dz(k) * dzb(k + 1));
    tridiagonal_[k][1] = c(N_, k) * rho(N_, k) + tridiagonal_[k][0] +
                         tridiagonal_[k][2] + dt * alpha2 / (0.25 * dx(N_));
    tridiagonal_[k][1] *= -1;

    Lambda_z = (lambda(N_, k + 0.5 /*(k + 1) in dock*/) / dzb(k + 1) *
                    (prev_iter(N_ + 1, k + 1) - prev_iter(N_ + 1, k)) -
                lambda(N_, k - 0.5 /*(k) in dock*/) / dzb(k) *
                    (prev_iter(N_ + 1, k) - prev_iter(N_ + 1, k - 1))) /
               dz(k);

    column[k] = c(N_, k) * rho(N_, k) * semi_prev_iter(N_ + 1, k) +
                dt / (0.25 * dx(N_)) * alpha2 * semi_prev_iter(N_ + 1, k) -
                dt * Lambda_z;
    column[k] *= -1;
  }
}
