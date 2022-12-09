#include "row_iteration_solver.h"
#include "tridiagonal.h"

#include <cassert>

RowIterationSolver::RowIterationSolver(PropertiesManager *properties)
    : properties_(properties),
      N_(properties_->GetGridWidth()),
      M_(properties_->GetGridHeight()),
      tridiagonal_(N_ + 2, 3),
      next_(M_ + 2, N_ + 2) {}

Matrix RowIterationSolver::CalculateNextIteration(const Matrix &prev_iter,
                                                  const Matrix &prev_layer) {
  for (int k = 0; k < M_ + 2; ++k) {
    BuildTridiagonal(prev_iter, prev_layer, next_[k], k);
    assert(TridiagonalAlgorithm(tridiagonal_, next_[k]));
  }

  return next_.Transposed();
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
    return properties_->GetThermalConductivity(lambda_ind, 0.5, prev_layer) /
           (0.5 * properties_->GetDeltaZ(1)) *
           (prev_iter[i][1] - prev_iter[i][0]);
    //        return properties_->GetThermalConductivity(lambda_ind, 0.25,
    //        prev_layer) / (0.5 * properties_->GetDeltaZ(1)) * (prev_iter[i][1]
    //        - prev_iter[i][0]);
  };

  // (0, 0) node
  long double R =
      properties_->GetHeatCapacity(0, 0, prev_layer) *
          properties_->GetDensity(0, 0, prev_layer) +
      properties_->GetDeltaT() / 0.25 *
          (properties_->GetThermalConductivity(/*0.25*/ 0.5, 0, prev_layer) /
               (0.5 * properties_->GetDeltaX(1) * properties_->GetDeltaX(1)) +
           properties_->GetAlpha1() / properties_->GetDeltaX(1) +
           properties_->GetAlpha3() / properties_->GetDeltaZ(1));

  tridiagonal_[0][0] = 1;
  tridiagonal_[0][1] =
      -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(1)) *
      properties_->GetThermalConductivity(/*0.25*/ 0.5, 0, prev_layer) /
      (0.5 * properties_->GetDeltaX(1)) / R;
  tridiagonal_[0][2] = 0;

  row[0] = (properties_->GetHeatCapacity(0, 0, prev_layer) *
                properties_->GetDensity(0, 0, prev_layer) * prev_layer[0][0] +
            properties_->GetDeltaT() / 0.25 *
                (properties_->GetAlpha1() * properties_->GetOutTemperature() /
                     properties_->GetDeltaX(1) +
                 properties_->GetAlpha3() * properties_->GetOutTemperature() /
                     properties_->GetDeltaZ(1) +
                 Lambda(0) / properties_->GetDeltaZ(1))) /
           R;

  // (N + 1, 0) node
  R = properties_->GetHeatCapacity(N_, 0, prev_layer) *
          properties_->GetDensity(N_, 0, prev_layer) +
      properties_->GetDeltaT() / 0.25 *
          (properties_->GetThermalConductivity(N_ - /*0.25*/ 0.5, 0,
                                               prev_layer) /
               (0.5 * properties_->GetDeltaX(N_) * properties_->GetDeltaX(N_)) +
           properties_->GetAlpha2() / properties_->GetDeltaX(N_) +
           properties_->GetAlpha3() / properties_->GetDeltaZ(1));

  tridiagonal_[N_ + 1][1] = 1;
  tridiagonal_[N_ + 1][0] =
      -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(N_)) *
      properties_->GetThermalConductivity(N_ - /*0.25*/ 0.5, 0, prev_layer) /
      (0.5 * properties_->GetDeltaX(N_)) / R;
  tridiagonal_[N_ + 1][2] = 0;

  row[N_ + 1] =
      (properties_->GetHeatCapacity(N_, 0, prev_layer) *
           properties_->GetDensity(N_, 0, prev_layer) * prev_layer[N_ + 1][0] +
       properties_->GetDeltaT() / 0.25 *
           (properties_->GetAlpha2() * properties_->GetOutTemperature() /
                properties_->GetDeltaX(N_) +
            properties_->GetAlpha3() * properties_->GetOutTemperature() /
                properties_->GetDeltaZ(1) +
            Lambda(N_ + 1) / properties_->GetDeltaZ(1))) /
      R;

  // (i, 0) node, i = 1..N
  for (int i = 1; i <= N_; ++i) {
    tridiagonal_[i][0] =
        properties_->GetDeltaT() *
        properties_->GetThermalConductivity(i - 0.5, 0, prev_layer) /
        (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaBackX(i));
    tridiagonal_[i][2] =
        properties_->GetDeltaT() *
        properties_->GetThermalConductivity(i + 0.5, 0, prev_layer) /
        (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaBackX(i + 1));
    tridiagonal_[i][1] = properties_->GetHeatCapacity(i, 0, prev_layer) *
                            properties_->GetDensity(i, 0, prev_layer) +
                        tridiagonal_[i][0] + tridiagonal_[i][2] +
                        properties_->GetDeltaT() * properties_->GetAlpha3() /
                            (0.25 * properties_->GetDeltaZ(1));
    tridiagonal_[i][1] *= -1;

    row[i] = properties_->GetHeatCapacity(i, 0, prev_layer) *
                 properties_->GetDensity(i, 0, prev_layer) * prev_layer[i][0] +
             properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(1)) *
                 (properties_->GetAlpha3() * properties_->GetOutTemperature() +
                  Lambda(i));
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
    return (properties_->GetThermalConductivity(lambda_x, lambda_z,
                                                prev_layer) *
                (prev_iter[i][k + 1] - prev_iter[i][k]) /
                properties_->GetDeltaBackZ(k + 1) -
            properties_->GetThermalConductivity(lambda_x, lambda_z - 1,
                                                prev_layer) *
                (prev_iter[i][k] - prev_iter[i][k - 1]) /
                properties_->GetDeltaBackZ(k)) /
           properties_->GetDeltaZ(k);
  };

  // (0, k) node
  long double R =
      properties_->GetHeatCapacity(0, k, prev_layer) *
          properties_->GetDensity(0, k, prev_layer) +
      properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(1)) *
          (properties_->GetThermalConductivity(/*0.25*/ 0.5, k, prev_layer) /
               (0.5 * properties_->GetDeltaX(1)) +
           properties_->GetAlpha1());

  tridiagonal_[0][0] = 1;
  tridiagonal_[0][1] =
      -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(1)) *
      properties_->GetThermalConductivity(/*0.25*/ 0.5, k, prev_layer) /
      (0.5 * properties_->GetDeltaX(1)) / R;
  tridiagonal_[0][2] = 0;

  row[0] = (properties_->GetHeatCapacity(0, k, prev_layer) *
                properties_->GetDensity(0, k, prev_layer) * prev_layer[0][k] +
            properties_->GetDeltaT() * properties_->GetAlpha1() *
                properties_->GetOutTemperature() /
                (0.25 * properties_->GetDeltaX(1)) +
            properties_->GetDeltaT() * Lambda_zz(0)) /
           R;

  // (N + 1, k) node
  R = properties_->GetHeatCapacity(N_, k, prev_layer) *
          properties_->GetDensity(N_, k, prev_layer) +
      properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(N_)) *
          (properties_->GetThermalConductivity(N_ - /*0.25*/ 0.5, k,
                                               prev_layer) /
               (0.5 * properties_->GetDeltaX(N_)) +
           properties_->GetAlpha2());

  tridiagonal_[N_ + 1][1] = 1;
  tridiagonal_[N_ + 1][0] =
      -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(N_)) *
      properties_->GetThermalConductivity(N_ - /*0.25*/ 0.5, k, prev_layer) /
      (0.5 * properties_->GetDeltaX(N_)) / R;
  tridiagonal_[N_ + 1][2] = 0;

  row[N_ + 1] =
      (properties_->GetHeatCapacity(N_, k, prev_layer) *
           properties_->GetDensity(N_, k, prev_layer) * prev_layer[N_ + 1][k] +
       properties_->GetDeltaT() * properties_->GetAlpha2() *
           properties_->GetOutTemperature() /
           (0.25 * properties_->GetDeltaX(N_)) +
       properties_->GetDeltaT() * Lambda_zz(N_ + 1)) /
      R;

  // (i, k), i = 1..N, nodes
  for (int i = 1; i <= N_; ++i) {
    // Probably without / 0.5
    tridiagonal_[i][0] =
        properties_->GetDeltaT() *
        properties_->GetThermalConductivity(i - 0.5, k, prev_layer) /
        (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaBackX(i));
    tridiagonal_[i][2] =
        properties_->GetDeltaT() *
        properties_->GetThermalConductivity(i + 0.5, k, prev_layer) /
        (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaBackX(i + 1));

    tridiagonal_[i][1] = properties_->GetHeatCapacity(i, k, prev_layer) *
                            properties_->GetDensity(i, k, prev_layer) +
                        tridiagonal_[i][0] + tridiagonal_[i][2];
    tridiagonal_[i][1] *= -1;

    row[i] =
        properties_->GetHeatCapacity(i, k, prev_layer) *
            properties_->GetDensity(i, k, prev_layer) * prev_layer[i][k] +
        properties_->GetDeltaT() * (Lambda_zz(i) + properties_->GetHeatX(i, k) +
                                    properties_->GetHeatZ(i, k));
    row[i] *= -1;
  }
}

void RowIterationSolver::TopRow(const Matrix &prev_iter,
                                const Matrix &prev_layer, Vector &row) {
  auto Lambda_z = [&](int i) -> long double {
    // according to doc
    int lambda_ind = i == N_ + 1 ? N_ : i;
    return properties_->GetThermalConductivity(lambda_ind, M_ - /*0.25*/ 0.5,
                                               prev_layer) *
           (prev_iter[i][M_ + 1] - prev_iter[i][M_]) /
           (0.5 * properties_->GetDeltaZ(M_));
  };

  // (0, M + 1) node
  long double R =
      properties_->GetHeatCapacity(0, M_, prev_layer) *
          properties_->GetDensity(0, M_, prev_layer) +
      properties_->GetDeltaT() / 0.25 *
          (properties_->GetThermalConductivity(/*0.25*/ 0.5, M_, prev_layer) /
               (0.5 * properties_->GetDeltaX(1) * properties_->GetDeltaX(1)) +
           properties_->GetAlpha1() / properties_->GetDeltaX(1) +
           properties_->GetAlpha4() / properties_->GetDeltaZ(M_));

  tridiagonal_[0][0] = 1;
  tridiagonal_[0][1] =
      -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(1)) *
      properties_->GetThermalConductivity(/*0.25*/ 0.5, M_, prev_layer) /
      (0.5 * properties_->GetDeltaX(1)) / R;
  tridiagonal_[0][2] = 0;

  row[0] =
      (properties_->GetHeatCapacity(0, M_, prev_layer) *
           properties_->GetDensity(0, M_, prev_layer) * prev_layer[0][M_ + 1] +
       properties_->GetDeltaT() / 0.25 *
           (-Lambda_z(0) / properties_->GetDeltaZ(M_) +
            properties_->GetAlpha1() * properties_->GetOutTemperature() /
                properties_->GetDeltaX(1) +
            properties_->GetAlpha4() * properties_->GetOutTemperature() /
                properties_->GetDeltaZ(M_))) /
      R;

  // (N + 1, M + 1) node
  R = properties_->GetHeatCapacity(N_, M_, prev_layer) *
          properties_->GetDensity(N_, M_, prev_layer) +
      properties_->GetDeltaT() / 0.25 *
          (properties_->GetThermalConductivity(N_ - /*0.25*/ 0.5, M_,
                                               prev_layer) /
               (0.5 * properties_->GetDeltaX(N_) * properties_->GetDeltaX(N_)) +
           properties_->GetAlpha2() / properties_->GetDeltaX(N_) +
           properties_->GetAlpha4() / properties_->GetDeltaZ(M_));

  tridiagonal_[N_ + 1][0] =
      -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(N_)) *
      properties_->GetThermalConductivity(N_ - /*0.25*/ 0.5, M_, prev_layer) /
      (0.5 * properties_->GetDeltaX(N_)) / R;
  tridiagonal_[N_ + 1][1] = 1;
  tridiagonal_[N_ + 1][2] = 0;

  row[N_ + 1] =
      (properties_->GetHeatCapacity(N_, M_, prev_layer) *
           properties_->GetDensity(N_, M_, prev_layer) *
           prev_layer[N_ + 1][M_ + 1] +
       properties_->GetDeltaT() / 0.25 *
           (-Lambda_z(N_ + 1) / properties_->GetDeltaZ(M_) +
            properties_->GetAlpha2() * properties_->GetOutTemperature() /
                properties_->GetDeltaX(N_) +
            properties_->GetAlpha4() * properties_->GetOutTemperature() /
                properties_->GetDeltaZ(M_))) /
      R;

  // (i, M + 1), i = 1..(Nx * Ib - 1) + (Nx * If + 2)..N, nodes and part of (Nx
  // * Ib, M + 1) and (Nx * If + 1, M + 1) nodes
  for (int i = 1; i <= N_; ++i) {
    if (i == properties_->GetToolStartI() + 1) {
      i = properties_->GetToolFinishI();
      continue;
    }

    tridiagonal_[i][0] =
        properties_->GetDeltaT() *
        properties_->GetThermalConductivity(i - 0.5, M_, prev_layer) /
        (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaBackX(i));
    tridiagonal_[i][2] =
        properties_->GetDeltaT() *
        properties_->GetThermalConductivity(i + 0.5, M_, prev_layer) /
        (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaBackX(i + 1));
    tridiagonal_[i][1] = properties_->GetHeatCapacity(i, M_, prev_layer) *
                            properties_->GetDensity(i, M_, prev_layer) +
                        tridiagonal_[i][0] + tridiagonal_[i][2] +
                        properties_->GetDeltaT() * properties_->GetAlpha4() /
                            (0.25 * properties_->GetDeltaZ(M_));
    tridiagonal_[i][1] *= -1;

    row[i] = properties_->GetHeatCapacity(i, M_, prev_layer) *
                 properties_->GetDensity(i, M_, prev_layer) *
                 prev_layer[i][M_ + 1] +
             properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(M_)) *
                 (-Lambda_z(i) +
                  properties_->GetAlpha4() * properties_->GetOutTemperature());
    row[i] *= -1;
  }

  // correct (Nx * Ib, M + 1) node
  {
    int i = properties_->GetToolStartI();
    row[i] -= properties_->GetDeltaT() * properties_->GetHeatX(i, M_);
  }

  // correct (Nx * If + 1, M + 1) node
  {
    int i = properties_->GetToolFinishI() + 1;
    row[i] -= properties_->GetDeltaT() * properties_->GetHeatX(i, M_);
  }

  // (i, M + 1), i = (Nx * Ib + 2)..(Nx * If - 1), nodes and part of (Nx * Ib +
  // 1, M + 1) and (Nx * If, M + 1) nodes
  for (int i = properties_->GetToolStartI() + 1;
       i <= properties_->GetToolFinishI(); ++i) {
    // Probably withouw / 0.5
    tridiagonal_[i][0] =
        properties_->GetDeltaT() *
        properties_->GetThermalConductivity(i - 0.5, M_, prev_layer) /
        (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaBackX(i));
    tridiagonal_[i][2] =
        properties_->GetDeltaT() *
        properties_->GetThermalConductivity(i + 0.5, M_, prev_layer) /
        (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaBackX(i + 1));

    tridiagonal_[i][1] =
        properties_->GetHeatCapacity(i, M_, prev_layer) *
            properties_->GetDensity(i, M_, prev_layer) +
        tridiagonal_[i][0] + tridiagonal_[i][2] +
        properties_->GetDeltaT() / properties_->GetToolWaveHeight() *
            properties_->GetThermalConductivity(i, M_, prev_layer) /
            (properties_->GetToolHeight() - properties_->GetToolPenetration());
    tridiagonal_[i][1] *= -1;

    row[i] = properties_->GetHeatCapacity(i, M_, prev_layer) *
                 properties_->GetDensity(i, M_, prev_layer) *
                 prev_layer[i][M_ + 1] +
             properties_->GetDeltaT() / properties_->GetToolWaveHeight() *
                 (properties_->GetThermalConductivity(i, M_, prev_layer) *
                      properties_->GetToolInitTemperature() /
                      (properties_->GetToolHeight() -
                       properties_->GetToolPenetration()) -
                  Lambda_z(i));
    row[i] *= -1;
  }

  // correct (Nx * Ib + 1, M + 1) node
  {
    int i = properties_->GetToolStartI() + 1;
    tridiagonal_[i][1] += tridiagonal_[i][0];
    tridiagonal_[i][0] *=
        0.25 * properties_->GetDeltaZ(M_) / properties_->GetToolWaveHeight();
    tridiagonal_[i][1] -= tridiagonal_[i][0];

    tridiagonal_[i][1] -=
        properties_->GetDeltaT() / properties_->GetDeltaX(i) *
        (properties_->GetToolWaveHeight() - 0.25 * properties_->GetDeltaZ(M_)) /
        properties_->GetToolWaveHeight() * properties_->GetAlpha4();

    row[i] -=
        properties_->GetDeltaT() /
        (properties_->GetDeltaX(i) * properties_->GetToolWaveHeight()) *
        (0.25 * properties_->GetDeltaZ(M_) * properties_->GetHeatOutputX() +
         (properties_->GetToolWaveHeight() -
          0.25 * properties_->GetDeltaZ(M_)) *
             properties_->GetAlpha4() * properties_->GetOutTemperature());
  }

  // correct (Nx * If, M + 1) node
  {
    int i = properties_->GetToolFinishI();
    tridiagonal_[i][1] += tridiagonal_[i][2];
    tridiagonal_[i][2] *=
        0.25 * properties_->GetDeltaZ(M_) / properties_->GetToolWaveHeight();
    tridiagonal_[i][1] -= tridiagonal_[i][2];

    tridiagonal_[i][1] -=
        properties_->GetDeltaT() / properties_->GetDeltaX(i) *
        (properties_->GetToolWaveHeight() - 0.25 * properties_->GetDeltaZ(M_)) /
        properties_->GetToolWaveHeight() * properties_->GetAlpha4();

    row[i] -=
        properties_->GetDeltaT() /
        (properties_->GetDeltaX(i) * properties_->GetToolWaveHeight()) *
        (0.25 * properties_->GetDeltaZ(M_) * properties_->GetHeatOutputX() +
         (properties_->GetToolWaveHeight() -
          0.25 * properties_->GetDeltaZ(M_)) *
             properties_->GetAlpha4() * properties_->GetOutTemperature());
  }
}
