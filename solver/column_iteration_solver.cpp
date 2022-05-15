#include "column_iteration_solver.h"
#include "tridiagonal.h"

ColumnIterationSolver::ColumnIterationSolver(PropertiesManager* properties)
    : properties_(properties) {}

Matrix ColumnIterationSolver::CalculateNextIteration(const Matrix& prev_iter,
                                                     const Matrix& semi_prev_iter) {
    N_ = prev_iter.GetRowCount() - 2;
    M_ = prev_iter.GetColumnCount() - 2;

    Matrix next(N_ + 2, M_ + 2);

    LeftColumn(prev_iter, semi_prev_iter, next[0]);
    for (int i = 1; i <= N_; ++i) {
        MiddleColumn(prev_iter, semi_prev_iter, next[i], i);
    }
    RightColumn(prev_iter, semi_prev_iter, next[N_ + 1]);

    return next;
}

void ColumnIterationSolver::LeftColumn(const Matrix& prev_iter,
                                       const Matrix& semi_prev_iter,
                                       Vector& column) {
    Matrix tridiagonal(M_ + 2, 3);

    // (0, 0) node
    long double R = properties_->GetHeatCapacity(0, 0) * properties_->GetDensity(0, 0) +
                    properties_->GetDeltaT() / 0.25 * (properties_->GetThermalConductivity(0, 0.25) / (0.5 * properties_->GetDeltaZ(1) * properties_->GetDeltaZ(1)) +
                                                       properties_->GetAlpha3() / properties_->GetDeltaZ(1) +
                                                       properties_->GetAlpha1() / properties_->GetDeltaX(1));
    long double Lambda_z = properties_->GetThermalConductivity(0, 0.25) / (0.5 * properties_->GetDeltaZ(1)) * (prev_iter[0][1] - prev_iter[0][0]);
    tridiagonal[0][0] = 1;
    tridiagonal[0][1] = -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(1)) * properties_->GetThermalConductivity(0, 0.25) / (0.5 * properties_->GetDeltaZ(1)) / R;
    column[0] = (properties_->GetHeatCapacity(0, 0) * properties_->GetDensity(0, 0) * semi_prev_iter[0][0] -
                 properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(1)) * Lambda_z* +
                 semi_prev_iter[0][0] * properties_->GetDeltaT() / 0.25 * (properties_->GetAlpha3() / properties_->GetDeltaZ(1) +
                                                                           properties_->GetAlpha1() / properties_->GetDeltaX(1))) / R;

    // (0, M + 1) node
    R = properties_->GetHeatCapacity(0, M_) * properties_->GetDensity(0, M_) +
        properties_->GetDeltaT() / 0.25 * (properties_->GetThermalConductivity(0, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_) * properties_->GetDeltaZ(M_)) +
                                           properties_->GetAlpha1() / properties_->GetDeltaX(1) +
                                           properties_->GetAlpha4() / properties_->GetDeltaZ(M_));
    Lambda_z = properties_->GetThermalConductivity(0, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_)) * (prev_iter[0][M_ + 1] - prev_iter[0][M_]);

    tridiagonal[M_ + 1][1] = 1;
    tridiagonal[M_ + 1][0] = -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(M_)) * properties_->GetThermalConductivity(0, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_)) / R;
    column[M_ + 1] =  (properties_->GetHeatCapacity(0, M_) * properties_->GetDensity(0, M_) * semi_prev_iter[0][M_ + 1] +
                       properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(M_)) * Lambda_z +
                       semi_prev_iter[0][M_ + 1] * properties_->GetDeltaT() / 0.25 * (properties_->GetAlpha4() / properties_->GetDeltaZ(M_) +
                                                                                      properties_->GetAlpha1() / properties_->GetDeltaX(1))) / R;

    // (0, i), i = 1..M, nodes
    for (int i = 1; i <= M_; ++i) {
        tridiagonal[i][0] = properties_->GetDeltaT() * properties_->GetThermalConductivity(0, i - 0.5) / (properties_->GetDeltaZ(i) * properties_->GetDeltaBackZ(i));
        tridiagonal[i][2] = properties_->GetDeltaT() * properties_->GetThermalConductivity(0, i + 0.5) / (properties_->GetDeltaZ(i) * properties_->GetDeltaBackZ(i + 1));
        tridiagonal[i][1] = properties_->GetHeatCapacity(0, i) * properties_->GetDensity(0, i) +
                            tridiagonal[i][0] + tridiagonal[i][2] +
                            properties_->GetDeltaT() * properties_->GetAlpha1() / (0.25 * properties_->GetDeltaX(1));
        tridiagonal[i][1] *= -1;

        Lambda_z = (properties_->GetThermalConductivity(0, i + 0.5) / properties_->GetDeltaBackZ(i + 1) * (prev_iter[0][i + 1] - prev_iter[0][i]) -
                    properties_->GetThermalConductivity(0, i - 0.5) / properties_->GetDeltaBackZ(i) * (prev_iter[0][i] - prev_iter[0][i - 1])) / properties_->GetDeltaZ(i);

        column[i] = properties_->GetHeatCapacity(0, i) * properties_->GetDensity(0, i) * semi_prev_iter[0][i] +
                 properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(1)) * properties_->GetAlpha1() * semi_prev_iter[0][i] -
                 properties_->GetDeltaT() * Lambda_z;
        column[i] *= -1;
    }

    assert(TridiagonalAlgorithm(tridiagonal, column));
}

void ColumnIterationSolver::MiddleColumn(const Matrix& prev_iter,
                                         const Matrix& semi_prev_iter,
                                         Vector& column, int i) {
    Matrix tridiagonal(M_ + 2, 3);

    // (i, 0) node
    long double R = properties_->GetHeatCapacity(i, 0) * properties_->GetDensity(i, 0) +
                    properties_->GetDeltaT() / 0.25 * (properties_->GetThermalConductivity(i, 0.25) / (0.5 * properties_->GetDeltaZ(1) * properties_->GetDeltaZ(1)) +
                                                       properties_->GetAlpha3() / properties_->GetDeltaZ(1));
    long double Lambda_z = properties_->GetThermalConductivity(i, 0.25) / (0.5 * properties_->GetDeltaZ(1)) * (prev_iter[i][1] - prev_iter[i][0]);

    tridiagonal[0][0] = 1;
    tridiagonal[0][1] = -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(1)) * properties_->GetThermalConductivity(i, 0.25) / (0.5 * properties_->GetDeltaZ(1)) / R;
    column[0] = (properties_->GetHeatCapacity(i, 0) * properties_->GetDensity(i, 0) * semi_prev_iter[i][0] -
                 properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(1)) * Lambda_z +
                 properties_->GetDeltaT() * properties_->GetAlpha3() * semi_prev_iter[i][0] / (0.25 * properties_->GetDeltaZ(1))) / R;

    // (N + 1, M + 1) node
    // Common part for all nodes in [1, NxIb - 1]+[NxIf + 2, Nx] and (NxIb) and (NxIf + 1) and [NxIb + 2, NxIf - 1] and (NxIb + 1) and (NxIf)
    Lambda_z = properties_->GetThermalConductivity(i, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_)) * (prev_iter[i][M_ + 1] - prev_iter[i][M_]);
    tridiagonal[M_ + 1][1] = 1;
    tridiagonal[M_ + 1][0] = -properties_->GetDeltaT() * properties_->GetThermalConductivity(i, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_));
    R = properties_->GetHeatCapacity(i, M_) * properties_->GetDensity(i, M_);

    // Correct
    if (i <= properties_->GetToolStartI() || properties_->GetToolFinishI() + 1 <= i) {
        R += properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(M_)) * (properties_->GetThermalConductivity(i, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_)) +
                                                                               properties_->GetAlpha4());

        tridiagonal[M_ + 1][0] /= 0.25 * properties_->GetDeltaZ(M_) * R;

        column[M_ + 1] =  (properties_->GetHeatCapacity(i, M_) * properties_->GetDensity(i, M_) * semi_prev_iter[i][M_ + 1] + /* (-) in doc ?*/
                           properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(M_)) * Lambda_z +
                           semi_prev_iter[i][M_ + 1] * properties_->GetDeltaT() * properties_->GetAlpha4() / (0.25 * properties_->GetDeltaZ(M_))) / R;
    }
    if (i > properties_->GetToolStartI() + 1 && i < properties_->GetToolFinishI()) {
        R += properties_->GetDeltaT() / properties_->GetToolWaveHeight() * (properties_->GetThermalConductivity(i, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_)) +
                                                                            properties_->GetThermalConductivity(i, M_) / (properties_->GetToolHeight() -
                                                                                                                          properties_->GetToolPenetration()));
        tridiagonal[M_ + 1][0] /= properties_->GetToolWaveHeight() * R;

        column[M_ + 1] = (properties_->GetHeatCapacity(i, M_) * properties_->GetDensity(i, M_) * semi_prev_iter[i][M_ + 1] +
                          properties_->GetDeltaT() / properties_->GetToolWaveHeight() * (
                                properties_->GetThermalConductivity(i, M_) * semi_prev_iter[i][M_ + 1] / (properties_->GetToolHeight() - properties_->GetToolPenetration()) +
                                Lambda_z
                )) / R;
    }
    if (i == properties_->GetToolStartI() + 1 || i == properties_->GetToolFinishI()) {
        R += properties_->GetDeltaT() / properties_->GetToolWaveHeight() * (
                    properties_->GetThermalConductivity(i, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_)) +
                    (properties_->GetToolWaveHeight() - 0.25 * properties_->GetDeltaZ(M_)) * properties_->GetAlpha4() / properties_->GetDeltaX(i) +
                    properties_->GetThermalConductivity(i, M_) / (properties_->GetToolHeight() - properties_->GetToolPenetration())
            );

        tridiagonal[M_ + 1][0] /= properties_->GetToolWaveHeight() * R;

        column[M_ + 1] = (properties_->GetHeatCapacity(i, M_) * properties_->GetDensity(i, M_) * semi_prev_iter[i][M_ + 1] +
                          properties_->GetDeltaT() / properties_->GetToolWaveHeight() * (
                                properties_->GetThermalConductivity(i, M_) * semi_prev_iter[i][M_ + 1] / (properties_->GetToolHeight() - properties_->GetToolPenetration()) +
                                (properties_->GetToolWaveHeight() - 0.25 * properties_->GetDeltaZ(M_)) * properties_->GetAlpha4() * semi_prev_iter[i][M_ + 1] / properties_->GetDeltaX(i) +
                                Lambda_z
                )) / R;
    }

    // (i, k), k = 1..M_, nodes
    for (int k = 1; k <= M_; ++k) {
        // Probably with / 0.5 (in doc)
        tridiagonal[k][0] = properties_->GetDeltaT() * properties_->GetThermalConductivity(i, k - 0.5) / (properties_->GetDeltaZ(k) * properties_->GetDeltaBackZ(k));
        tridiagonal[k][2] = properties_->GetDeltaT() * properties_->GetThermalConductivity(i, k + 0.5) / (properties_->GetDeltaZ(k) * properties_->GetDeltaBackZ(k + 1));
        tridiagonal[k][1] = properties_->GetHeatCapacity(i, k) * properties_->GetDensity(i, k) +
                            tridiagonal[k][0] + tridiagonal[k][2];
        tridiagonal[k][1] *= -1;

        Lambda_z = (properties_->GetThermalConductivity(i, k + 0.5) / properties_->GetDeltaBackZ(k + 1) * (prev_iter[i][k + 1] - prev_iter[i][k]) -
                    properties_->GetThermalConductivity(i, k - 0.5) / properties_->GetDeltaBackZ(k) * (prev_iter[i][k] - prev_iter[i][k - 1])) / properties_->GetDeltaZ(k);

        column[k] = properties_->GetHeatCapacity(i, k) * properties_->GetDensity(i, k) * semi_prev_iter[i][k] -
                    properties_->GetDeltaT() * Lambda_z;
        column[k] *= -1;
    }

    assert(TridiagonalAlgorithm(tridiagonal, column));
}

void ColumnIterationSolver::RightColumn(const Matrix& prev_iter,
                                        const Matrix& semi_prev_iter,
                                        Vector& column) {
    Matrix tridiagonal(M_ + 2, 3);

    // (N + 1, 0) node
    long double R = properties_->GetHeatCapacity(N_, 0) * properties_->GetDensity(N_, 0) +
                    properties_->GetDeltaT() / 0.25 * (properties_->GetThermalConductivity(N_, 0.25) / (0.5 * properties_->GetDeltaZ(1) * properties_->GetDeltaZ(1)) +
                                                       properties_->GetAlpha3() / properties_->GetDeltaZ(1) +
                                                       properties_->GetAlpha2() / properties_->GetDeltaX(N_));
    long double Lambda_z = properties_->GetThermalConductivity(N_, 0.25) / (0.5 * properties_->GetDeltaZ(1)) * (prev_iter[N_ + 1][1] - prev_iter[N_ + 1][0]);
    tridiagonal[0][0] = 1;
    tridiagonal[0][1] = -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(1)) * properties_->GetThermalConductivity(N_, 0.25) / (0.5 * properties_->GetDeltaZ(1)) / R;
    column[0] = (properties_->GetHeatCapacity(N_, 0) * properties_->GetDensity(N_, 0) * semi_prev_iter[N_ + 1][0] -
                 properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(1)) * Lambda_z +
                 semi_prev_iter[N_ + 1][0] * properties_->GetDeltaT() / 0.25 * (properties_->GetAlpha3() / properties_->GetDeltaZ(1) +
                                                                                properties_->GetAlpha2() / properties_->GetDeltaX(N_))) / R;

    // (N + 1, M + 1) node
    R = properties_->GetHeatCapacity(N_, M_) * properties_->GetDensity(N_, M_) +
        properties_->GetDeltaT() / 0.25 * (properties_->GetThermalConductivity(N_, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_) * properties_->GetDeltaZ(M_)) +
                                           properties_->GetAlpha2() / properties_->GetDeltaX(N_) +
                                           properties_->GetAlpha4() / properties_->GetDeltaZ(M_));
    Lambda_z = properties_->GetThermalConductivity(N_, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_)) * (prev_iter[N_ + 1][M_ + 1] - prev_iter[N_ + 1][M_]);
    tridiagonal[M_ + 1][1] = 1;
    tridiagonal[M_ + 1][0] = -properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(M_)) * properties_->GetThermalConductivity(N_, M_ - 0.25) / (0.5 * properties_->GetDeltaZ(M_)) / R;
    column[M_ + 1] =  (properties_->GetHeatCapacity(N_, M_) * properties_->GetDensity(N_, M_) * semi_prev_iter[0][M_ + 1] + /* (-) in doc ?*/
                       properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(M_)) * Lambda_z +
                       semi_prev_iter[N_ + 1][M_ + 1] * properties_->GetDeltaT() / 0.25 * (properties_->GetAlpha4() / properties_->GetDeltaZ(M_) +
                                                                                           properties_->GetAlpha2() / properties_->GetDeltaX(N_))) / R;

    // (N + 1, i), i = 1..M, nodes
    for (int i = 1; i <= M_; ++i) {
        tridiagonal[i][0] = properties_->GetDeltaT() * properties_->GetThermalConductivity(N_, i - 0.5) / (properties_->GetDeltaZ(i) * properties_->GetDeltaBackZ(i));
        tridiagonal[i][2] = properties_->GetDeltaT() * properties_->GetThermalConductivity(N_, i + 0.5) / (properties_->GetDeltaZ(i) * properties_->GetDeltaBackZ(i + 1));
        tridiagonal[i][1] = properties_->GetHeatCapacity(N_, i) * properties_->GetDensity(N_, i) +
                            tridiagonal[i][0] + tridiagonal[i][2] +
                            properties_->GetDeltaT() * properties_->GetAlpha2() / (0.25 * properties_->GetDeltaX(N_));
        tridiagonal[i][1] *= -1;

        Lambda_z = (properties_->GetThermalConductivity(N_, i + 0.5) / properties_->GetDeltaBackZ(i + 1) * (prev_iter[N_ + 1][i + 1] - prev_iter[N_ + 1][i]) -
                    properties_->GetThermalConductivity(N_, i - 0.5) / properties_->GetDeltaBackZ(i) * (prev_iter[N_ + 1][i] - prev_iter[N_ + 1][i - 1])) / properties_->GetDeltaZ(i);

        column[i] = properties_->GetHeatCapacity(N_, i) * properties_->GetDensity(N_, i) * semi_prev_iter[N_ + 1][i] +
                 properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(N_)) * properties_->GetAlpha2() * semi_prev_iter[N_ + 1][i] -
                 properties_->GetDeltaT() * Lambda_z;
        column[i] *= -1;
    }

    assert(TridiagonalAlgorithm(tridiagonal, column));
}
