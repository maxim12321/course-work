#include "row_iteration_solver.h"
#include "tridiagonal.h"

RowIterationSolver::RowIterationSolver(PropertiesManager* properties)
    : properties_(properties) {}

Matrix RowIterationSolver::CalculateNextLayer(const Matrix& previous) {
    int N = previous.GetRowCount() - 2;
    int M = previous.GetColumnCount() - 2;
    Matrix next(N + 2, M + 2);

    BottomRow(previous, next[0]);
    for (int i = 1; i <= N; ++i) {
        MiddleRow(previous, next[i]);
    }
    TopRow(previous, next[N + 1]);

    return next;
}

void RowIterationSolver::BottomRow(const Matrix& previous, Vector& row) {
    Matrix tridiagonal(previous.GetColumnCount(), 3);

    auto Lambda = [&](int i) -> long double{
      return properties_->GetThermalConductivity(i, 0.25) / (0.5 * properties_->GetDeltaZ(1)) * (previous[i][1] - previous[i][0]);
    };

    // (0, 0) node
    long double R = properties_->GetHeatCapacity(0, 0) * properties_->GetDensity(0, 0) +
                    properties_->GetDeltaT() / 0.25 * (properties_->GetThermalConductivity(0.25, 0) / (0.5 * properties_->GetDeltaX(1) * properties_->GetDeltaX(1)) +
                                                       properties_->GetAlpha1() / properties_->GetDeltaX(1) +
                                                       properties_->GetAlpha3() / properties_->GetDeltaZ(1));

    tridiagonal[0][0] = 1;
    tridiagonal[1][0] = - properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(1)) * properties_->GetThermalConductivity(0.25, 0) / (0.5 * properties_->GetDeltaX(1)) / R;

    row[0] = (properties_->GetHeatCapacity(0, 0) * properties_->GetDensity(0, 0) * properties_->GetInitTemperature(0, 0) +
                properties_->GetDeltaT() / 0.25 * (properties_->GetAlpha1() * properties_->GetOutTemperature() / properties_->GetDeltaX(1) +
                                                   properties_->GetAlpha3() * properties_->GetOutTemperature() / properties_->GetDeltaZ(1) +
                                                   Lambda(0) / properties_->GetDeltaZ(1))) / R;

    int N = previous.GetColumnCount() - 2;
    // (N + 1, 0) node
    R = properties_->GetHeatCapacity(N, 0) * properties_->GetDensity(N, 0) +
        properties_->GetDeltaT() / 0.25 * (properties_->GetThermalConductivity(N - 0.25, 0) / (0.5 * properties_->GetDeltaX(N) * properties_->GetDeltaX(N)) +
                                           properties_->GetAlpha2() / properties_->GetDeltaX(N) +
                                           properties_->GetAlpha3() / properties_->GetDeltaZ(N));

    tridiagonal[N + 1][0] = 1;
    tridiagonal[N][0] = - properties_->GetDeltaT() / (0.25 * properties_->GetDeltaX(N)) * properties_->GetThermalConductivity(N - 0.25, 0) / (0.5 * properties_->GetDeltaX(N)) / R;

    row[N + 1] = (properties_->GetHeatCapacity(N, 0) * properties_->GetDensity(N, 0) * properties_->GetInitTemperature(N + 1, 0) +
                    properties_->GetDeltaT() / 0.25 * (properties_->GetAlpha2() * properties_->GetOutTemperature() / properties_->GetDeltaX(N) +
                                                       properties_->GetAlpha3() * properties_->GetOutTemperature() / properties_->GetDeltaZ(1) +
                                                       Lambda(N + 1) / properties_->GetDeltaZ(1))) / R;

    // (i, 0) node, i = 1..N
    for (int i = 1; i <= N; ++i) {
        tridiagonal[i][0] = properties_->GetDeltaT() * properties_->GetThermalConductivity(i + 0.5, 0) / (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaX(i + 1));
        tridiagonal[i][2] = properties_->GetDeltaT() * properties_->GetThermalConductivity(i - 0.5, 0) / (0.5 * properties_->GetDeltaX(i) * properties_->GetDeltaX(i));
        tridiagonal[i][1] = -properties_->GetHeatCapacity(i, 0) * properties_->GetDensity(i, 0) -
                            tridiagonal[i][0] - tridiagonal[i][2] -
                            properties_->GetDeltaT() * properties_->GetAlpha3() / (0.25 * properties_->GetDeltaZ(1));

        row[i] = properties_->GetHeatCapacity(i, 0) * properties_->GetDensity(i, 0) * properties_->GetInitTemperature(i, 0) +
                 properties_->GetDeltaT() / (0.25 * properties_->GetDeltaZ(1)) * (properties_->GetAlpha3() * properties_->GetOutTemperature() + Lambda(i));
    }

    assert(TridiagonalAlgorithm(tridiagonal, row));
}

void RowIterationSolver::TopRow(const Matrix &previous, Vector &row) {
    Matrix tridiagonal(previous.GetColumnCount(), 3);
    // initializing
    assert(TridiagonalAlgorithm(tridiagonal, row));
}
void RowIterationSolver::MiddleRow(const Matrix &previous, Vector &row) {
    Matrix tridiagonal(previous.GetColumnCount(), 3);
    // initializing
    assert(TridiagonalAlgorithm(tridiagonal, row));
}
