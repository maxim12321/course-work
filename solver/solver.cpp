#include "solver.h"

Solver::Solver(PropertiesManager* properties, Callback callback)
    : properties_(properties),
      on_layer_ready_(std::move(callback)),
      row_solver_(properties),
      column_solver_(properties) {
    Initialize();
}

void Solver::Start() {
    // Put this in a loop later
    CalculateNextLayer();
}

void Solver::Initialize() {
    nx_ = 10;
    nz_ = 10;

    current_temp_ = properties_->InitializeGrids(nx_, nz_);
    previous_temp_ = current_temp_;
}

void Solver::CalculateNextLayer() {
    int iteration = 0;
    double eps = 0;
    do {
        Matrix semi_next = row_solver_.CalculateNextIteration(current_temp_, previous_temp_);
        current_temp_ = column_solver_.CalculateNextIteration(current_temp_, semi_next);
        ++iteration;

        for (int i = 0; i < current_temp_.GetRowCount(); ++i) {
            for (int j = 0; j < current_temp_.GetColumnCount(); ++j) {
                double delta = std::fabs(current_temp_[i][j] - previous_temp_[i][j]);
                eps = std::max(eps, delta);
            }
        }
    } while (eps > properties_->GetEpsilon1() || iteration < properties_->GetMaxIterations());

    previous_temp_ = current_temp_;
    on_layer_ready_(current_temp_);
}
