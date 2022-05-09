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
    // Using eps1, eps2 and sizes of objects compute grid size
    nx_ = 10;
    nz_ = 10;

    temperature_ = properties_->InitializeGrids(nx_, nz_);
}

void Solver::CalculateNextLayer() {
    Matrix semi_next = row_solver_.CalculateNextLayer(temperature_);
    Matrix next = column_solver_.CalculateNextLayer(temperature_, semi_next);

    temperature_ = next;
    on_layer_ready_(next);
}
