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
    temperature_ = Matrix(10, 10);
}

void Solver::CalculateNextLayer() {
    Matrix semi_next = row_solver_.CalculateNextLayer(temperature_);
    Matrix next = column_solver_.CalculateNextLayer(temperature_, semi_next);

    temperature_ = next;
    on_layer_ready_(next);
}
