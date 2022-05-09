#pragma once

#include <functional>

#include "column_iteration_solver.h"
#include "properties_manager.h"
#include "row_iteration_solver.h"
#include "utils/matrix.h"

class Solver {
public:
    using Callback = std::function<void(const Matrix&)>;

    Solver(PropertiesManager* properties, Callback callback);

    // Separate thread for this later?
    void Start();

private:
    // Initialize fields with values from PropertiesManager...
    void Initialize();

    // Calculate new values for *temperature_* based on current values,
    // invoke *on_layer_ready_* callback
    void CalculateNextLayer();

private:
    PropertiesManager* properties_;
    Callback on_layer_ready_;

    Matrix temperature_;

    RowIterationSolver row_solver_;
    ColumnIterationSolver column_solver_;
};
