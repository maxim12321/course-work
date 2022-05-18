#include "solver.h"

#include <QDebug>
#include <iostream>

Solver::Solver(int cells_x, int cells_z, PropertiesManager* properties, Callback callback)
    : properties_(properties),
      on_layer_ready_(std::move(callback)),
      nx_(cells_x),
      nz_(cells_z),
      row_solver_(properties),
      column_solver_(properties) {
    Initialize();
}

void Solver::Start() {
    // Put this in a loop later
    for (size_t i = 0; i < properties_->GetMaxIterations(); ++i) {
        CalculateNextLayer();
    }
    on_layer_ready_(current_temp_.Transposed());
}

void Solver::Initialize() {
    current_temp_ = properties_->InitializeGrids(nx_, nz_);
    previous_temp_ = current_temp_;
}

void Solver::CalculateNextLayer() {
    for (int iteration = 0; iteration < properties_->GetMaxIterations(); ++iteration) {
        Matrix semi_next = row_solver_.CalculateNextIteration(current_temp_, previous_temp_);
        Matrix next_temp = column_solver_.CalculateNextIteration(current_temp_, semi_next);
        ++iteration;

//        std::cout << "\n------------- Semi-next -------------\n";
//        std::cout << semi_next.Transposed() << std::endl;

//        std::cout << "\n------------- Next temp -------------\n";
//        std::cout << next_temp.Transposed() << std::endl << std::endl;

        if (HasConverged(current_temp_, next_temp)) {
            qDebug() << "Converged!";
            current_temp_ = next_temp;
            break;
        }
        current_temp_ = next_temp;
    }

    previous_temp_ = current_temp_;
//    on_layer_ready_(current_temp_.Transposed());
}

bool Solver::HasConverged(const Matrix& current, const Matrix& next) {
    for (int i = 0; i < current_temp_.GetRowCount(); ++i) {
        for (int j = 0; j < current_temp_.GetColumnCount(); ++j) {
            double delta = std::fabs(current[i][j] - next[i][j]);

            if (delta > properties_->GetEpsilon1() * current[i][j] + properties_->GetEpsilon2()) {
                return false;
            }
        }
    }
    return true;
}
