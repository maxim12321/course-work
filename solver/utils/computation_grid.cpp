#include "computation_grid.h"

ComputationGrid::ComputationGrid() {}

ComputationGrid::ComputationGrid(int n, int m, int j_plate_finish, int i_tool_start, int i_tool_finish)
    : Matrix(n, m), j_plate_finish_(j_plate_finish), i_tool_start_(i_tool_start), i_tool_finish_(i_tool_finish) {
}

Vector ComputationGrid::GetRow(int i) {
    int size = i < i_tool_start_ || i >= i_tool_finish_ ? j_plate_finish_ : columns_;
    return GetVector(i, 0, 1, size);
}

Vector ComputationGrid::GetColumn(int i) {
    int start = i < j_plate_finish_ ? 0 : i_tool_start_;
    int size = i < j_plate_finish_ ? rows_ : i_tool_finish_ - i_tool_start_;
    return GetVector(start, i, columns_, size);
}
