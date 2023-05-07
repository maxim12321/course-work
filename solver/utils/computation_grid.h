#pragma once

#include "matrix.h"

class ComputationGrid : public Matrix {
public:
    ComputationGrid();
    ComputationGrid(int n, int m, int j_plate_finish, int i_tool_start, int i_tool_finish);

    Vector GetRow(int i);
    Vector GetColumn(int i);

private:
    int j_plate_finish_;
    int i_tool_start_;
    int i_tool_finish_;
}; 