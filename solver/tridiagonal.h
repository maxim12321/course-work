#pragma once

#include "utils/matrix.h"
#include "utils/vector.h"

bool TridiagonalAlgorithm(Matrix &matrix, Vector &right);

bool Gauss(Matrix &matrix, Vector &right);

bool ReverseGauss(const Matrix &matrix, Vector &right);
