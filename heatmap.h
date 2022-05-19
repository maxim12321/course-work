#pragma once

#include <qcustomplot/qcustomplot.h>

#include "solver/utils/matrix.h"

class Heatmap : public QCustomPlot {
public:
    Heatmap(int width, int height, QWidget* parent);
    void SetValues(const Matrix& data);

private:
    static constexpr int kMaxTemperature = 120;

private:
    int width_;
    int height_;
    QCPColorMap *color_map_;
};
