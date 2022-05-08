#pragma once

#include <qcustomplot/qcustomplot.h>

class Heatmap : public QCustomPlot {
public:
    Heatmap(int width, int height, QWidget* parent);
    void SetValues(const QVector<QVector<qreal>>& data);

private:
    int width_;
    int height_;
    QCPColorMap *color_map_;
};
