#pragma once

#include <QVector>
#include <../qcustomplot/qcustomplot.h>

class Heatmap : public QCustomPlot {
public:
    Heatmap(int width, int height, QWidget* parent);
    void SetValues(const QVector<QVector<double>>& data);
    void Resize(int width, int height);

private:
    static constexpr int kMaxTemperature = 25;
    static constexpr int kMinTemperature = 20;

private:
    int width_;
    int height_;
    QCPColorMap *color_map_;
};
