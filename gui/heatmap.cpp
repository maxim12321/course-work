#include "heatmap.h"

#include <iostream>

Heatmap::Heatmap(int width, int height, QWidget* parent) : QCustomPlot(parent), width_(width), height_(height) {
    setInteractions(QCP::iNone);
    axisRect()->setupFullAxesBox(true);
    xAxis->setTickLabels(false);
    yAxis->setTickLabels(false);

    color_map_ = new QCPColorMap(xAxis, yAxis);

    QCPRange range(20, kMaxTemperature);
    color_scale_ = new QCPColorScale(this);
    color_scale_->setDataRange(range);

    plotLayout()->addElement(0, 1, color_scale_);
    color_scale_->setType(QCPAxis::atRight);
    color_map_->setColorScale(color_scale_);
    color_scale_->axis()->setLabel("Temperature, Celsius");

    color_map_->setGradient(QCPColorGradient::gpPolar);
    // we could have also created a QCPColorGradient instance and added own colors to
    // the gradient, see the documentation of QCPColorGradient for what's possible.

    QCPMarginGroup *marginGroup = new QCPMarginGroup(this);
    axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
    color_scale_->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);


    rescaleAxes();
}

void Heatmap::SetValues(const QVector<QVector<double>>& data) {
    double min_value = data[0][0] - 273;
    double max_value = min_value;
    for (int x = 0; x < width_; ++x) {
        for (int y = 0; y < height_; ++y) {
            double value = data[y][x] - 273;
            color_map_->data()->setCell(x, y, value);
            min_value = std::min(min_value, value);
            max_value = std::max(max_value, value);
        }
    }
    QCPRange range(min_value, max_value);
    color_scale_->setDataRange(range);
    replot();
}

void Heatmap::Resize(int width, int height) {
    color_map_->data()->setSize(width, height);
    width_ = width;
    height_ = height;
}
