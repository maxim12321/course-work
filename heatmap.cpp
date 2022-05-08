#include "heatmap.h"

Heatmap::Heatmap(int width, int height, QWidget* parent) : QCustomPlot(parent), width_(width), height_(height) {
    setInteractions(QCP::iNone);
    axisRect()->setupFullAxesBox(true);
    xAxis->setTickLabels(false);
    yAxis->setTickLabels(false);

    color_map_ = new QCPColorMap(xAxis, yAxis);
    color_map_->data()->setSize(width_, height_);

    QCPColorScale *colorScale = new QCPColorScale(this);
    plotLayout()->addElement(0, 1, colorScale);
    colorScale->setType(QCPAxis::atRight);
    color_map_->setColorScale(colorScale);
    colorScale->axis()->setLabel("Temperature, Celsius");

    color_map_->setGradient(QCPColorGradient::gpPolar);
    // we could have also created a QCPColorGradient instance and added own colors to
    // the gradient, see the documentation of QCPColorGradient for what's possible.

    QCPMarginGroup *marginGroup = new QCPMarginGroup(this);
    axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
    colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);

    rescaleAxes();
}

void Heatmap::SetValues(const QVector<QVector<qreal>>& data) {
    for (int x = 0; x < width_; ++x) {
        for (int y = 0; y < height_; ++y) {
            color_map_->data()->setCell(x, y, data[y][x]);
        }
    }
    color_map_->rescaleDataRange();
}
