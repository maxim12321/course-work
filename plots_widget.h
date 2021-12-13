#pragma once

#include <QtCharts/QChartView>
#include <QtCharts/QChart>

class PlotsWidget : public QChartView {
public:
    PlotsWidget();
    void MakePlot(const QList<QPointF>& data);

private:
    QChart* chart_;
};
