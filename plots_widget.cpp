#include "plots_widget.h"

#include <QtCharts/QLineSeries>
#include <QGraphicsView>
#include <QtCharts/QChart>

PlotsWidget::PlotsWidget() : chart_(new QChart()) {
    chart_->setTitle("Line chart");
    chart_->createDefaultAxes();
    setChart(chart_);
//    chart_->axes(Qt::Horizontal).first()->setRange(0, 2);
//    chart_->axes(Qt::Vertical).first()->setRange(270, 370);
    // Add space to label to add space between labels and axis
//    QValueAxis *axisY = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).first());
//    axisY->setLabelFormat("%.1f  ");
}

void PlotsWidget::MakePlot(const QList<QPointF> &data) {
    QLineSeries* series = new QLineSeries();
    series->append(data);
    chart_ = new QChart();
    chart_->removeAllSeries();
    chart_->addSeries(series);
//    chart_->axes(Qt::Horizontal).first()->setRange(0, 2);
//    chart_->axes(Qt::Vertical).first()->setRange(270, 370);
}
