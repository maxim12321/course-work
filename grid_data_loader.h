#pragma once

#include <QString>
#include <QVector>
#include <QPointF>

class GridDataLoader {
public:
    GridDataLoader(int width, int height);

    QVector<QVector<qreal>> LoadGridData() const;
    QVector<QPointF> LoadPlotData() const;

private:
    const QString kGridFileName = "./solver_input/Te_C.dat";
    const QString kPlotFileName = "./solver_input/TeS1_zI.dat";

private:
    int width_;
    int height_;
};
